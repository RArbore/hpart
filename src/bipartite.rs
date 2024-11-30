use std::ops::{Add, AddAssign};

use bitvec::prelude::*;

type Index = u32;

/// The data structure internal to this hypergraph partitioner is a bipartite
/// graph between pins and nets.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Bipartite<V: Add + AddAssign + Copy + Sized, E: Add + AddAssign + Copy + Sized> {
    // V, E, and A arrays as described in section 4.5 of Schlag '2015.
    pub(crate) v: Vec<(Index, Index)>,
    pub(crate) e: Vec<(Index, Index)>,
    pub(crate) a: Vec<Index>,

    // Use an explicit bit vector to represent whether vertices are enabled or
    // disabled.
    pub(crate) v_enabled: BitVec,

    // Capacities and weights of pins and nets.
    pub(crate) c: Vec<V>,
    pub(crate) w: Vec<E>,
}

/// Contractions produce mementos that can be applied in reverse to perform
/// uncontractions.
pub(crate) struct Memento {
    // Stores the pair of vertices contracted, as well as the slice into A of u
    // prior to contraction.
    pub(crate) u: Index,
    pub(crate) v: Index,
    pub(crate) u_idx: Index,
    pub(crate) u_len: Index,
}

impl<V: Add + AddAssign + Copy + Sized, E: Add + AddAssign + Copy + Sized> Bipartite<V, E> {
    /// Implements Algorithm 2: Contract from Schlag '2015.
    pub(crate) fn contract(&mut self, u: Index, v: Index) -> Memento {
        assert_ne!(u, v);
        let u_idx = self.v[u as usize].0;
        let u_len = self.v[u as usize].1;
        let memento = Memento { u, v, u_idx, u_len };

        let c_v = self.c[v as usize];
        self.c[u as usize] += c_v;
        let mut copy = true;
        let incident: Vec<_> = self.incident_nets(v).collect();
        for e in incident {
            let e_idx = self.e[e as usize].0;
            let e_len = self.e[e as usize].1;
            let l = e_idx + e_len - 1;
            let mut tau = l;
            for i in e_idx..=l {
                if self.a[i as usize] == v {
                    let a_i = self.a[i as usize];
                    self.a[i as usize] = self.a[l as usize];
                    self.a[l as usize] = a_i;
                }

                if self.a[i as usize] == u {
                    tau = i;
                }
            }

            if tau == l {
                self.a[l as usize] = u;
                if copy {
                    let incident_u: Vec<_> = self.incident_nets(u).collect();
                    self.a.extend(incident_u);
                    self.v[u as usize].0 = self.a.len() as u32 - self.v[u as usize].1;
                    copy = false;
                }
                self.a.push(e);
                self.v[u as usize].1 += 1;
            } else {
                self.e[e as usize].1 -= 1;
            }
        }

        self.v_enabled.set(v as usize, false);
        memento
    }

    /// Implements Algorithm 3: Uncontract from Schlag '2015.
    pub(crate) fn uncontract(&mut self, memento: Memento) {
        todo!()
    }

    fn incident_nets(&self, v: Index) -> impl Iterator<Item = Index> + '_ {
        let v_idx = self.v[v as usize].0 as usize;
        let v_len = self.v[v as usize].1 as usize;
        self.a[v_idx..v_idx + v_len].iter().map(|x| *x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_contract() {
        // Construct the hypergraph shown in Figure 2 in Schlag '2015.
        let mut contract = Bipartite {
            v: vec![(0, 1), (1, 2), (3, 1), (4, 1), (5, 1), (6, 1)],
            e: vec![(7, 4), (11, 3)],
            a: vec![0, 0, 1, 0, 0, 1, 1, 0, 1, 2, 3, 1, 4, 5],
            v_enabled: bitvec![usize, Lsb0; 1; 6],
            c: vec![1; 6],
            w: vec![1; 2],
        };

        contract.contract(0, 1);

        let correct = Bipartite {
            v: vec![(14, 2), (1, 2), (3, 1), (4, 1), (5, 1), (6, 1)],
            e: vec![(7, 3), (11, 3)],
            a: vec![0, 0, 1, 0, 0, 1, 1, 0, 3, 2, 1, 5, 4, 0, 0, 1],
            v_enabled: bitvec![1, 0, 1, 1, 1, 1],
            c: vec![2, 1, 1, 1, 1, 1],
            w: vec![1; 2],
        };

        assert_eq!(contract, correct);
    }
}
