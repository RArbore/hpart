use bitvec::prelude::*;

pub(crate) type Index = u32;

/// The data structure internal to this hypergraph partitioner is a bipartite
/// graph between pins and nets.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Bipartite {
    // V, E, and A arrays as described in section 4.5 of Schlag '2015.
    pub(crate) v: Vec<(Index, Index)>,
    pub(crate) e: Vec<(Index, Index)>,
    pub(crate) a: Vec<Index>,

    // Use an explicit bit vector to represent whether vertices are enabled or
    // disabled.
    pub(crate) v_enabled: BitVec,

    // Capacities and weights of pins and nets.
    pub(crate) c: Vec<f32>,
    pub(crate) w: Vec<f32>,
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

impl Bipartite {
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
                    self.v[u as usize].0 = self.a.len() as Index - self.v[u as usize].1;
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
    pub(crate) fn uncontract(&mut self, m: Memento) {
        self.v_enabled.set(m.v as usize, true);
        let mut b = bitvec![usize, Lsb0; 0; self.e.len()];
        for e in self.incident_nets(m.v) {
            b.set(e as usize, true);
        }
        for i in m.u_idx..m.u_idx + m.u_len {
            b.set(self.a[i as usize] as usize, false);
        }

        if self.v[m.u as usize].1 > m.u_len {
            let inc_u: Vec<_> = self.incident_nets(m.u).collect();
            for e in inc_u {
                if b[e as usize] {
                    let e_idx = self.e[e as usize].0;
                    let e_len = self.e[e as usize].1;
                    for p_idx in e_idx..e_idx + e_len {
                        if self.a[p_idx as usize] == m.u {
                            self.a[p_idx as usize] = m.v;
                            break;
                        }
                    }
                }
            }
        }

        self.v[m.u as usize].0 = m.u_idx;
        self.v[m.u as usize].1 = m.u_len;
        self.c[m.u as usize] -= self.c[m.v as usize];
        let inc_v: Vec<_> = self.incident_nets(m.v).collect();
        for e in inc_v {
            if !b[e as usize] {
                self.e[e as usize].1 += 1;
            }
        }
    }

    pub(crate) fn incident_nets(&self, v: Index) -> impl Iterator<Item = Index> + '_ {
        let v_idx = self.v[v as usize].0 as usize;
        let v_len = self.v[v as usize].1 as usize;
        self.a[v_idx..v_idx + v_len].iter().map(|x| *x)
    }

    pub(crate) fn incident_pins(&self, v: Index) -> impl Iterator<Item = Index> + '_ {
        self.incident_nets(v)
            .map(|e| {
                let e_idx = self.e[e as usize].0 as usize;
                let e_len = self.e[e as usize].1 as usize;
                self.a[e_idx..e_idx + e_len].iter().map(|x| *x)
            })
            .flatten()
    }
}

/// A bipartition is just the vertices in each half. Technically, only one of
/// these is necessary to derive the other, but it's convenient to have both.
pub(crate) type Bipartition = (Vec<Index>, Vec<Index>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_contract() {
        // Construct the hypergraph shown in Figure 2 in Schlag '2015.
        let original = Bipartite {
            v: vec![(0, 1), (1, 2), (3, 1), (4, 1), (5, 1), (6, 1)],
            e: vec![(7, 4), (11, 3)],
            a: vec![0, 0, 1, 0, 0, 1, 1, 0, 1, 2, 3, 1, 4, 5],
            v_enabled: bitvec![usize, Lsb0; 1; 6],
            c: vec![1.0; 6],
            w: vec![1.0; 2],
        };

        let mut contract = original.clone();
        let m = contract.contract(0, 1);

        let correct = Bipartite {
            v: vec![(14, 2), (1, 2), (3, 1), (4, 1), (5, 1), (6, 1)],
            e: vec![(7, 3), (11, 3)],
            a: vec![0, 0, 1, 0, 0, 1, 1, 0, 3, 2, 1, 5, 4, 0, 0, 1],
            v_enabled: bitvec![1, 0, 1, 1, 1, 1],
            c: vec![2.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            w: vec![1.0; 2],
        };
        assert_eq!(contract, correct);

        contract.uncontract(m);

        let correct = Bipartite {
            v: vec![(0, 1), (1, 2), (3, 1), (4, 1), (5, 1), (6, 1)],
            e: vec![(7, 4), (11, 3)],
            a: vec![0, 0, 1, 0, 0, 1, 1, 0, 3, 2, 1, 5, 4, 1, 0, 1],
            v_enabled: bitvec![usize, Lsb0; 1; 6],
            c: vec![1.0; 6],
            w: vec![1.0; 2],
        };
        assert_eq!(contract, correct);
    }
}
