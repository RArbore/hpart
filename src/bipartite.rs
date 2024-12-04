use std::iter::zip;

use bitvec::prelude::*;

pub type Index = u32;

/// The data structure internal to this hypergraph partitioner is a bipartite
/// graph between pins and nets.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Bipartite {
    // V, E, and A arrays as described in section 4.5 of Schlag '2015.
    v: Vec<(Index, Index)>,
    e: Vec<(Index, Index)>,
    a: Vec<Index>,

    // Use an explicit bit vector to represent whether vertices are enabled or
    // disabled.
    v_enabled: BitVec,
    num_disabled: usize,

    // Capacities and weights of pins and nets.
    c: Vec<f32>,
    w: Vec<f32>,
}

/// Contractions produce mementos that can be applied in reverse to perform
/// uncontractions.
#[derive(Clone, Copy)]
pub(crate) struct Memento {
    // Stores the pair of vertices contracted, as well as the slice into A of u
    // prior to contraction.
    pub(crate) u: Index,
    pub(crate) v: Index,
    u_idx: Index,
    u_len: Index,
}

impl Bipartite {
    pub(crate) fn new(capacities: &[f32], weights: &[f32], nets: &[&[Index]]) -> Self {
        let num_v = capacities.len();
        let num_e = weights.len();
        assert_eq!(num_e, nets.len());
        let mut bipartite = Bipartite {
            v: vec![],
            e: vec![],
            a: vec![],
            v_enabled: bitvec![usize, Lsb0; 1; num_v],
            num_disabled: 0,
            c: Vec::from(capacities),
            w: Vec::from(weights),
        };

        for net in nets {
            bipartite
                .e
                .push((bipartite.a.len() as Index, net.len() as Index));
            bipartite.a.extend(*net);
        }
        for v_idx in 0..num_v {
            let old_a = bipartite.a.len();
            for (e_idx, net) in nets.into_iter().enumerate() {
                if net.into_iter().any(|p| *p == v_idx as Index) {
                    bipartite.a.push(e_idx as Index);
                }
            }
            let num_nets = bipartite.a.len() - old_a;
            bipartite.v.push((old_a as Index, num_nets as Index));
        }

        bipartite
    }

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
        self.num_disabled += 1;
        memento
    }

    /// Implements Algorithm 3: Uncontract from Schlag '2015.
    pub(crate) fn uncontract(&mut self, m: Memento) {
        self.v_enabled.set(m.v as usize, true);
        self.num_disabled -= 1;
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

    pub(crate) fn incident_nets(
        &self,
        v: Index,
    ) -> impl ExactSizeIterator<Item = Index> + Clone + '_ {
        let v_idx = self.v[v as usize].0 as usize;
        let v_len = self.v[v as usize].1 as usize;
        self.a[v_idx..v_idx + v_len].iter().map(|x| *x)
    }

    pub(crate) fn pins_in_net(
        &self,
        e: Index,
    ) -> impl ExactSizeIterator<Item = Index> + Clone + '_ {
        let e_idx = self.e[e as usize].0 as usize;
        let e_len = self.e[e as usize].1 as usize;
        self.a[e_idx..e_idx + e_len].iter().map(|x| *x)
    }

    pub(crate) fn incident_pins(&self, v: Index) -> impl Iterator<Item = Index> + Clone + '_ {
        self.incident_nets(v)
            .map(|e| self.pins_in_net(e))
            .flatten()
            .filter(move |u| *u != v)
    }

    pub(crate) fn total_capacity(&self) -> f32 {
        self.c
            .iter()
            .enumerate()
            .filter(|(idx, _)| self.v_enabled[*idx])
            .map(|(_, c)| c)
            .sum()
    }

    pub(crate) fn pins(&self) -> impl Iterator<Item = Index> + Clone + '_ {
        (0..self.v.len())
            .filter(|idx| self.v_enabled[*idx])
            .map(|idx| idx as Index)
    }

    pub(crate) fn nets(&self) -> impl Iterator<Item = Index> + Clone {
        (0..self.e.len()).map(|idx| idx as Index)
    }

    pub(crate) fn enabled(&self, v: Index) -> bool {
        self.v_enabled[v as usize]
    }

    pub(crate) fn capacity(&self, v: Index) -> f32 {
        assert!(self.v_enabled[v as usize]);
        self.c[v as usize]
    }

    pub(crate) fn weight(&self, e: Index) -> f32 {
        self.w[e as usize]
    }

    pub(crate) fn num_pins(&self) -> usize {
        self.v.len() - self.num_disabled
    }

    pub(crate) fn pin_index_space_size(&self) -> usize {
        self.v.len()
    }

    pub(crate) fn num_nets(&self) -> usize {
        self.e.len()
    }

    /// Evaluates a bipartition. Returns the bipartition's imbalance (capacity of
    /// the largest cluster) and cost (weight of the cut edges).
    pub(crate) fn evaluate_bipartition(&self, p: &Bipartition) -> (f32, f32) {
        let mut cut_weight = 0.0;
        for e in self.nets() {
            let clusters = self.pins_in_net(e).map(|v| p[v as usize]);
            if zip(clusters.clone(), clusters.skip(1)).any(|(l1, l2)| l1 != l2) {
                cut_weight += self.w[e as usize];
            }
        }

        let mut capacity1 = 0.0;
        let mut capacity2 = 0.0;
        for v in self.pins() {
            let c = self.c[v as usize];
            if p[v as usize] {
                capacity1 += c;
            } else {
                capacity2 += c;
            }
        }
        let imbalance = capacity1.max(capacity2);

        (imbalance, cut_weight)
    }

    /// Calculate the maximum size for a bipartition, given epsilon.
    pub(crate) fn size_constraint(&self, epsilon: f32) -> f32 {
        (1.0 + epsilon) * (self.total_capacity() as f32 / 2.0)
    }
}

/// A bipartition is an assignment of a bool to each vertex.
pub(crate) type Bipartition = Vec<bool>;

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
            num_disabled: 0,
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
            num_disabled: 1,
            c: vec![2.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            w: vec![1.0; 2],
        };
        assert_eq!(contract, correct);

        contract.uncontract(m);

        let correct = Bipartite {
            v: vec![(0, 1), (1, 2), (3, 1), (4, 1), (5, 1), (6, 1)],
            e: vec![(7, 4), (11, 3)],
            a: vec![0, 0, 1, 0, 0, 1, 1, 0, 3, 2, 1, 5, 4, 1, 0, 1],
            num_disabled: 0,
            v_enabled: bitvec![usize, Lsb0; 1; 6],
            c: vec![1.0; 6],
            w: vec![1.0; 2],
        };
        assert_eq!(contract, correct);
    }
}
