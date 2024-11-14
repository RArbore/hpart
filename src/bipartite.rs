use bitvec::prelude::*;

type Index = u32;

/// The data structure internal to this hypergraph partitioner is a bipartite
/// graph between pins and nets.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Bipartite<V, E> {
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

impl<V, E> Bipartite<V, E> {
    /// Implements Algorithm 2: Contract from Schlag '2015.
    pub(crate) fn contract(&mut self, u: Index, v: Index) -> Memento {
        todo!()
    }

    /// Implements Algorithm 3: Uncontract from Schlag '2015.
    pub(crate) fn uncontract(&mut self, memento: Memento) {
        todo!()
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
