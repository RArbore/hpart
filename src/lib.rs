mod bipartite;
mod coarsen;
mod initial;
mod uncoarsen;

use rand::prelude::*;

use bipartite::{Bipartite, Bipartition, Index};
use coarsen::coarsen;
use initial::initial_partitioning;
use uncoarsen::uncoarsen;

pub fn bipartition(
    capacities: &[f32],
    weights: &[f32],
    nets: &[&[Index]],
    epsilon: f32,
) -> Bipartition {
    let mut h = Bipartite::new(capacities, weights, nets);
    let size_constraint = h.size_constraint(epsilon);
    let mementos = coarsen(&mut h);
    let mut bipart = initial_partitioning(&h, epsilon);
    uncoarsen(&mut h, size_constraint, mementos, &mut bipart);
    bipart
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    #[test]
    fn bipartition_random_hypergraph() {
        let num_v = 100;
        let num_e = 300;

        let capacities: Vec<_> = (0..num_v).map(|_| random::<f32>()).collect();
        let weights: Vec<_> = (0..num_e).map(|_| random::<f32>()).collect();
        let nets: Vec<Vec<_>> = (0..num_e)
            .map(|_| {
                let pins: BTreeSet<_> = (0..num_v).map(|_| random::<Index>() % num_v).collect();
                pins.into_iter().collect()
            })
            .collect();
        let nets_ref: Vec<&[_]> = nets.iter().map(|vec| &**vec).collect();
        bipartition(&capacities, &weights, &nets_ref, 0.03);
    }
}
