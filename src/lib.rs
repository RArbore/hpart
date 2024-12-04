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
) -> (Bipartition, (f32, f32)) {
    let mut h = Bipartite::new(capacities, weights, nets);
    let size_constraint = h.size_constraint(epsilon);
    let mementos = coarsen(&mut h);
    let mut bipart = initial_partitioning(&h, epsilon);
    uncoarsen(&mut h, size_constraint, mementos, &mut bipart);
    let eval = h.evaluate_bipartition(&bipart);
    (bipart, eval)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use ordered_float::OrderedFloat;
    use rayon::prelude::*;

    use super::*;

    #[test]
    fn bipartition_random_hypergraph() {
        let num_v = 200;
        let num_e = 80;
        let max_net_size = 4;

        let capacities: Vec<_> = (0..num_v).map(|_| random::<f32>()).collect();
        let weights: Vec<_> = (0..num_e).map(|_| random::<f32>()).collect();
        let nets: Vec<Vec<_>> = (0..num_e)
            .map(|_| {
                let pins: BTreeSet<_> = (0..max_net_size)
                    .map(|_| random::<Index>() % num_v)
                    .collect();
                pins.into_iter().collect()
            })
            .collect();
        let nets_ref: Vec<&[_]> = nets.iter().map(|vec| &**vec).collect();
        let best_eval = (0..32)
            .into_par_iter()
            .map(|_| bipartition(&capacities, &weights, &nets_ref, 0.03).1)
            .min_by_key(|(_, cost)| OrderedFloat(*cost));
        println!(
            "{:?} {:?} {:?}",
            best_eval,
            capacities.iter().sum::<f32>(),
            weights.iter().sum::<f32>()
        );
    }
}
