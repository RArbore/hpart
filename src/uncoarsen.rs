use std::collections::BinaryHeap;

use bitvec::prelude::*;
use ordered_float::OrderedFloat;

use crate::bipartite::*;

const C: usize = 20;

pub(crate) fn uncoarsen(
    h: &mut Bipartite,
    size_constraint: f32,
    mementos: Vec<Memento>,
    partition: &mut Bipartition,
) {
    let mut gain_pq = BinaryHeap::new();
    let mut gain_vec = vec![0.0; h.pin_index_space_size()];

    for memento in mementos.into_iter().rev() {
        let u = memento.u;
        let v = memento.v;
        h.uncontract(memento);
        partition[v as usize] = partition[u as usize];

        let border_u = h
            .incident_pins(u)
            .any(|p| partition[p as usize] != partition[u as usize]);
        if border_u {
            let gain = gain(h, u, partition);
            gain_pq.push((OrderedFloat(gain), u));
            gain_vec[u as usize] = gain;
        }
        let border_v = h
            .incident_pins(v)
            .any(|p| partition[p as usize] != partition[v as usize]);
        if border_v {
            let gain = gain(h, v, partition);
            gain_pq.push((OrderedFloat(gain), v));
            gain_vec[v as usize] = gain;
        }
        if !border_u && !border_v {
            continue;
        }

        let mut steps = vec![];
        let mut best_step = 0;
        let mut best_gain = 0.0;
        let mut current_gain = 0.0;
        let mut non_increase_count = 0;
        while let Some((OrderedFloat(g), v)) = gain_pq.pop() {
            if g != gain_vec[v as usize] {
                gain_pq.push((OrderedFloat(gain_vec[v as usize]), v));
                continue;
            }
            if g <= 0.0 {
                non_increase_count += 1;
            } else {
                non_increase_count = 0;
            }
            if non_increase_count >= C {
                break;
            }

            current_gain += g;
            partition[v as usize] = !partition[v as usize];
            steps.push(v);
            if current_gain >= best_gain {
                best_step = steps.len();
                best_gain = current_gain;
            }
            for p in h.incident_pins(v) {
                gain_vec[p as usize] = gain(h, p, partition);
            }
        }
        for idx in (best_step..steps.len()).rev() {
            let step = steps[idx];
            partition[step as usize] = !partition[step as usize];
        }
    }
}

fn gain(h: &Bipartite, v: Index, partition: &Bipartition) -> f32 {
    let b_v = partition[v as usize];
    let upside: f32 = h
        .incident_nets(v)
        .filter(|e| {
            h.pins_in_net(*e)
                .filter(|n| partition[*n as usize] == !b_v)
                .count()
                == h.pins_in_net(*e).len() - 1
        })
        .map(|e| h.weight(e))
        .sum();
    let downside: f32 = h
        .incident_nets(v)
        .filter(|e| {
            h.pins_in_net(*e)
                .filter(|n| partition[*n as usize] == b_v)
                .count()
                == h.pins_in_net(*e).len()
        })
        .map(|e| h.weight(e))
        .sum();
    upside - downside
}
