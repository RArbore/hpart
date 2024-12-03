use std::collections::BinaryHeap;

use bitvec::prelude::*;
use ordered_float::OrderedFloat;

use crate::bipartite::*;

pub(crate) fn coarsen(h: &mut Bipartite) -> Vec<Memento> {
    let c_max = S * h.total_capacity() / T as f32;

    let mut pq = BinaryHeap::new();
    for u in h.pins() {
        if h.capacity(u) > c_max {
            continue;
        }

        let v = h
            .incident_pins(u)
            .max_by_key(|v| OrderedFloat(rate(h, u, *v)))
            .unwrap();
        pq.push((OrderedFloat(rate(h, u, v)), (u, v)));
    }
    let mut removed = bitvec![usize, Lsb0; 0; h.pin_index_space_size()];
    let mut invalid = bitvec![usize, Lsb0; 0; h.pin_index_space_size()];

    let mut mementos = vec![];
    while h.num_pins() >= T {
        let Some((_, (u, v))) = pq.pop() else {
            continue;
        };

        if removed[u as usize] {
            continue;
        } else if invalid[u as usize] {
            let v = h
                .incident_pins(u)
                .max_by_key(|v| OrderedFloat(rate(h, u, *v)))
                .unwrap();
            pq.push((OrderedFloat(rate(h, u, v)), (u, v)));
            invalid.set(u as usize, false);
            continue;
        }

        mementos.push(h.contract(u, v));
        removed.set(v as usize, true);

        for n in h.incident_pins(u) {
            invalid.set(n as usize, true);
        }
    }

    mementos
}

/// Constants from Section 5 from Schlag '2015.
const T: usize = 320;
const S: f32 = 3.25;

fn rate(h: &Bipartite, u: Index, v: Index) -> f32 {
    let inv_c = 1.0 / (h.capacity(u) * h.capacity(v));
    let mut heavy_edge = 0.0;
    for e in h.nets() {
        if h.pins_in_net(e).any(|p| p == u) && h.pins_in_net(e).any(|p| p == u) {
            heavy_edge += h.weight(e) / (h.pins_in_net(e).len() - 1) as f32;
        }
    }
    inv_c * heavy_edge
}
