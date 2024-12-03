use std::collections::BinaryHeap;

use bitvec::prelude::*;
use ordered_float::OrderedFloat;

use crate::bipartite::*;

pub(crate) fn coarsen(h: &mut Bipartite) -> Vec<Memento> {
    let c_max = S * h.total_capacity() / T as f32;

    let mut pq = BinaryHeap::new();
    for u in 0..h.v.len() as Index {
        if h.c[u as usize] > c_max {
            continue;
        }

        let v = h
            .incident_pins(u)
            .max_by_key(|v| OrderedFloat(rate(h, u, *v)))
            .unwrap();
        pq.push((OrderedFloat(rate(h, u, v)), (u, v)));
    }
    let mut removed = bitvec![usize, Lsb0; 0; h.v.len()];
    let mut invalid = bitvec![usize, Lsb0; 0; h.v.len()];

    let mut mementos = vec![];
    while h.v.len() >= T {
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
    let inv_c = 1.0 / (h.c[u as usize] * h.c[v as usize]);
    let mut heavy_edge = 0.0;
    for e in 0..h.e.len() {
        let e_idx = h.e[e].0 as usize;
        let e_len = h.e[e].1 as usize;
        let e_adj = &h.a[e_idx..e_idx + e_len];
        if e_adj.contains(&u) && e_adj.contains(&v) {
            heavy_edge += h.w[e] / (e_adj.len() - 1) as f32;
        }
    }
    inv_c * heavy_edge
}
