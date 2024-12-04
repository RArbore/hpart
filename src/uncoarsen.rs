use std::collections::BinaryHeap;

use bitvec::prelude::*;
use ordered_float::OrderedFloat;

use crate::bipartite::*;

pub(crate) fn uncoarsen(
    h: &mut Bipartite,
    size_constraint: f32,
    mementos: Vec<Memento>,
    partition: &mut Bipartition,
) {
    //let mut b1_pq = BinaryHeap::new();
    //let mut b2_pq = BinaryHeap::new();
    //let mut b1_total = 0.0;
    //let mut b2_total = 0.0;
    //let mut active = bitvec![usize, Lsb0; 0; h.pin_index_space_size()];
    //let mut marked = bitvec![usize, Lsb0; 0; h.pin_index_space_size()];
    //for v in h.pins() {
    //    if partition[v as usize] {
    //        b1_total += h.capacity(v);
    //    } else {
    //        b2_total += h.capacity(v);
    //    }
    //}
    //
    //let activate = |v: Index,
    //                h: &Bipartite,
    //                partition: &Bipartition,
    //                b1_pq: &mut BinaryHeap<(OrderedFloat<f32>, Index)>,
    //                b2_pq: &mut BinaryHeap<(OrderedFloat<f32>, Index)>,
    //                active: &mut BitVec| {
    //    let b_v = partition[v as usize];
    //    let g = OrderedFloat(gain(h, v, !b_v, partition));
    //    if b_v {
    //        b2_pq.push((g, v));
    //    } else {
    //        b1_pq.push((g, v));
    //    }
    //    active.set(v as usize, true);
    //};

    for memento in mementos.into_iter().rev() {
        let u = memento.u;
        let v = memento.v;
        h.uncontract(memento);
        partition[v as usize] = partition[u as usize];

        //let border_u = h
        //    .incident_pins(u)
        //    .any(|p| partition[p as usize] != partition[u as usize]);
        //if border_u {
        //    activate(u, h, partition, &mut b1_pq, &mut b2_pq, &mut active);
        //}
        //let border_v = h
        //    .incident_pins(v)
        //    .any(|p| partition[p as usize] != partition[v as usize]);
        //if border_v {
        //    activate(v, h, partition, &mut b1_pq, &mut b2_pq, &mut active);
        //}
        //if !border_u && !border_v {
        //    continue;
        //}
    }
}

fn gain(h: &Bipartite, v: Index, i: bool, partition: &Bipartition) -> f32 {
    let b_v = partition[v as usize];
    let upside: f32 = h
        .incident_nets(v)
        .filter(|e| {
            h.pins_in_net(*e)
                .filter(|n| partition[*n as usize] == i)
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