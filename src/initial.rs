use std::collections::VecDeque;
use std::mem::replace;

use bitvec::prelude::*;
use ordered_float::OrderedFloat;
use rand::prelude::*;

use crate::bipartite::*;

const NUM_ITERS: usize = 20;

/// Compute an initial partitioning on a coarsened hypergraph. Runs a portfolio
/// of algorithms and returns the best found bipartition.
pub(crate) fn initial_partitioning(h: &Bipartite, epsilon: f32) -> Bipartition {
    let random_parts = (0..NUM_ITERS).map(|_| random_partitioning(h));
    let bfs_half_parts = (0..NUM_ITERS).map(|_| bfs_half(h));
    let sclap_parts = (0..NUM_ITERS).map(|_| size_constrained_label_propagation(h, epsilon));
    let all_partitions: Vec<_> = random_parts
        .chain(bfs_half_parts)
        .chain(sclap_parts)
        .collect();

    let size_constraint = h.size_constraint(epsilon);
    all_partitions
        .into_iter()
        .min_by_key(|part| {
            let (imbalance, weight) = h.evaluate_bipartition(part);
            (
                OrderedFloat(imbalance.max(size_constraint)),
                OrderedFloat(weight),
            )
        })
        .unwrap()
}

fn random_partitioning(h: &Bipartite) -> Bipartition {
    (0..h.pin_index_space_size()).map(|_| random()).collect()
}

fn bfs_half(h: &Bipartite) -> Bipartition {
    let start_v = random_vertex(h);
    let mut queue = VecDeque::new();
    queue.push_back(start_v);
    let mut visited = bitvec![usize, Lsb0; 0; h.pin_index_space_size()];
    visited.set(start_v as usize, true);

    let mut bipartition = vec![false; h.pin_index_space_size()];
    let mut num_visited = 0;
    while let Some(pop) = queue.pop_front() {
        for neighbor in h.incident_pins(pop) {
            if !visited[neighbor as usize] {
                visited.set(neighbor as usize, true);
                queue.push_back(neighbor);
            }
        }

        if num_visited < h.num_pins() / 2 {
            bipartition[pop as usize] = true;
            num_visited += 1;
        } else {
            break;
        }
    }

    bipartition
}

const TAU: usize = 5;

fn size_constrained_label_propagation(h: &Bipartite, epsilon: f32) -> Bipartition {
    let (v1, v2) = pseudo_peripheral_vertices(h);
    let mut labels = vec![None; h.pin_index_space_size()];
    let mut c1 = 0.0;
    let mut c2 = 0.0;
    let mut n = 0;
    let assign = |v: Index,
                  l: bool,
                  labels: &mut Vec<Option<bool>>,
                  c1: &mut f32,
                  c2: &mut f32,
                  n: &mut usize| {
        let old = replace(&mut labels[v as usize], Some(l));
        let c = h.capacity(v);
        match (old, l) {
            (Some(false), true) => {
                *c1 += c;
                *c2 -= c;
            }
            (Some(true), false) => {
                *c1 -= c;
                *c2 += c;
            }
            (None, true) => {
                *c1 += c;
                *n += 1;
            }
            (None, false) => {
                *c2 += c;
                *n += 1;
            }
            _ => {}
        }
    };
    assign(v1, true, &mut labels, &mut c1, &mut c2, &mut n);
    assign(v2, false, &mut labels, &mut c1, &mut c2, &mut n);

    let mut n1: Vec<_> = h.incident_pins(v1).collect();
    let mut n2: Vec<_> = h.incident_pins(v2).collect();
    let mut rng = thread_rng();
    n1.shuffle(&mut rng);
    n2.shuffle(&mut rng);
    n1.truncate(TAU);
    n2.truncate(TAU);
    for v in n1 {
        assign(v, true, &mut labels, &mut c1, &mut c2, &mut n);
    }
    for v in n2 {
        assign(v, false, &mut labels, &mut c1, &mut c2, &mut n);
    }

    let size_constraint = h.size_constraint(epsilon);
    while n < h.num_pins() {
        for v in h.pins() {
            if labels[v as usize].is_none() {
                let mut w1 = 0.0;
                let mut w2 = 0.0;
                for n in h.incident_nets(v) {
                    let w = h.weight(n);
                    if h.pins_in_net(n).any(|p| labels[p as usize] == Some(true)) {
                        w1 += w;
                    }
                    if h.pins_in_net(n).any(|p| labels[p as usize] == Some(false)) {
                        w2 += w;
                    }
                }
                let valid1 = c1 + h.capacity(v) < size_constraint;
                let valid2 = c1 + h.capacity(v) < size_constraint;
                assign(
                    v,
                    valid1 && w1 >= w2 || !valid2,
                    &mut labels,
                    &mut c1,
                    &mut c2,
                    &mut n,
                );
            }
        }
    }

    labels.into_iter().map(|l| l.unwrap_or(false)).collect()
}

fn pseudo_peripheral_vertices(h: &Bipartite) -> (Index, Index) {
    let last_bfs = |start_v| {
        let mut queue = VecDeque::new();
        queue.push_back(start_v);
        let mut visited = bitvec![usize, Lsb0; 0; h.pin_index_space_size()];
        visited.set(start_v as usize, true);

        let mut last = start_v;
        while let Some(pop) = queue.pop_front() {
            for neighbor in h.incident_pins(pop) {
                if !visited[neighbor as usize] {
                    visited.set(neighbor as usize, true);
                    queue.push_back(neighbor);
                }
            }
            last = pop;
        }
        last
    };

    let v1 = random_vertex(h);
    let v2 = last_bfs(v1);
    let v3 = last_bfs(v2);
    (v2, v3)
}

fn random_vertex(h: &Bipartite) -> Index {
    loop {
        let v: Index = random::<Index>() % h.pin_index_space_size() as Index;
        if h.enabled(v) {
            return v;
        }
    }
}
