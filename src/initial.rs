use std::collections::VecDeque;
use std::mem::replace;

use bitvec::prelude::*;
use rand::prelude::*;

use crate::bipartite::*;

fn random_partioning(h: &Bipartite) -> Bipartition {
    (0..h.v.len()).map(|_| random()).collect()
}

fn bfs_half(h: &Bipartite) -> Bipartition {
    let start_v: Index = random::<Index>() % h.v.len() as Index;
    let mut queue = VecDeque::new();
    queue.push_back(start_v);
    let mut visited = bitvec![usize, Lsb0; 0; h.v.len()];
    visited.set(start_v as usize, true);

    let mut bipartition = vec![false; h.v.len()];
    let mut num_visited = 0;
    while let Some(pop) = queue.pop_front() {
        for neighbor in h.incident_pins(pop) {
            if !visited[neighbor as usize] {
                visited.set(neighbor as usize, true);
                queue.push_back(neighbor);
            }
        }

        if num_visited < h.v.len() / 2 {
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
    let mut labels = vec![None; h.v.len()];
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
        let c = h.c[v as usize];
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

    let size_constraint = (1.0 + epsilon) * (h.total_capacity() as f32 / 2.0);
    while n < h.v.len() {
        for v_idx in 0..h.v.len() {
            if labels[v_idx].is_none() {
                let mut w1 = 0.0;
                let mut w2 = 0.0;
                for n in h.incident_nets(v_idx as Index) {
                    let w = h.w[n as usize];
                    if h.pins_in_net(n).any(|p| labels[p as usize] == Some(true)) {
                        w1 += w;
                    }
                    if h.pins_in_net(n).any(|p| labels[p as usize] == Some(false)) {
                        w2 += w;
                    }
                }
                let valid1 = c1 + h.c[v_idx] < size_constraint;
                let valid2 = c1 + h.c[v_idx] < size_constraint;
                assign(
                    v_idx as Index,
                    valid1 && w1 >= w2 || !valid2,
                    &mut labels,
                    &mut c1,
                    &mut c2,
                    &mut n,
                );
            }
        }
    }

    labels.into_iter().map(|l| l.unwrap()).collect()
}

fn pseudo_peripheral_vertices(h: &Bipartite) -> (Index, Index) {
    let last_bfs = |start_v| {
        let mut queue = VecDeque::new();
        queue.push_back(start_v);
        let mut visited = bitvec![usize, Lsb0; 0; h.v.len()];
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

    let v1: Index = random::<Index>() % h.v.len() as Index;
    let v2 = last_bfs(v1);
    let v3 = last_bfs(v2);
    (v2, v3)
}
