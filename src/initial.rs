use std::collections::VecDeque;

use bitvec::prelude::*;

use crate::bipartite::*;

pub(crate) fn random_partioning(h: &Bipartite) -> Bipartition {
    let mut bipartition = (vec![], vec![]);

    for v in 0..h.v.len() {
        if todo!() {
            bipartition.0.push(v as Index);
        } else {
            bipartition.1.push(v as Index);
        }
    }

    bipartition
}

pub(crate) fn bfs_half(h: &Bipartite) -> Bipartition {
    let start_v: Index = todo!();
    let mut queue = VecDeque::new();
    queue.push_back(start_v);
    let mut visited = bitvec![usize, Lsb0; 0; h.v.len()];
    visited.set(start_v as usize, true);

    let mut bipartition = (vec![], vec![]);
    while let Some(pop) = queue.pop_front() {
        for neighbor in h.incident_pins(pop) {
            if !visited[neighbor as usize] {
                visited.set(neighbor as usize, true);
                queue.push_back(neighbor);
            }
        }

        if bipartition.0.len() < h.v.len() / 2 {
            bipartition.0.push(pop);
        } else {
            bipartition.1.push(pop);
        }
    }

    bipartition
}

const TAU: usize = 5;

pub(crate) fn size_constrained_label_propagation(h: &Bipartite) -> Bipartition {
    let mut bipartition = (vec![], vec![]);
    let (v1, v2) = pseudo_peripheral_vertices(h);
    let mut labels = vec![None; h.v.len()];
    
    todo!();

    for (v, l) in labels.into_iter().enumerate() {
        if l.unwrap() {
            bipartition.0.push(v as Index);
        } else{
            bipartition.1.push(v as Index);
        }
    }
    bipartition
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

    let v1: Index = todo!();
    let v2 = last_bfs(v1);
    let v3 = last_bfs(v2);
    (v1, v2)
}
