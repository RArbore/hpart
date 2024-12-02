use crate::bipartite::*;

pub(crate) fn coarsen(h: &mut Bipartite) {
    todo!()
}

/// Constants from Section 5 from Schlag '2015.
const T: usize = 320;
const S: f32 = 3.25;

fn needs_coarsening(h: &Bipartite) -> bool {
    todo!()
}

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
