extern crate nalgebra as na;
use na::{Vector2, Matrix2};
use crate::single::*;

static GRAVITY: f64 = 9.8;

pub fn forward_double_kin(state: &Vector2<f64>, l: &Vector2<f64>) -> [Vector2<f64>; 2] {
    let pos0 = forward_single_kin(state[0], l[0]);
    let pos1 = Vector2::<f64>::from([pos0[0] + l[1] * state.sum().cos(),
                                     pos0[1] + l[1] * state.sum().sin()]);
    return [pos0, pos1];
}

pub fn forward_double_dyn(state: &Vector2<f64>,
                          state_dot: &Vector2<f64>,
                          l: &Vector2<f64>,
                          lgc: &Vector2<f64>,
                          mass: &Vector2<f64>,
                          inertia: &Vector2<f64>) -> Vector2<f64> {
    let c1 = state[0].cos();
    let c2 = state[1].cos();
    let s2 = state[1].sin();
    let c12 = (state[0] + state[1]).cos();
    let m00 = mass[0] * lgc[0].powi(2) + inertia[0] + mass[1] * (l[0].powi(2) + lgc[1].powi(2) + 2.0 * l[0] * lgc[1] * c2) + inertia[1];
    let m01 = mass[1] * (lgc[1].powi(2) + l[0] * lgc[1] * c2) + inertia[1];
    let m10 = m01.clone();
    let m11 = mass[1] * lgc[1].powi(2) + inertia[1];
    let htmp = -mass[1] * l[0] * lgc[1] * s2;
    let h0 = htmp * (state_dot[1].powi(2) + 2.0 * state_dot[0] * state_dot[1]);
    let h1 = -htmp * state_dot[0].powi(2);
    let g0 = mass[0] * GRAVITY * lgc[0] * c1 + mass[1] * GRAVITY * (l[0] * c1 + lgc[1] * c12);
    let g1 = mass[1] * GRAVITY * lgc[1] * c12;
    let inv_mmat_op = Matrix2::<f64>::from([[m00, m01], [m10, m11]]).try_inverse();
    if let Some(inv_mmat) = inv_mmat_op {
        let state_2dot = -inv_mmat * Vector2::<f64>::from([h0 + g0, h1 + g1]);
        return state_2dot;
    } else {
        return Vector2::<f64>::from([0.0, 0.0]);
    }
}
