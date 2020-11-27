extern crate nalgebra as na;
use na::{Vector2, Vector3, Matrix3};
use crate::double::*;

static GRAVITY: f64 = 9.8;

pub fn forward_triple_kin(state: &Vector3<f64>, l: &Vector3<f64>) -> [Vector2<f64>; 3] {
    let pos01 = forward_double_kin(&Vector2::<f64>::from([state[0], state[1]]),
                                   &Vector2::<f64>::from([l[0], l[1]]));
    let pos2 = Vector2::<f64>::from([pos01[1][0] + l[2] * state.sum().cos(),
                                     pos01[1][1] + l[2] * state.sum().sin()]);
    return [pos01[0], pos01[1], pos2];
}

pub fn forward_triple_dyn(state: &Vector3<f64>,
                          state_dot: &Vector3<f64>,
                          l: &Vector3<f64>,
                          lgc: &Vector3<f64>,
                          mass: &Vector3<f64>,
                          inertia: &Vector3<f64>) -> Vector3<f64> {
    let c1 = state[0].cos();
    let c2 = state[1].cos();
    let c3 = state[2].cos();
    let c12 = (state[0] + state[1]).cos();
    let c23 = (state[1] + state[2]).cos();
    let c123 = state.sum().cos();
    let s2 = state[1].sin();
    let s3 = state[2].sin();
    let s23 = (state[1] + state[2]).sin();

    let m00_1 = mass[0] * lgc[0].powi(2) + inertia[0];
    let m00_2 = mass[1] * (l[0].powi(2) + lgc[1].powi(2) + 2.0 * l[0] * lgc[1] * c2) + inertia[1];
    let m00_3 = mass[2] * (l[0].powi(2) + l[1].powi(2) + lgc[2].powi(2) + 2.0 * l[0] * l[1] * c2 + 2.0 * l[1] * lgc[2] * c3 + 2.0 * l[0] * lgc[2] * c23) + inertia[2];
    let m00 = m00_1 + m00_2 + m00_3;

    let m01_1 = mass[1] * (lgc[1].powi(2) + l[0] * lgc[1] * c2) + inertia[1];
    let m01_2 = mass[2] * (l[1].powi(2) + lgc[2].powi(2) + l[0] * l[1] * c2 + 2.0 * l[1] * lgc[2] * c3 + l[0] * lgc[2] * c23) + inertia[2];
    let m01 = m01_1 + m01_2;

    let m02 = mass[2] * (lgc[2].powi(2) + l[1] * lgc[2] * c3 + l[0] * lgc[2] * c23) + inertia[2];
    let m10 = m01.clone();
    let m11 = mass[1] * lgc[1].powi(2) + inertia[1] + mass[2] * (l[1].powi(2) + lgc[2].powi(2) + 2.0 * l[1] * lgc[2] * c3) + inertia[2];
    let m12 = mass[2] * (lgc[2].powi(2) + l[1] * lgc[2] * c3) + inertia[2];
    let m20 = m02.clone();
    let m21 = m12.clone();
    let m22 = mass[2] * lgc[2].powi(2) + inertia[2];

    let htmp2 = -mass[1] * l[0] * lgc[1] * s2;
    let htmp312 = -mass[2] * l[0] * l[1] * s2;
    let htmp323 = -mass[2] * l[1] * lgc[2] * s3;
    let htmp313 = -mass[2] * l[0] * lgc[2] * s23;
    let h0_1 = htmp2 * (2.0 * state_dot[0] + state_dot[1]) * state_dot[1] + htmp312 * (2.0 * state_dot[0] + state_dot[1]) * state_dot[1];
    let h0_2 = htmp323 * (2.0 * (state_dot[0] + state_dot[1]) + state_dot[2]) * state_dot[2] + htmp313 * (state_dot[1] + state_dot[2]) * (2.0 * state_dot[0] + state_dot[1] + state_dot[2]);
    let h0 = h0_1 + h0_2;

    let h1_1 = -htmp2 * state_dot[0].powi(2) - htmp312 * state_dot[0].powi(2);
    let h1_2 = htmp323 * state_dot[2] * (2.0 * (state_dot[0] + state_dot[1]) + state_dot[2]) - htmp313 * state_dot[0].powi(2);
    let h1 = h1_1 + h1_2;

    let h2 = -htmp323 * (state_dot[0] + state_dot[1]).powi(2) - htmp313 * state_dot[0].powi(2);

    let g0 = mass[0] * GRAVITY * lgc[0] * c1 + mass[1] * GRAVITY * (l[0] * c1 + lgc[1] * c12) + mass[2] * GRAVITY * (l[0] * c1 + l[1] * c12 + lgc[2] * c123);
    let g1 = mass[1] * GRAVITY * lgc[1] * c12 + mass[2] * GRAVITY * (l[1] * c12 + lgc[2] * c123);
    let g2 = mass[2] * GRAVITY * lgc[2] * c123;

    let inv_mmat_op = Matrix3::<f64>::from([[m00, m01, m02],
                                            [m10, m11, m12],
                                            [m20, m21, m22]]).try_inverse();
    if let Some(inv_mmat) = inv_mmat_op {
        let state_2dot = -inv_mmat * Vector3::<f64>::from([h0 + g0, h1 + g1, h2 + g2]);
        return state_2dot;
    } else {
        return Vector3::<f64>::zeros();
    }
}
