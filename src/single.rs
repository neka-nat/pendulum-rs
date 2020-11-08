extern crate nalgebra as na;

use na::Vector2;

static GRAVITY: f64 = 9.8;

pub fn forward_single_kin(state: f64, l: f64) -> Vector2<f64> {
    let pos = Vector2::<f64>::from([l * state.cos(), l * state.sin()]);
    return pos;
}

pub fn forward_single_dyn(state: f64, l: f64) -> f64 {
    let th2dot = -(GRAVITY / l) * state.sin();
    return th2dot;
}
