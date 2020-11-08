extern crate pendulum_rs as pn;
use pn::single::*;

pub fn main() {
    let pos = forward_single_kin(0.0, 1.0);
    println!("{:?}", pos);
}
