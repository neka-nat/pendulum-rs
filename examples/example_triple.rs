extern crate pendulum_rs as pn;
extern crate nalgebra as na;
use plotters::prelude::*;
use na::Vector3;
use pn::triple::*;

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut th = Vector3::<f64>::from([0.0, 0.0, 0.0]);
    let mut thdot = Vector3::<f64>::from([0.0, 0.0, 0.0]);
    let l: Vector3<f64> = Vector3::<f64>::from([1.0, 1.0, 1.0]);
    let lgc: Vector3<f64> = Vector3::<f64>::from([0.5, 0.5, 0.5]);
    let mass: Vector3<f64> = Vector3::<f64>::from([0.5, 0.5, 0.5]);
    let inertia: Vector3<f64> = Vector3::<f64>::from([mass[0] * l[0].powi(2),
                                                      mass[1] * l[1].powi(2),
                                                      mass[2] * l[2].powi(2)]) / 12.0;
    let dt = 0.005;
    let root = BitMapBackend::gif("results/animation_triple.gif", (600, 600), 10)?
        .into_drawing_area();
    for _ in 0..800 {
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .build_cartesian_2d(-3.5..3.5, -3.5..3.5)?;

        let th2dot = forward_triple_dyn(&th, &thdot, &l, &lgc, &mass, &inertia);
        let pos = forward_triple_kin(&th, &l);
        let pend_vertices = vec![(0.0, 0.0), (pos[0][0], pos[0][1]), (pos[1][0], pos[1][1]), (pos[2][0], pos[2][1])];
        th = th + thdot * dt;
        thdot = thdot + th2dot * dt;
        println!("{0} {1} {2}", th, thdot, th2dot);


        chart.draw_series(LineSeries::new(pend_vertices.clone(), &BLACK))?;
        chart.draw_series(
            pend_vertices
                .iter()
                .map(|(x, y)| Circle::new((*x, *y), 5, RED.filled())),
        )?;

        root.present()?;
    }

    Ok(())
}

