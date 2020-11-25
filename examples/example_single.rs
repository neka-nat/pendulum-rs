extern crate pendulum_rs as pn;
use plotters::prelude::*;
use pn::single::*;

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut th = std::f64::consts::PI / 2.0;
    let mut thdot = 0.0;
    let dt = 0.01;
    let root = BitMapBackend::gif("results/animation_single.gif", (800, 600), 1000)?
        .into_drawing_area();
    for _ in 0..100 {
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .build_cartesian_2d(-2.0..2.0, -1.5..1.5)?;

        let th2dot = forward_single_dyn(th, 1.0);
        let pos = forward_single_kin(th, 1.0);
        let pend_vertices = vec![(0.0, 0.0), (pos[1], -pos[0])];
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
