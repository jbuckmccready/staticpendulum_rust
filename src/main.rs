extern crate time;
extern crate image;
mod integrators;
mod differential_systems;
use differential_systems::*;
use integrators::*;
use time::{PreciseTime};
use std::fs::File;
use std::path::Path;

pub struct IntegrationResult {
    pub converge_result: i32,
    pub converge_time: f64,
    pub step_count: u32
}

#[inline]
fn integrate_point<T: Integrator<PendulumSystem>>(integrator: &T, pendulum_system: &PendulumSystem, point: &mut [f64]) -> IntegrationResult {
    fn is_near_attractor(attr_x: f64, attr_y: f64, curr_x: f64, curr_y: f64, tol: f64) -> bool {
        return (attr_x - tol < curr_x) && (curr_x < attr_x + tol) &&
                (attr_y  - tol < curr_y) && (curr_y < attr_y + tol);
    }

    fn is_near_mid(curr_x: f64, curr_y: f64, tol: f64) -> bool {
        return (-tol < curr_x) && (curr_x < tol) && (-tol < curr_y) && (curr_y < tol);
    }

    let pos_tol = 0.5;
    let mid_tol = 0.1;
    let time_tol = 5.0;

    let mut curr_time = 0.0;
    let mut step_size = 0.01;
    let mut trial_count = 0;
    let mut initial_time_found = 0.0;
    let mut current_attractor: i32 = -2;
    let mut step_count = 0u32;
    let mut near_attractor = false;
    let attractor_count = pendulum_system.attractors.len();
    //println!("curr_time: {}", curr_time);
    unsafe {
    'outer: while trial_count < 1000 {
        step_count += integrator.do_step(&pendulum_system, &mut *point, &mut curr_time, &mut step_size) as u32;
       // println!("curr_time: {}", curr_time);
        trial_count += 1;
        'inner: for i in 0..attractor_count {
            near_attractor = is_near_attractor(pendulum_system.attractors[i].x_position,
                                               pendulum_system.attractors[i].y_position,
                                               *point.get_unchecked(0), *point.get_unchecked(1), pos_tol);
            if near_attractor {
                //println!("near attractor {}!", i);
                if current_attractor == i as i32 {
                    if curr_time - initial_time_found > time_tol {
                        //println!("Converge time: {}", curr_time);
                        //println!("Converge index: {}", i);
                        return IntegrationResult { converge_result: current_attractor, converge_time: curr_time, step_count: step_count };
                        //break 'outer;
                    }
                } else {
                    current_attractor = i as i32;
                    initial_time_found = curr_time;
                    //println!("{}", current_attractor);
                }
                break 'inner;
            }
        }

        if !near_attractor && is_near_mid(*point.get_unchecked(0), *point.get_unchecked(1), mid_tol) {
            if current_attractor == -1 {
                if curr_time - initial_time_found > time_tol {
                    //println!("Converge time: {}", curr_time);
                    //println!("Converge index: {}", -1);
                    return IntegrationResult { converge_result: current_attractor, converge_time: curr_time, step_count: step_count };
                    //break 'outer;
                }
            } else {
                current_attractor = -1;
                initial_time_found = curr_time;
            }

        }

        near_attractor = false;
    }
    }
    return IntegrationResult { converge_result: current_attractor, converge_time: curr_time, step_count: step_count };
    //return current_attractor;
}


fn main() {
    let test_sys = PendulumSystem {
        height: 0.05,
        mass: 1.0,
        gravity: 9.8,
        drag: 0.2,
        length: 10.0,
        attractors: vec![
            Attractor {x_position: -0.5, y_position: 0.86602540378, force_coefficient: 1.0},
            Attractor {x_position: -0.5, y_position: -0.86602540378, force_coefficient: 1.0},
            Attractor {x_position: 1.0, y_position: 0.0, force_coefficient: 1.0}
        ]
    };

    let test_integrator = CashKarp54 { rel_tol: 1.0e-6, abs_tol: 1.0e-6, max_step_size: 0.1 };
    //let mut x: Vec<f64> = vec![0.001,0.001,0.0,0.0];
    //integrate_point(&test_integrator, &test_sys, &mut x);
    let res = 0.1;
    let extent = 7.5;
    let dim = (extent / res) as i32;
    let mut imgbuf = image::ImageBuffer::new(2*dim as u32 + 1, 2*dim as u32 + 1);

    println!("integrating {0} points", (2*dim)*(2*dim));
    let start = PreciseTime::now();
    for x in -dim..dim {
        for y in -dim..dim {
            let mut point = vec![(x as f64) * res, (y as f64) * res, 0.0, 0.0];
            let result = integrate_point(&test_integrator, &test_sys, &mut point);
            let x_pixel = (x + dim) as u32;
            let y_pixel = (dim - y) as u32;
            match result.converge_result {
                -1 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([255, 255, 255])),
                0 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([255, 140, 0])),
                1 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([30, 144, 255])),
                2 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([178, 34, 34])),
                _ => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([0, 0, 0]))
            }
            // point_map.push(vec![(x as f64) * res, (y as f64) * res, 0.0, 0.0]);
        }
    }
    let time = start.to(PreciseTime::now());

    let ref mut fout = File::create(&Path::new("fractal.png")).unwrap();

    let _    = image::ImageRgb8(imgbuf).save(fout, image::PNG);

    println!("timing: {}", time.num_milliseconds());
}
