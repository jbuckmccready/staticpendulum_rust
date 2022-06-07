use std::time::Instant;

use itertools::iproduct;
use rayon::prelude::*;
mod differential_systems;
mod integrators;
use differential_systems::*;
use integrators::*;

pub struct IntegrationResult {
    pub converge_result: i32,
    pub converge_time: f64,
    pub step_count: u32,
}

#[inline]
fn integrate_point<T: Integrator<PendulumSystem>>(
    integrator: &T,
    pendulum_system: &PendulumSystem,
    point: &mut [f64; 4],
) -> IntegrationResult {
    fn is_near_attractor(attr_x: f64, attr_y: f64, curr_x: f64, curr_y: f64, tol: f64) -> bool {
        (attr_x - tol < curr_x)
            && (curr_x < attr_x + tol)
            && (attr_y - tol < curr_y)
            && (curr_y < attr_y + tol)
    }

    fn is_near_mid(curr_x: f64, curr_y: f64, tol: f64) -> bool {
        (-tol < curr_x) && (curr_x < tol) && (-tol < curr_y) && (curr_y < tol)
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
    while trial_count < 1000 {
        step_count += integrator.do_step(
            &pendulum_system,
            &mut *point,
            &mut curr_time,
            &mut step_size,
        ) as u32;
        trial_count += 1;
        'inner: for i in 0..attractor_count {
            near_attractor = is_near_attractor(
                pendulum_system.attractors[i].x_position,
                pendulum_system.attractors[i].y_position,
                point[0],
                point[1],
                pos_tol,
            );
            if near_attractor {
                if current_attractor == i as i32 {
                    if curr_time - initial_time_found > time_tol {
                        return IntegrationResult {
                            converge_result: current_attractor,
                            converge_time: curr_time,
                            step_count,
                        };
                    }
                } else {
                    current_attractor = i as i32;
                    initial_time_found = curr_time;
                }
                break 'inner;
            }
        }

        if !near_attractor && is_near_mid(point[0], point[1], mid_tol) {
            if current_attractor == -1 {
                if curr_time - initial_time_found > time_tol {
                    return IntegrationResult {
                        converge_result: current_attractor,
                        converge_time: curr_time,
                        step_count,
                    };
                }
            } else {
                current_attractor = -1;
                initial_time_found = curr_time;
            }
        }
        near_attractor = false;
    }
    IntegrationResult {
        converge_result: current_attractor,
        converge_time: curr_time,
        step_count,
    }
}

fn main() {
    let test_sys = PendulumSystem::new(
        0.05,
        1.0,
        9.8,
        0.2,
        10.0,
        vec![
            Attractor {
                x_position: -0.5,
                y_position: 0.866_025_403_78,
                force_coefficient: 1.0,
            },
            Attractor {
                x_position: -0.5,
                y_position: -0.866_025_403_78,
                force_coefficient: 1.0,
            },
            Attractor {
                x_position: 1.0,
                y_position: 0.0,
                force_coefficient: 1.0,
            },
        ],
    );

    let test_integrator = CashKarp54 {
        rel_tol: 1.0e-6,
        abs_tol: 1.0e-6,
        max_step_size: 0.1,
    };
    let res = 0.1;
    let extent = 7.5;
    let dim = (extent / res) as i32;
    let mut imgbuf = image::ImageBuffer::new(2 * dim as u32 + 1, 2 * dim as u32 + 1);

    println!("integrating {0} points", (2 * dim) * (2 * dim));
    let start = Instant::now();
    let mut points: Vec<[f64; 4]> = iproduct!(-dim..dim, -dim..dim)
        .map(|(x, y)| [f64::from(x) * res, f64::from(y) * res, 0.0, 0.0])
        .collect();

    let results: Vec<IntegrationResult> = points
        .par_iter_mut()
        .map(|point| integrate_point(&test_integrator, &test_sys, point))
        .collect();

    let mut index = 0;
    for x in -dim..dim {
        for y in -dim..dim {
            let x_pixel = (x + dim) as u32;
            let y_pixel = (dim - y) as u32;
            let result = &results[index];
            match result.converge_result {
                -1 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([255u8, 255, 255])),
                0 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([255, 140, 0])),
                1 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([30, 144, 255])),
                2 => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([178, 34, 34])),
                _ => imgbuf.put_pixel(x_pixel, y_pixel, image::Rgb([0, 0, 0])),
            }
            index += 1;
        }
    }

    let time = Instant::now() - start;

    imgbuf.save("fractal.png").unwrap();

    println!("timing: {}ms", time.as_millis());
}
