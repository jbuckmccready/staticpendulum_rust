use differential_systems::DifferentialSystem;
use std::f64;

pub trait Integrator<T: DifferentialSystem> {
    fn do_step(&self, system: &T, x: &mut [f64], t: &mut f64, h: &mut f64) -> i32;
}

pub struct CashKarp54 {
    pub rel_tol: f64,
    pub abs_tol: f64,
    pub max_step_size: f64
}

impl<T: DifferentialSystem> Integrator<T> for CashKarp54 {
    #[inline]
    fn do_step(&self, system: &T, x: &mut [f64], t: &mut f64, h: &mut f64) -> i32 {

        unsafe {
            const C2: f64 = 1.0 / 5.0;
            const C3: f64 = 3.0 / 10.0;
            const C4: f64 = 3.0 / 5.0;
            // const C5: f64 = 1.0; 
            const C6: f64 = 7.0 / 8.0;

            const B5O1: f64 = 37.0 / 378.0;
            const B5O2: f64 = 0.0;
            const B5O3: f64 = 250.0 / 621.0;
            const B5O4: f64 = 125.0 / 594.0;
            const B5O5: f64 = 0.0;
            const B5O6: f64 = 512.0 / 1771.0;

            const B4O1: f64 = 2825.0 / 27648.0;
            const B4O2: f64 = 0.0;
            const B4O3: f64 = 18575.0 / 48384.0;
            const B4O4: f64 = 13525.0 / 55296.0;
            const B4O5: f64 = 277.0 / 14336.0;
            const B4O6: f64 = 1.0 / 4.0;

            const BDIFF1: f64 = B5O1 - B4O1;
            const BDIFF2: f64 = B5O2 - B4O2;
            const BDIFF3: f64 = B5O3 - B4O3;
            const BDIFF4: f64 = B5O4 - B4O4;
            const BDIFF5: f64 = B5O5 - B4O5;
            const BDIFF6: f64 = B5O6 - B4O6;

            const A21: f64 = 1.0 / 5.0;
            const A31: f64 = 3.0 / 40.0;
            const A32: f64 = 9.0 / 40.0;
            const A41: f64 = 3.0 / 10.0;
            const A42: f64 = -9.0 / 10.0;
            const A43: f64 = 6.0 / 5.0;
            const A51: f64 = -11.0 / 54.0;
            const A52: f64 = 5.0 / 2.0;
            const A53: f64 = -70.0 / 27.0;
            const A54: f64 = 35.0 / 27.0;
            const A61: f64 = 1631.0 / 55296.0;
            const A62: f64 = 175.0 / 512.0;
            const A63: f64 = 575.0 / 13824.0;
            const A64: f64 = 44275.0 / 110592.0;
            const A65: f64 = 253.0 / 4096.0;

            let state_size = x.len();
            let mut temp_state = vec![0.0; state_size];

            let mut k1: Vec<f64> = vec![0.0; state_size];
            system.evaluate(&x, &mut k1, *t); // fill k1

            for i in 0..state_size {
                *temp_state.get_unchecked_mut(i) = *x.get_unchecked(i) + *h * A21 * *k1.get_unchecked(i);
            }
            let mut k2: Vec<f64> = vec![0.0; state_size];
            system.evaluate(&temp_state, &mut k2, *t + C2 * *h); // fill k2

            for i in 0..state_size {
                *temp_state.get_unchecked_mut(i) = *x.get_unchecked(i) + *h * (A31 * *k1.get_unchecked(i) + A32 * *k2.get_unchecked(i));
            }
            let mut k3: Vec<f64> = vec![0.0; state_size];
            system.evaluate(&temp_state, &mut k3, *t + C3 * *h); // fill k3

            for i in 0..state_size {
                *temp_state.get_unchecked_mut(i) = *x.get_unchecked(i) + *h * (A41 * *k1.get_unchecked(i) + A42 * *k2.get_unchecked(i) + A43 * *k3.get_unchecked(i));
            }
            let mut k4: Vec<f64> = vec![0.0; state_size];
            system.evaluate(&temp_state, &mut k4, *t + C4 * *h); // fill k4

            for i in 0..state_size {
                *temp_state.get_unchecked_mut(i) = *x.get_unchecked(i) + *h * (A51 * *k1.get_unchecked(i) + A52 * *k2.get_unchecked(i) + A53 * *k3.get_unchecked(i) + A54 * *k4.get_unchecked(i));
            }
            let mut k5: Vec<f64> = vec![0.0; state_size];
            system.evaluate(&temp_state, &mut k5, *t + *h); // fill k5

            for i in 0..state_size {
                *temp_state.get_unchecked_mut(i) = *x.get_unchecked(i) + *h * (A61 * *k1.get_unchecked(i) + A62 * *k2.get_unchecked(i) + A63 * *k3.get_unchecked(i) +
                                             A64 * *k4.get_unchecked(i) + A65 * *k5.get_unchecked(i));
            }
            let mut k6: Vec<f64> = vec![0.0; state_size];
            system.evaluate(&temp_state, &mut k6, *t + C6 * *h); // fill k6

            let mut order5_solution: Vec<f64> = vec![0.0; state_size];
            for i in 0..state_size {
                *order5_solution.get_unchecked_mut(i) = *h * (B5O1 * *k1.get_unchecked(i) + B5O2 * *k2.get_unchecked(i) + B5O3 * *k3.get_unchecked(i) +
                                           B5O4 * *k4.get_unchecked(i) + B5O5 * *k5.get_unchecked(i) + B5O6 * *k6.get_unchecked(i));
            }

            for i in 0..state_size {
                *temp_state.get_unchecked_mut(i) = *h * (BDIFF1 * *k1.get_unchecked(i) + BDIFF2 * *k2.get_unchecked(i) + BDIFF3 * *k3.get_unchecked(i) +
                                      BDIFF4 * *k4.get_unchecked(i) + BDIFF5 * *k5.get_unchecked(i) + BDIFF6 * *k6.get_unchecked(i));
            }
            let mut potential_solution: Vec<f64> = vec![0.0; state_size];

            for i in 0..state_size {
                *potential_solution.get_unchecked_mut(i) = *x.get_unchecked(i) + *order5_solution.get_unchecked(i)
            }

            let mut error_values: Vec<f64> = vec![0.0; state_size];
            for i in 0..state_size {
                *error_values.get_unchecked_mut(i) = (*temp_state.get_unchecked(i) / (self.abs_tol + self.rel_tol * *potential_solution.get_unchecked(i))).abs();
            }

            let max_error: f64 = error_values.into_iter().fold(f64::NAN, |a, b| if a > b { a } else { b });

            if max_error > 1.0 {
                *h *= (0.9 * max_error.powf(-0.25)).max(0.2);
                return 0;
            }

            *t += *h;
            for i in 0..state_size {
                *x.get_unchecked_mut(i) = *potential_solution.get_unchecked(i);
            }

            if max_error < 0.5 {
                *h = (*h * (0.9 * max_error.powf(-0.20)).min(5.0)).min(self.max_step_size);
            }
            return 1;
        }
    }
}
