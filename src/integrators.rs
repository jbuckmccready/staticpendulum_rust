use differential_systems::DifferentialSystem;
use std::f64;

pub trait Integrator<T: DifferentialSystem> {
    fn do_step(&self, system: &T, x: &mut [f64; 4], t: &mut f64, h: &mut f64) -> i32;
}

pub struct CashKarp54 {
    pub rel_tol: f64,
    pub abs_tol: f64,
    pub max_step_size: f64
}

impl<T: DifferentialSystem> Integrator<T> for CashKarp54 {
    fn do_step(&self, system: &T, x: &mut [f64; 4], t: &mut f64, h: &mut f64) -> i32 {

        const C2: f64 = 1.0 / 5.0;
        const C3: f64 = 3.0 / 10.0;
        const C4: f64 = 3.0 / 5.0;
        const C5: f64 = 1.0;
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
        const A64: f64 = 44275.0 / 110_592.0;
        const A65: f64 = 253.0 / 4096.0;

        let mut temp_state = [0.0; 4];

        let mut k1 = [0.0; 4];
        system.evaluate(&x, &mut k1, *t); // fill k1

        izip!((*x).iter(), k1.iter(), &mut temp_state)
            .for_each(|(xi, k1i, ts)| { *ts = xi + *h *A21 * k1i });
        let mut k2 = [0.0; 4];
        system.evaluate(&temp_state, &mut k2, *t + C2 * *h); // fill k2

        izip!((*x).iter(), k1.iter(), k2.iter(), &mut temp_state)
            .for_each(|(xi, k1i, k2i, ts)| { *ts = xi + *h * (A31 * k1i + A32 * k2i) });
        let mut k3 = [0.0; 4];
        system.evaluate(&temp_state, &mut k3, *t + C3 * *h); // fill k3

        izip!((*x).iter(), k1.iter(), k2.iter(), k3.iter(), &mut temp_state)
            .for_each(|(xi, k1i, k2i, k3i, ts)| { *ts = xi + *h * (A41 * k1i + A42 * k2i + A43 * k3i) });
        let mut k4 = [0.0; 4];
        system.evaluate(&temp_state, &mut k4, *t + C4 * *h); // fill k4

        izip!((*x).iter(), k1.iter(), k2.iter(), k3.iter(), k4.iter(), &mut temp_state)
            .for_each(|(xi, k1i, k2i, k3i, k4i, ts)| { *ts = xi + *h * (A51 * k1i + A52 * k2i + A53 * k3i + A54 * k4i) });
        let mut k5 = [0.0; 4];
        system.evaluate(&temp_state, &mut k5, *t + C5 * *h); // fill k5

        izip!((*x).iter(), k1.iter(), k2.iter(), k3.iter(), k4.iter(), k5.iter(), &mut temp_state)
            .for_each(|(xi, k1i, k2i, k3i, k4i, k5i, ts)| { *ts = xi + *h * (A61 * k1i + A62 * k2i + A63 * k3i + A64 * k4i + A65 * k5i) });
        let mut k6 = [0.0; 4];
        system.evaluate(&temp_state, &mut k6, *t + C6 * *h); // fill k6

        let mut order5_solution = [0.0; 4];
        izip!(k1.iter(), k2.iter(), k3.iter(), k4.iter(), k5.iter(), k6.iter(), &mut order5_solution)
            .for_each(|(k1i, k2i, k3i, k4i, k5i, k6i, o5s)|{ *o5s = *h * (B5O1 * k1i + B5O2 * k2i + B5O3 * k3i + B5O4 * k4i + B5O5 * k5i + B5O6 * k6i) });

        izip!(k1.iter(), k2.iter(), k3.iter(), k4.iter(), k5.iter(), k6.iter(), &mut temp_state)
            .for_each(|(k1i, k2i, k3i, k4i, k5i, k6i, ts)| { *ts = *h * (BDIFF1 * k1i + BDIFF2 * k2i + BDIFF3 * k3i + BDIFF4 * k4i + BDIFF5 * k5i + BDIFF6 * k6i) });

        let mut potential_solution = [0.0; 4];
        izip!((*x).iter(), order5_solution.iter(), &mut potential_solution)
            .for_each(|(xi, o5s, ps)| { *ps = xi + o5s} );

        let mut error_values = [0.0; 4];
        izip!(temp_state.iter(), potential_solution.iter(), &mut error_values)
            .for_each(|(ts, ps, ev)| { *ev = (ts / (self.abs_tol + self.rel_tol * ps)).abs() });

        let max_error: f64 = error_values.iter().fold(f64::NAN, |a, b| if a > *b { a } else { *b });

        if max_error > 1.0 {
            *h *= (0.9 * max_error.powf(-0.25)).max(0.2);
            return 0;
        }

        *t += *h;
        *x = potential_solution;

        if max_error < 0.5 {
            *h = (*h * (0.9 * max_error.powf(-0.20)).min(5.0)).min(self.max_step_size);
        }

        1
    }
}
