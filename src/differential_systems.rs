pub trait DifferentialSystem {
    fn evaluate(&self, x: &[f64], dxdt: &mut [f64], t: f64) -> ();
}

pub struct Attractor {
    pub x_position: f64,
    pub y_position: f64,
    pub force_coefficient: f64
}

pub struct PendulumSystem {
    pub height: f64,
    pub mass: f64,
    pub gravity: f64,
    pub drag: f64,
    pub length: f64,
    pub attractors: Vec<Attractor>
}

impl DifferentialSystem for PendulumSystem {
    #[inline]
    fn evaluate(&self, x: &[f64], dxdt: &mut [f64], t: f64) -> () {
        unsafe {
            let _ = t; // system not dependent on t
            let x_squared = *x.get_unchecked(0) * *x.get_unchecked(0);
            let y_squared = *x.get_unchecked(1) * *x.get_unchecked(1);
            let length_squared = self.length * self.length;
            let norm_squared = x_squared + y_squared;
            let sqrt_term = (1.0 - norm_squared / length_squared).sqrt();
            let gravity_value = -self.gravity / self.length * sqrt_term;

            let mut x_attraction_force = 0.0;
            let mut y_attraction_force = 0.0;

            let value1 = self.height + self.length * (1.0 - sqrt_term);
            let value2 = value1 * value1;

            for attractor in &self.attractors {
                let value3 = *x.get_unchecked(0) - attractor.x_position;
                let value4 = *x.get_unchecked(1) - attractor.y_position;
                let value5 = value3 * value3;
                let value6 = value4 * value4;
                let value7 = -attractor.force_coefficient / (value5 + value6 + value2).powf(1.5);
                x_attraction_force += value3 * value7;
                y_attraction_force += value4 * value7;
            }

            //let value3_0 = *x.get_unchecked(0) - self.attractors.get_unchecked(0).x_position;
            //let value4_0 = *x.get_unchecked(1) - self.attractors.get_unchecked(0).y_position;
            //let value5_0 = value3_0 * value3_0;
            //let value6_0 = value4_0 * value4_0;
            //let value7_0 = -self.attractors.get_unchecked(0).force_coefficient / (value5_0 + value6_0 + value2).powf(1.5);

            //let value3_1 = *x.get_unchecked(0) - self.attractors.get_unchecked(1).x_position;
            //let value4_1 = *x.get_unchecked(1) - self.attractors.get_unchecked(1).y_position;
            //let value5_1 = value3_1 * value3_1;
            //let value6_1 = value4_1 * value4_1;
            //let value7_1 = -self.attractors.get_unchecked(1).force_coefficient / (value5_1 + value6_1 + value2).powf(1.5);

            //let value3_2 = *x.get_unchecked(0) - self.attractors.get_unchecked(2).x_position;
            //let value4_2 = *x.get_unchecked(1) - self.attractors.get_unchecked(2).y_position;
            //let value5_2 = value3_2 * value3_2;
            //let value6_2 = value4_2 * value4_2;
            //let value7_2 = -self.attractors.get_unchecked(2).force_coefficient / (value5_2 + value6_2 + value2).powf(1.5);

            //let x_attraction_force = value3_0 * value7_0 + value3_1 * value7_1 + value3_2 * value7_2;
            //let y_attraction_force = value4_0 * value7_0 + value4_1 * value7_1 + value4_2 * value7_2;

            *dxdt.get_unchecked_mut(0) = *x.get_unchecked(2);
            *dxdt.get_unchecked_mut(1) = *x.get_unchecked(3);
            *dxdt.get_unchecked_mut(2) = (*x.get_unchecked(0) * gravity_value) + (-self.drag * *x.get_unchecked(2) + x_attraction_force) / self.mass;
            *dxdt.get_unchecked_mut(3) = (*x.get_unchecked(1) * gravity_value) + (-self.drag * *x.get_unchecked(3) + y_attraction_force) / self.mass;
        }
    }
}

