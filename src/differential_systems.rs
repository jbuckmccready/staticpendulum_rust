pub trait DifferentialSystem<const D: usize> {
    fn evaluate(&self, x: &[f64; D], dxdt: &mut [f64; D], t: f64);
}

pub struct Attractor {
    pub x_position: f64,
    pub y_position: f64,
    pub force_coefficient: f64,
}

pub struct PendulumSystem {
    pub height: f64,
    pub mass: f64,
    pub gravity: f64,
    pub drag: f64,
    pub length: f64,
    pub attractors: Vec<Attractor>,
    length_squared: f64,
    gravity_over_length: f64,
}

impl PendulumSystem {
    pub fn new(
        height: f64,
        mass: f64,
        gravity: f64,
        drag: f64,
        length: f64,
        attractors: Vec<Attractor>,
    ) -> PendulumSystem {
        PendulumSystem {
            height,
            mass,
            gravity,
            drag,
            length,
            attractors,
            length_squared: length * length,
            gravity_over_length: gravity / length,
        }
    }
}

impl DifferentialSystem<4> for PendulumSystem {
    fn evaluate(&self, x: &[f64; 4], dxdt: &mut [f64; 4], t: f64) {
        let _ = t; // system not dependent on t
        let x_squared = x[0] * x[0];
        let y_squared = x[1] * x[1];
        let norm_squared = x_squared + y_squared;
        let sqrt_term = (1.0 - norm_squared / self.length_squared).sqrt();
        let gravity_value = -self.gravity_over_length * sqrt_term;

        let mut x_attraction_force = 0.0;
        let mut y_attraction_force = 0.0;

        let value1 = self.height + self.length * (1.0 - sqrt_term);
        let value2 = value1 * value1;

        for attractor in &self.attractors {
            let value3 = x[0] - attractor.x_position;
            let value4 = x[1] - attractor.y_position;
            let value5 = value3 * value3;
            let value6 = value4 * value4;
            let value7 = -attractor.force_coefficient / (value5 + value6 + value2).powf(1.5);
            x_attraction_force += value3 * value7;
            y_attraction_force += value4 * value7;
        }

        dxdt[0] = x[2];
        dxdt[1] = x[3];
        dxdt[2] = x[0] * gravity_value + (-self.drag * x[2] + x_attraction_force) / self.mass;
        dxdt[3] = x[1] * gravity_value + (-self.drag * x[3] + y_attraction_force) / self.mass;
    }
}
