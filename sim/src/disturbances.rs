use rand_distr::{Normal, Distribution};
use nalgebra as na;
use rand;

pub struct Wind_torque {
    pub tau: na::Vector3<f64>,
    pub az: Normal<f64>,
    pub el: Normal<f64>,
    pub roll: Normal<f64>,
    pub en: bool,
}

impl Wind_torque {
    pub fn new() -> Wind_torque {
        Wind_torque {
            tau: na::Vector3::new(0.0, 0.0, 0.0),
            az: Normal::<f64>::new(0.0, 0.034).unwrap(),
            el: Normal::<f64>::new(0.0, 2.83).unwrap(),
            roll: Normal::<f64>::new(0.0, 2.83).unwrap(),
            en: false,
        }
    }

    pub fn generate(&mut self) {
        if self.en{
            self.tau[0] = self.az.sample(&mut rand::thread_rng());
            self.tau[1] = self.roll.sample(&mut rand::thread_rng());
            self.tau[2] = self.el.sample(&mut rand::thread_rng());
        }
    }
}
