//! This module assembles a struct of structs,
//!  containing one for each type of measurement 
//! generation necessary, ex gyros, star cameras etc
pub mod gyros;
use crate::initialization::Params;

pub struct Meas {
    pub gyros_bs: gyros::Gyro_bs,
}

impl Meas {
    pub fn new() -> Meas {
        Meas {
            gyros_bs: gyros::Gyro_bs::new(),
        }
    }

    pub fn read_measurements(&mut self, bp: &Params) {
        // bp.omega_m was the "measurement" before this function, 
        // hence it is actually omega_k
        self.gyros_bs.omega_k = bp.omega_m;
        self.gyros_bs.generate_measurement();
    }
}