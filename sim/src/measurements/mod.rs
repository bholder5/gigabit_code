//! This module assembles a struct of structs,
//!  containing one for each type of measurement 
//! generation necessary, ex gyros, star cameras etc
pub mod gyros;
use crate::initialization::Params;
use adcs::control::state as st;
use nalgebra as na;

pub struct Meas {
    pub gyros_bs: gyros::Gyro_bs,
    pub gps: st::gps::Gps,
    pub cbh: na::Rotation3<f64>
}

impl Meas {
    pub fn new() -> Meas {
        Meas {
            gyros_bs: gyros::Gyro_bs::new(),
            gps: st::gps::Gps::new(),
            cbh: na::Rotation3::<f64>::identity()
        }
    }

    pub fn read_measurements(&mut self, bp: &Params) {
        // bp.omega_m was the "measurement" before this function, 
        // hence it is actually omega_k
        self.gps.get_greenwhich_apparent_sidereal_time();

        let last = self.gps.gast + self.gps.lon;
        let arg_z = (last + 180.0).to_radians();

        let arg_y = (self.gps.lat - 90.0).to_radians();

        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), arg_y);
        let rotz = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::z_axis(), arg_z);
        let ceh = rotz * roty;

        let w_di = na::Vector3::<f64>::new(0.0, 0.0, 15.04108/206265.0);
        self.gyros_bs.omega_k = bp.omega_m + (self.cbh * ceh.transpose() * w_di);
        self.gyros_bs.generate_measurement();
    }
}