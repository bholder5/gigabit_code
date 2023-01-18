//! This module assembles a struct of structs,
//!  containing one for each type of measurement 
//! generation necessary, ex gyros, star cameras etc
pub mod gyros;
use crate::initialization::Params;
use adcs::control::state as st;
use crate::flex_sim;
use nalgebra as na;
#[derive(Clone)]
pub struct Meas {
    pub gyros_bs: gyros::Gyro_bs,
    pub gps: st::gps::Gps,
    pub cbh: na::Rotation3<f64>,
    pub roll: f64,
    pub pitch: f64,
    /// yaw_p is azimuth with opposite sign, not the position of OF WRT pivot
    pub yaw_p: f64,
}

impl Meas {
    pub fn new() -> Meas {
        Meas {
            gyros_bs: gyros::Gyro_bs::new(),
            gps: st::gps::Gps::new(),
            cbh: na::Rotation3::<f64>::identity(),
            roll: 0.0,
            pitch: 0.0,
            yaw_p: 0.0,
        }
    }

    pub fn read_measurements(&mut self, bp: &Params, sim_state: &st::State, flex: &flex_sim::Flex_model) -> () {
        // READ GYROS/GENERATE GYRO BORESIGHT MEASUREMENT
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
        // need to add in flexible affects.
        self.gyros_bs.omega_k = bp.omega_m + (self.cbh * ceh.transpose() * w_di);

        if flex.flex_enable{
            self.gyros_bs.omega_k = flex.g1_out + self.gyros_bs.omega_k;
        }

        self.gyros_bs.generate_measurement();

        // ENCODERS

        if flex.flex_enable{
            self.roll = bp.x[16] + flex.c_out[1];
            self.pitch = bp.x[17] + flex.c_out[3];
        } else {
            self.roll = bp.x[16];
            self.pitch = bp.x[17];
        }
        
        self.yaw_p = -sim_state.hor.az;

    }
}