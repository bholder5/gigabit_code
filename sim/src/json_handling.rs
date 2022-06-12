//! This module contains all of the initialization of BITs
//! parameters for use in the simulation
//!
//! provides the main struct formulation that carries all of
//! the parameters, defines matrix and vector types necessary,
//! then a function which is called to then initialize all of
//! the values. Changing BIT's parameters that have to do with the
//! bit simulation parameters happens here
//!
//! Examples would be:
//! -generating sensors measurements from the simulation state
//!

use adcs::control::Ctrl;
use adcs::estimation::Estimator;
use crate::initialization::Params;

// use std::fmt;
pub use std::env;

pub use std::ffi::OsString;
pub use std::fs::{File, OpenOptions};
use std::io::BufReader;
pub use std::path::Path;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use serde::Deserialize;

#[derive(Debug, Deserialize, Default)]
pub struct Gain_json {
    kp_roll: f64,
    kp_pitch: f64,
    kp_yaw: f64,
    ki_roll: f64,
    ki_pitch: f64,
    ki_yaw: f64,
    kip_roll: f64,
    kip_pitch: f64,
    kip_yaw: f64,
}
#[derive(Debug, Deserialize, Default)]
pub struct Read_json {
    fine: Gain_json,
    yaw_des: f64,
    pitch_des: f64,
    roll_des: f64,
    new_targ: bool,
    latency: bool,
    ctrl_from_est: bool,
    est_reset: bool,
}

pub fn read_gains(ctrl: &mut Ctrl, bp: &mut Params, est: &mut Estimator) -> () {
    let path = "sim/src/gains.json";

    #[allow(unused_assignments)]
    let mut u = Read_json::default();

    let _data = match File::open(path) {
        Ok(file) => {
            let reader = BufReader::new(file);
            let _gains = match serde_json::from_reader(reader) {
                Ok(gains) => {
                    u = gains;
                    //rewrite control gain vecs;
                    let kp_vec = [u.fine.kp_roll, u.fine.kp_pitch, u.fine.kp_yaw];
                    let ki_vec = [u.fine.ki_roll, u.fine.ki_pitch, u.fine.ki_yaw];
                    let kip_vec = [u.fine.kip_roll, u.fine.kip_pitch, u.fine.kip_yaw];

                    ctrl.fine_gains
                        .update_gain_matrices(&kp_vec, &kip_vec, &ki_vec);

                    if u.new_targ {
                        ctrl.state.gmb_d.roll = u.roll_des;
                        ctrl.state.gmb_d.pitch = u.pitch_des;
                        ctrl.state.gmb_d.yaw = u.yaw_des;
                        ctrl.state.gmb_d.calculate_rotation_matrix();
                        ctrl.state.update_desired_eq_from_gmb();
                        u.new_targ = false;
                        // println!("\n\n\n READ GAINS UPDATE\n\n\n{:?}", ctrl.state.gmb_d);
                    }

                    bp.latency = u.latency;
                    bp.ctrl_from_est = u.ctrl_from_est;
                    est.reset = u.est_reset;

                }
                Err(err) => {
                    error!("error reading json: {}", err);
                }
            };
        }
        Err(err) => {
            error!("error opening file: {}", err);
        }
    };
}
