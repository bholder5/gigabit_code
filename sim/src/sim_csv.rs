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

use crate::flex_sim::Flex_model;
use crate::initialization::Params;
use adcs::control::Ctrl;
use adcs::control::state::State;
use adcs::estimation::Estimator;
use crate::measurements::Meas;
use crate::flex_sim as fx;
use adcs::control::flex_control as fc;

use std::error::Error;
use csv::WriterBuilder;
pub use std::env;

pub use std::ffi::OsString;
pub use std::fs::{File, OpenOptions};
pub use std::path::Path;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use serde::Deserialize;
use serde::Serialize;
use chrono::{Utc, NaiveDateTime};

#[derive(Debug, Serialize, Deserialize)]
pub struct Record {
    t: f64,
    utc: f64,
    rate_des_roll: f64,
    rate_des_pitch: f64,
    rate_des_yaw: f64,
    omega_roll: f64,
    omega_pitch: f64,
    omega_yaw: f64,
    omega_gmm_i_roll: f64,
    omega_gmm_i_pitch: f64,
    omega_gmm_i_yaw: f64,
    ra: f64,
    dec: f64,
    fr: f64,
    ra_d: f64,
    dec_d: f64,
    fr_d: f64,
    
}

pub fn push_record(t: &f64, bp: &Params, est: &Estimator, ctrl: &Ctrl, meas: &Meas, sim_st: &State, flex: &fx::Flex_model, flex_c: &fc::DampingControl) -> Result<(), Box<dyn Error>> {
    // let file_path = "/home/b/sim_data/out.csv";
    let file_path = "/home/bholder/data/out.csv";
    // let file_path = "/home/brad/data/out.csv";

    if Path::new(&file_path.clone()).exists() {
        let file = OpenOptions::new()
            .append(true)
            .open(&file_path.clone())
            .unwrap();
        let mut wtr = WriterBuilder::new().has_headers(false).from_writer(file);

        wtr.serialize(Record {
            t: t.clone(),
            utc: meas.gps.utc.timestamp() as f64,
            rate_des_roll: ctrl.error.rate_des[0],
            rate_des_pitch: ctrl.error.rate_des[1],
            rate_des_yaw: ctrl.error.rate_des[2],
            omega_roll: ctrl.state.omega[0],
            omega_pitch: ctrl.state.omega[1],
            omega_yaw: ctrl.state.omega[2],
            omega_gmm_i_roll: ctrl.error._d_theta[0],
            omega_gmm_i_pitch: ctrl.error._d_theta[1],
            omega_gmm_i_yaw: ctrl.error._d_theta[2],
            //
            ra: sim_st.eq_k.ra,
            dec: sim_st.eq_k.dec,
            fr: sim_st.eq_k.fr,
            ra_d: ctrl.state.eq_d.ra,
            dec_d: ctrl.state.eq_d.dec,
            fr_d: ctrl.state.eq_d.fr,
        })?;
        wtr.flush()?;
    } else {
        let file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(&file_path.clone())
            .unwrap();
        let mut wtr = WriterBuilder::new().has_headers(true).from_writer(file);
        wtr.serialize(Record {
            t: t.clone(),
            utc: meas.gps.utc.timestamp() as f64,
            rate_des_roll: ctrl.error.rate_des[0],
            rate_des_pitch: ctrl.error.rate_des[1],
            rate_des_yaw: ctrl.error.rate_des[2],
            omega_roll: ctrl.state.omega[0],
            omega_pitch: ctrl.state.omega[1],
            omega_yaw: ctrl.state.omega[2],
            omega_gmm_i_roll: ctrl.error._d_theta[0],
            omega_gmm_i_pitch: ctrl.error._d_theta[1],
            omega_gmm_i_yaw: ctrl.error._d_theta[2],
            ra: sim_st.eq_k.ra,
            dec: sim_st.eq_k.dec,
            fr: sim_st.eq_k.fr,
            ra_d: ctrl.state.eq_d.ra,
            dec_d: ctrl.state.eq_d.dec,
            fr_d: ctrl.state.eq_d.fr,
        })?;
        wtr.flush()?;
    }
    Ok(())
}

#[allow(dead_code)]
pub fn read_last_state(mut _t: f64, bp: &mut Params, est: &mut Estimator, ctrl: &mut Ctrl, meas: &mut Meas, sim_st: &mut State, flex: &fx::Flex_model) -> () {
    let file_path = "/home/bholder/data/out.csv";
    // let file_path = "/home/brad/data/out.csv";

    if Path::new(&file_path.clone()).exists() {
        let file = OpenOptions::new()
            .read(true)
            .open(&file_path.clone())
            .unwrap();

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(file);
        let iter = rdr.deserialize();

        if let Some(result) = iter.last() {
            let rec: Record = result.unwrap();
            println!("{:?}", rec);
            _t = rec.t.clone();
            meas.gps.utc = chrono::DateTime::<Utc>::from_utc(NaiveDateTime::from_timestamp_opt(rec.utc as i64, 0).unwrap(), Utc);
            ctrl.error.rate_des[0] = rec.rate_des_roll;
            ctrl.error.rate_des[1] = rec.rate_des_pitch;
            ctrl.error.rate_des[2] = rec.rate_des_yaw;
            ctrl.state.omega[0] = rec.omega_roll;//bore
            ctrl.state.omega[1] = rec.omega_pitch;//pitch
            ctrl.state.omega[2] = rec.omega_yaw;//cross_pitch
            ctrl.error._d_theta[0] = rec.omega_gmm_i_roll;
            ctrl.error._d_theta[1] = rec.omega_gmm_i_pitch;
            ctrl.error._d_theta[2] = rec.omega_gmm_i_yaw;
            sim_st.eq_k.ra = rec.ra;
            sim_st.eq_k.dec = rec.dec;
            sim_st.eq_k.fr = rec.fr;
            ctrl.state.eq_d.ra = rec.ra_d;
            ctrl.state.eq_d.dec = rec.dec_d;
            ctrl.state.eq_d.fr = rec.fr_d;
            info!(
                "Previous bp.x read in and being using in simulation {:?}",
                rec
            );
        } else {
            info!("CSV file not availabe/no state to read");
        }
    }
}
