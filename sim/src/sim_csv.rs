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

use crate::initialization::Params;
use adcs::control::Ctrl;
use adcs::estimation::Est;
// import crate for matrix math

// use std::fmt;
use std::error::Error;
// use std::io;
use csv::WriterBuilder;
pub use std::env;

pub use std::ffi::OsString;
pub use std::fs::{File, OpenOptions};
use std::io::BufReader;
pub use std::path::Path;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use serde::Deserialize;
use serde::Serialize;

#[derive(Debug, Serialize, Deserialize)]
pub struct Record {
    t: f64,
    dth1: f64,
    dth2: f64,
    dth3: f64,
    dth4: f64,
    dth5: f64,
    dth6: f64,
    dth7: f64,
    dth8: f64,
    dth9: f64,
    th1: f64,
    th2: f64,
    th3: f64,
    th4: f64,
    th5: f64,
    th6: f64,
    th7: f64,
    th8: f64,
    th9: f64,
    h1: f64,
    h2: f64,
    h3: f64,
    // Ctrl
    err_roll: f64,
    err_pitch: f64,
    err_yaw: f64,
    err_fine_yaw: f64,
    err_fine_pitch: f64,
    err_fine_roll: f64,
    err_roll_comb: f64,
    err_pitch_comb: f64,
    err_yaw_comb: f64,
    err_rate_roll: f64,
    err_rate_pitch: f64,
    err_rate_yaw: f64,
    err_rate_sum_roll: f64,
    err_rate_sum_pitch: f64,
    err_rate_sum_yaw: f64,
    omega_roll: f64,
    omega_pitch: f64,
    omega_yaw: f64,
    ra: f64,
    dec: f64,
    fr: f64,
    ra_d: f64,
    dec_d: f64,
    fr_d: f64,
    roll: f64,
    pitch: f64,
    yaw: f64,
    roll_d: f64,
    pitch_d: f64,
    yaw_d: f64,
    tau_roll: f64,
    tau_pitch: f64,
    tau_yaw: f64,
}

pub fn push_record(t: &f64, bp: &Params, _est: &Est, ctrl: &Ctrl) -> Result<(), Box<dyn Error>> {
    // let file_path = "/home/brad/sim_data/out.csv";
    let file_path = "/media/brad/linux_storage/sim_data/out.csv";

    if Path::new(&file_path.clone()).exists() {
        let file = OpenOptions::new()
            .append(true)
            .open(&file_path.clone())
            .unwrap();
        let mut wtr = WriterBuilder::new().has_headers(false).from_writer(file);

        wtr.serialize(Record {
            t: t.clone(),
            dth1: bp.x[0],
            dth2: bp.x[1],
            dth3: bp.x[2],
            dth4: bp.x[3],
            dth5: bp.x[4],
            dth6: bp.x[5],
            dth7: bp.x[6],
            dth8: bp.x[7],
            dth9: bp.x[8],
            th1: bp.x[9],
            th2: bp.x[10],
            th3: bp.x[11],
            th4: bp.x[12],
            th5: bp.x[13],
            th6: bp.x[14],
            th7: bp.x[15],
            th8: bp.x[16],
            th9: bp.x[17],
            h1: bp.x[18],
            h2: bp.x[19],
            h3: bp.x[20],
            // ctrl
            err_roll_comb: ctrl.error.err_comb_th[0],
            err_pitch_comb: ctrl.error.err_comb_th[1],
            err_yaw_comb: ctrl.error.err_comb_th[2],
            err_roll: ctrl.error.err_gmb_th.roll,
            err_pitch: ctrl.error.err_gmb_th.pitch,
            err_yaw: ctrl.error.err_gmb_th.yaw,
            err_fine_yaw: ctrl.error.err_b_th.yaw,
            err_fine_pitch: ctrl.error.err_b_th.pitch,
            err_fine_roll: ctrl.error.err_b_th.roll,
            err_rate_roll: ctrl.error.err_rate[0],
            err_rate_pitch: ctrl.error.err_rate[1],
            err_rate_yaw: ctrl.error.err_rate[2],
            err_rate_sum_roll: ctrl.error.err_rate_sum[0],
            err_rate_sum_pitch: ctrl.error.err_rate_sum[1],
            err_rate_sum_yaw: ctrl.error.err_rate_sum[2],
            omega_roll: ctrl.state.omega[0],
            omega_pitch: ctrl.state.omega[1],
            omega_yaw: ctrl.state.omega[2],
            //
            ra: ctrl.state.eq_k.ra,
            dec: ctrl.state.eq_k.dec,
            fr: ctrl.state.eq_k.fr,
            ra_d: ctrl.state.eq_d.ra,
            dec_d: ctrl.state.eq_d.dec,
            fr_d: ctrl.state.eq_d.fr,
            //
            roll: ctrl.state.gmb_k.roll,
            pitch: ctrl.state.gmb_k.pitch,
            yaw: ctrl.state.gmb_k.yaw,
            roll_d: ctrl.state.gmb_d.roll,
            pitch_d: ctrl.state.gmb_d.pitch,
            yaw_d: ctrl.state.gmb_d.yaw,
            //
            tau_roll: ctrl.fmot_roll.tau_applied,
            tau_pitch: ctrl.fmot_pitch.tau_applied,
            tau_yaw: -ctrl.rw.tau_applied,
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
            dth1: bp.x[0],
            dth2: bp.x[1],
            dth3: bp.x[2],
            dth4: bp.x[3],
            dth5: bp.x[4],
            dth6: bp.x[5],
            dth7: bp.x[6],
            dth8: bp.x[7],
            dth9: bp.x[8],
            th1: bp.x[9],
            th2: bp.x[10],
            th3: bp.x[11],
            th4: bp.x[12],
            th5: bp.x[13],
            th6: bp.x[14],
            th7: bp.x[15],
            th8: bp.x[16],
            th9: bp.x[17],
            h1: bp.x[18],
            h2: bp.x[19],
            h3: bp.x[20],
            // ctrl
            err_roll_comb: ctrl.error.err_comb_th[0],
            err_pitch_comb: ctrl.error.err_comb_th[1],
            err_yaw_comb: ctrl.error.err_comb_th[2],
            err_roll: ctrl.error.err_gmb_th.roll,
            err_pitch: ctrl.error.err_gmb_th.pitch,
            err_yaw: ctrl.error.err_gmb_th.yaw,
            err_fine_yaw: ctrl.error.err_b_th.yaw,
            err_fine_pitch: ctrl.error.err_b_th.pitch,
            err_fine_roll: ctrl.error.err_b_th.roll,
            err_rate_roll: ctrl.error.err_rate[0],
            err_rate_pitch: ctrl.error.err_rate[1],
            err_rate_yaw: ctrl.error.err_rate[2],
            err_rate_sum_roll: ctrl.error.err_rate_sum[0],
            err_rate_sum_pitch: ctrl.error.err_rate_sum[1],
            err_rate_sum_yaw: ctrl.error.err_rate_sum[2],
            omega_roll: ctrl.state.omega[0],
            omega_pitch: ctrl.state.omega[1],
            omega_yaw: ctrl.state.omega[2],
            ra: ctrl.state.eq_k.ra,
            dec: ctrl.state.eq_k.dec,
            fr: ctrl.state.eq_k.fr,
            ra_d: ctrl.state.eq_d.ra,
            dec_d: ctrl.state.eq_d.dec,
            fr_d: ctrl.state.eq_d.fr,
            //
            roll: ctrl.state.gmb_k.roll,
            pitch: ctrl.state.gmb_k.pitch,
            yaw: ctrl.state.gmb_k.yaw,
            roll_d: ctrl.state.gmb_d.roll,
            pitch_d: ctrl.state.gmb_d.pitch,
            yaw_d: ctrl.state.gmb_d.yaw,
            //
            tau_roll: ctrl.fmot_roll.tau_applied,
            tau_pitch: ctrl.fmot_pitch.tau_applied,
            tau_yaw: ctrl.rw.tau_applied,
        })?;
        wtr.flush()?;
    }
    Ok(())
}

#[allow(dead_code)]
pub fn read_last_state(mut t: f64, bp: &mut Params, _est: &Est, ctrl: &mut Ctrl) -> () {
    let file_path = "/media/brad/linux_storage/sim_data/out.csv";

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
            t = rec.t.clone();
            bp.x[0] = rec.dth1;
            bp.x[1] = rec.dth2;
            bp.x[2] = rec.dth3;
            bp.x[3] = rec.dth4;
            bp.x[4] = rec.dth5;
            bp.x[5] = rec.dth6;
            bp.x[6] = rec.dth7;
            bp.x[7] = rec.dth8;
            bp.x[8] = rec.dth9;
            bp.x[9] = rec.th1;
            bp.x[10] = rec.th2;
            bp.x[11] = rec.th3;
            bp.x[12] = rec.th4;
            bp.x[13] = rec.th5;
            bp.x[14] = rec.th6;
            bp.x[15] = rec.th7;
            bp.x[16] = rec.th8;
            bp.x[17] = rec.th9;
            bp.x[18] = rec.h1;
            bp.x[19] = rec.h2;
            bp.x[20] = rec.h3;
            // ctrl
            ctrl.error.err_comb_th[0] = rec.err_roll_comb;
            ctrl.error.err_comb_th[1] = rec.err_pitch_comb;
            ctrl.error.err_comb_th[2] = rec.err_yaw_comb;
            ctrl.error.err_gmb_th.roll = rec.err_roll;
            ctrl.error.err_gmb_th.pitch = rec.err_pitch;
            ctrl.error.err_gmb_th.yaw = rec.err_yaw;
            ctrl.error.err_b_th.yaw = rec.err_fine_yaw;
            ctrl.error.err_b_th.pitch = rec.err_fine_pitch;
            ctrl.error.err_b_th.roll = rec.err_fine_roll;
            ctrl.error.err_rate[0] = rec.err_rate_roll;
            ctrl.error.err_rate[1] = rec.err_rate_pitch;
            ctrl.error.err_rate[2] = rec.err_rate_yaw;
            ctrl.error.err_rate_sum[0] = rec.err_rate_sum_roll;
            ctrl.error.err_rate_sum[1] = rec.err_rate_sum_pitch;
            ctrl.error.err_rate_sum[2] = rec.err_rate_sum_yaw;
            ctrl.state.omega[0] = rec.omega_roll;
            ctrl.state.omega[1] = rec.omega_pitch;
            ctrl.state.omega[2] = rec.omega_yaw;
            ctrl.state.eq_k.ra = rec.ra;
            ctrl.state.eq_k.dec = rec.dec;
            ctrl.state.eq_k.fr = rec.fr;
            ctrl.state.eq_d.ra = rec.ra_d;
            ctrl.state.eq_d.dec = rec.dec_d;
            ctrl.state.eq_d.fr = rec.fr_d;
            ctrl.state.gmb_k.roll = rec.roll;
            ctrl.state.gmb_k.pitch = rec.pitch;
            ctrl.state.gmb_k.yaw = rec.yaw;
            ctrl.state.gmb_d.roll = rec.roll_d;
            ctrl.state.gmb_d.pitch = rec.pitch_d;
            ctrl.state.gmb_d.yaw = rec.yaw_d;
            ctrl.fmot_roll.tau_applied = rec.tau_roll;
            ctrl.fmot_pitch.tau_applied = rec.tau_pitch;
            ctrl.rw.tau_applied = rec.tau_yaw;
            info!(
                "Previous bp.x read in and being using in simulation {:?}",
                rec
            );
        } else {
            info!("CSV file not availabe/no state to read");
        }
    }
}
