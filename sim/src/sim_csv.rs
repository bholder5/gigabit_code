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
    rate_des_roll: f64,
    rate_des_pitch: f64,
    rate_des_yaw: f64,
    omega_roll: f64,
    omega_pitch: f64,
    omega_yaw: f64,
    omega_gmm_i_roll: f64,
    omega_gmm_i_pitch: f64,
    omega_gmm_i_yaw: f64,
    az: f64,
    el: f64,
    ir: f64,
    az_d: f64,
    el_d: f64,
    ir_d: f64,
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
    // estimation
    ra_hat: f64,
    dec_hat: f64,
    fr_hat: f64,
    gyro_bs_bias_roll: f64,
    gyro_bs_bias_pitch: f64,
    gyro_bs_bias_yaw: f64,
    gyro_bs_bias_roll_act: f64,
    gyro_bs_bias_pitch_act: f64,
    gyro_bs_bias_yaw_act: f64,
}

pub fn push_record(t: &f64, bp: &Params, est: &Estimator, ctrl: &Ctrl, meas: &Meas, sim_st: &State, flex: &fx::Flex_model) -> Result<(), Box<dyn Error>> {
    // let file_path = "/home/brad/sim_data/out.csv";
    let file_path = "/home/brad/data/out.csv";

    if Path::new(&file_path.clone()).exists() {
        let file = OpenOptions::new()
            .append(true)
            .open(&file_path.clone())
            .unwrap();
        let mut wtr = WriterBuilder::new().has_headers(false).from_writer(file);

        wtr.serialize(Record {
            t: t.clone(),
            utc: meas.gps.utc.timestamp() as f64,
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
            az: sim_st.hor.az,
            el: sim_st.hor.el,
            ir: sim_st.hor.ir,
            az_d: ctrl.state.hor_d.az,
            el_d: ctrl.state.hor_d.el,
            ir_d: ctrl.state.hor_d.ir,
            ra: sim_st.eq_k.ra,
            dec: sim_st.eq_k.dec,
            fr: sim_st.eq_k.fr,
            ra_d: ctrl.state.eq_d.ra,
            dec_d: ctrl.state.eq_d.dec,
            fr_d: ctrl.state.eq_d.fr,
            //
            roll: meas.roll,
            pitch: meas.pitch,
            yaw: meas.yaw_p,
            roll_d: ctrl.state.gmb_d.roll,
            pitch_d: ctrl.state.gmb_d.pitch,
            yaw_d: ctrl.state.gmb_d.yaw,
            //
            tau_roll: ctrl.fmot_roll.tau_applied,
            tau_pitch: ctrl.fmot_pitch.tau_applied,
            tau_yaw: -ctrl.rw.tau_applied,
            //Estimation
            ra_hat: est.eq_hat_k.ra,
            dec_hat: est.eq_hat_k.dec,
            fr_hat: est.eq_hat_k.fr,
            gyro_bs_bias_roll: est.gyros_bs.bias[0],
            gyro_bs_bias_pitch: est.gyros_bs.bias[1],
            gyro_bs_bias_yaw: est.gyros_bs.bias[2],
            gyro_bs_bias_roll_act: meas.gyros_bs.bias[0],
            gyro_bs_bias_pitch_act: meas.gyros_bs.bias[1],
            gyro_bs_bias_yaw_act: meas.gyros_bs.bias[2],
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
            rate_des_roll: ctrl.error.rate_des[0],
            rate_des_pitch: ctrl.error.rate_des[1],
            rate_des_yaw: ctrl.error.rate_des[2],
            omega_roll: ctrl.state.omega[0],
            omega_pitch: ctrl.state.omega[1],
            omega_yaw: ctrl.state.omega[2],
            omega_gmm_i_roll: ctrl.error._d_theta[0],
            omega_gmm_i_pitch: ctrl.error._d_theta[1],
            omega_gmm_i_yaw: ctrl.error._d_theta[2],
            az: sim_st.hor.az,
            el: sim_st.hor.el,
            ir: sim_st.hor.ir,
            az_d: ctrl.state.hor_d.az,
            el_d: ctrl.state.hor_d.el,
            ir_d: ctrl.state.hor_d.ir,
            ra: sim_st.eq_k.ra,
            dec: sim_st.eq_k.dec,
            fr: sim_st.eq_k.fr,
            ra_d: ctrl.state.eq_d.ra,
            dec_d: ctrl.state.eq_d.dec,
            fr_d: ctrl.state.eq_d.fr,
            //
            roll: meas.roll,
            pitch: meas.pitch,
            yaw: meas.yaw_p,
            roll_d: ctrl.state.gmb_d.roll,
            pitch_d: ctrl.state.gmb_d.pitch,
            yaw_d: ctrl.state.gmb_d.yaw,
            //
            tau_roll: ctrl.fmot_roll.tau_applied,
            tau_pitch: ctrl.fmot_pitch.tau_applied,
            tau_yaw: ctrl.rw.tau_applied,
            // Estimation
            ra_hat: est.eq_hat_k.ra,
            dec_hat: est.eq_hat_k.dec,
            fr_hat: est.eq_hat_k.fr,
            gyro_bs_bias_roll: est.gyros_bs.bias[0],
            gyro_bs_bias_pitch: est.gyros_bs.bias[1],
            gyro_bs_bias_yaw: est.gyros_bs.bias[2],
            gyro_bs_bias_roll_act: meas.gyros_bs.bias[0],
            gyro_bs_bias_pitch_act: meas.gyros_bs.bias[1],
            gyro_bs_bias_yaw_act: meas.gyros_bs.bias[2],
        })?;
        wtr.flush()?;
    }
    Ok(())
}

#[allow(dead_code)]
pub fn read_last_state(mut _t: f64, bp: &mut Params, est: &mut Estimator, ctrl: &mut Ctrl, meas: &mut Meas, sim_st: &mut State, flex: &fx::Flex_model) -> () {
    let file_path = "/home/brad/data/out.csv";

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
            meas.roll = rec.roll;
            meas.pitch = rec.pitch;
            meas.yaw_p = rec.yaw;
            ctrl.state.gmb_d.roll = rec.roll_d;
            ctrl.state.gmb_d.pitch = rec.pitch_d;
            ctrl.state.gmb_d.yaw = rec.yaw_d;
            ctrl.fmot_roll.tau_applied = rec.tau_roll;
            ctrl.fmot_pitch.tau_applied = rec.tau_pitch;
            ctrl.rw.tau_applied = rec.tau_yaw;
            sim_st.hor.az = rec.az;
            sim_st.hor.el = rec.el;
            sim_st.hor.ir = rec.ir;
            ctrl.state.hor_d.az = rec.az_d;
            ctrl.state.hor_d.el = rec.el_d;
            ctrl.state.hor_d.ir = rec.ir_d;
            est.eq_hat_k.ra = rec.ra_hat;
            est.eq_hat_k.dec = rec.dec_hat;
            est.eq_hat_k.fr = rec.fr_hat;
            est.gyros_bs.bias[0] = rec.gyro_bs_bias_roll; 
            est.gyros_bs.bias[1] = rec.gyro_bs_bias_pitch;
            est.gyros_bs.bias[2] = rec.gyro_bs_bias_yaw;
            meas.gyros_bs.bias[0] = rec.gyro_bs_bias_roll_act;
            meas.gyros_bs.bias[1] = rec.gyro_bs_bias_pitch_act;
            meas.gyros_bs.bias[2] = rec.gyro_bs_bias_yaw_act;
            info!(
                "Previous bp.x read in and being using in simulation {:?}",
                rec
            );
        } else {
            info!("CSV file not availabe/no state to read");
        }
    }
}
