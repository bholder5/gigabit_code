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
use adcs::control::servo;

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
    t: f32,
    utc: f32,
    dth1: f32,
    dth2: f32,
    dth3: f32,
    dth4: f32,
    dth5: f32,
    dth6: f32,
    dth7: f32,
    dth8: f32,
    dth9: f32,
    th1: f32,
    th2: f32,
    th3: f32,
    th4: f32,
    th5: f32,
    th6: f32,
    th7: f32,
    th8: f32,
    th9: f32,
    h1: f32,
    h2: f32,
    h3: f32,
    // flex
    flx1: f32,
    flx2: f32,
    flx3: f32,
    flx4: f32,
    flx5: f32,
    flx6: f32,
    flx7: f32,
    flx8: f32,
    flx9: f32,
    flx10: f32,
    // flex measurements
    flx_yaw_r: f32,
    flx_yaw_p: f32,
    flx_bow_r: f32,
    flx_bow_p: f32,
    flx_stern_r: f32,
    flx_stern_p: f32,
    flx_port_r: f32,
    flx_port_p: f32,
    flx_sb_r: f32,
    flx_sb_p: f32,
    // Ctrl
    err_roll: f32,
    err_pitch: f32,
    err_yaw: f32,
    err_fine_yaw: f32,
    err_fine_pitch: f32,
    err_fine_roll: f32,
    err_roll_comb: f32,
    err_pitch_comb: f32,
    err_yaw_comb: f32,
    err_rate_roll: f32,
    err_rate_pitch: f32,
    err_rate_yaw: f32,
    err_rate_sum_roll: f32,
    err_rate_sum_pitch: f32,
    err_rate_sum_yaw: f32,
    rate_des_roll: f32,
    rate_des_pitch: f32,
    rate_des_yaw: f32,
    omega_roll: f32,
    omega_pitch: f32,
    omega_yaw: f32,
    omega_gmm_i_roll: f32,
    omega_gmm_i_pitch: f32,
    omega_gmm_i_yaw: f32,
    az: f32,
    el: f32,
    ir: f32,
    az_d: f32,
    el_d: f32,
    ir_d: f32,
    ra: f32,
    dec: f32,
    fr: f32,
    ra_d: f32,
    dec_d: f32,
    fr_d: f32,
    roll: f32,
    pitch: f32,
    yaw: f32,
    roll_d: f32,
    pitch_d: f32,
    yaw_d: f32,
    tau_roll: f32,
    tau_pitch: f32,
    tau_yaw: f32,
    tau_flex_yaw: f32,
    tau_flex_bow: f32,
    tau_flex_stern: f32,
    tau_flex_port: f32,
    tau_flex_sb: f32,
    // estimation
    ra_hat: f32,
    dec_hat: f32,
    fr_hat: f32,
    gyro_bs_bias_roll: f32,
    gyro_bs_bias_pitch: f32,
    gyro_bs_bias_yaw: f32,
    gyro_bs_bias_roll_act: f32,
    gyro_bs_bias_pitch_act: f32,
    gyro_bs_bias_yaw_act: f32,
    // servo control
    yaw_z1: f32,
    yaw_z2: f32,
    roll_z1: f32,
    roll_z2: f32,
    pitch_z1: f32,
    pitch_z2: f32,
    yaw_torque: f32,
    roll_torque: f32,
    pitch_torque: f32,
}

pub fn push_record(t: &f64, bp: &Params, est: &Estimator, ctrl: &Ctrl, meas: &Meas, sim_st: &State, flex: &fx::Flex_model, flex_c: &fc::DampingControl, servo: &servo::ServoControl) -> Result<(), Box<dyn Error>> {
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
            t: t.clone()as f32,
            utc: meas.gps.utc.timestamp() as f32,
            dth1: bp.x[0]as f32,
            dth2: bp.x[1]as f32,
            dth3: bp.x[2]as f32,
            dth4: bp.x[3]as f32,
            dth5: bp.x[4]as f32,
            dth6: bp.x[5]as f32,
            dth7: bp.x[6]as f32,
            dth8: bp.x[7]as f32,
            dth9: bp.x[8]as f32,
            th1: bp.x[9]as f32,
            th2: bp.x[10]as f32,
            th3: bp.x[11]as f32,
            th4: bp.x[12]as f32,
            th5: bp.x[13]as f32,
            th6: bp.x[14]as f32,
            th7: bp.x[15]as f32,
            th8: bp.x[16]as f32,
            th9: bp.x[17]as f32,
            h1: bp.x[18] as f32,
            h2: bp.x[19] as f32,
            h3: bp.x[20] as f32,
            // flex
            flx1: flex.eta[0] as f32,
            flx2: flex.eta[1] as f32,
            flx3: flex.eta[2] as f32,
            flx4: flex.eta[3] as f32,
            flx5: flex.eta[4] as f32,
            flx6: flex.eta[5] as f32,
            flx7: flex.eta[6] as f32,
            flx8: flex.eta[7] as f32,
            flx9: flex.eta[8] as f32,
            flx10: flex.eta[9] as f32,
            //
            // flex measurements
            flx_yaw_r: flex.c_out[0] as f32,
            flx_yaw_p: flex.c_pos_out[0] as f32,
            flx_bow_r: flex.c_out[1] as f32,
            flx_bow_p: flex.c_pos_out[1] as f32,
            flx_stern_r: flex.c_out[2] as f32,
            flx_stern_p: flex.c_pos_out[2] as f32,
            flx_port_r: flex.c_out[3] as f32,
            flx_port_p: flex.c_pos_out[3] as f32,
            flx_sb_r: 0.0, //flex.c_out[4] as f32,
            flx_sb_p: 0.0, //flex.c_pos_out[4] as f32,
            // ctrl
            err_roll_comb: ctrl.error.err_comb_th[0] as f32,
            err_pitch_comb: ctrl.error.err_comb_th[1] as f32,
            err_yaw_comb: ctrl.error.err_comb_th[2] as f32,
            err_roll: ctrl.error.err_gmb_th.roll as f32,
            err_pitch: ctrl.error.err_gmb_th.pitch as f32,
            err_yaw: ctrl.error.err_gmb_th.yaw as f32,
            err_fine_yaw: ctrl.error.err_b_th.yaw as f32,
            err_fine_pitch: ctrl.error.err_b_th.pitch as f32,
            err_fine_roll: ctrl.error.err_b_th.roll as f32,
            err_rate_roll: ctrl.error.err_rate[0] as f32,
            err_rate_pitch: ctrl.error.err_rate[1] as f32,
            err_rate_yaw: ctrl.error.err_rate[2] as f32,
            err_rate_sum_roll: ctrl.error.err_rate_sum[0] as f32,
            err_rate_sum_pitch: ctrl.error.err_rate_sum[1] as f32,
            err_rate_sum_yaw: ctrl.error.err_rate_sum[2] as f32,
            rate_des_roll: ctrl.error.rate_des[0] as f32,
            rate_des_pitch: ctrl.error.rate_des[1] as f32,
            rate_des_yaw: ctrl.error.rate_des[2] as f32,
            omega_roll: ctrl.state.omega[0] as f32,
            omega_pitch: ctrl.state.omega[1] as f32,
            omega_yaw: ctrl.state.omega[2] as f32,
            omega_gmm_i_roll: ctrl.error._d_theta[0] as f32,
            omega_gmm_i_pitch: ctrl.error._d_theta[1] as f32,
            omega_gmm_i_yaw: ctrl.error._d_theta[2] as f32,
            //
            az: sim_st.hor.az as f32,
            el: sim_st.hor.el as f32,
            ir: sim_st.hor.ir as f32,
            az_d: ctrl.state.hor_d.az as f32,
            el_d: ctrl.state.hor_d.el as f32,
            ir_d: ctrl.state.hor_d.ir as f32,
            ra: sim_st.eq_k.ra as f32,
            dec: sim_st.eq_k.dec as f32,
            fr: sim_st.eq_k.fr as f32,
            ra_d: ctrl.state.eq_d.ra as f32,
            dec_d: ctrl.state.eq_d.dec as f32,
            fr_d: ctrl.state.eq_d.fr as f32,
            //
            roll: meas.roll as f32,
            pitch: meas.pitch as f32,
            yaw: meas.yaw_p as f32,
            roll_d: ctrl.state.gmb_d.roll as f32,
            pitch_d: ctrl.state.gmb_d.pitch as f32,
            yaw_d: ctrl.state.gmb_d.yaw as f32,
            //
            tau_roll: ctrl.fmot_roll.tau_applied as f32,
            tau_pitch: ctrl.fmot_pitch.tau_applied as f32,
            tau_yaw: -ctrl.rw.tau_applied as f32,
            tau_flex_yaw: flex_c.u[0] as f32,
            tau_flex_bow: flex_c.u[1] as f32,
            tau_flex_stern: flex_c.u[2] as f32,
            tau_flex_port: flex_c.u[3] as f32,
            tau_flex_sb: flex_c.u[4] as f32,
            //Estimation
            ra_hat: est.eq_hat_k.ra as f32,
            dec_hat: est.eq_hat_k.dec as f32,
            fr_hat: est.eq_hat_k.fr as f32,
            gyro_bs_bias_roll: est.gyros_bs.bias[0] as f32,
            gyro_bs_bias_pitch: est.gyros_bs.bias[1] as f32,
            gyro_bs_bias_yaw: est.gyros_bs.bias[2] as f32,
            gyro_bs_bias_roll_act: meas.gyros_bs.bias[0] as f32,
            gyro_bs_bias_pitch_act: meas.gyros_bs.bias[1] as f32,
            gyro_bs_bias_yaw_act: meas.gyros_bs.bias[2] as f32,
            // servo control
            yaw_z1: servo.yaw.z[0] as f32,
            yaw_z2: servo.yaw.z[1] as f32,
            roll_z1: servo.roll.z[0] as f32,
            roll_z2: servo.roll.z[1] as f32,
            pitch_z1: servo.pitch.z[0] as f32,
            pitch_z2: servo.pitch.z[1] as f32,
            yaw_torque: servo.tau_y_lim as f32,
            roll_torque: servo.tau_r_lim as f32,
            pitch_torque: servo.tau_p_lim as f32,
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
            t: t.clone() as f32,
            utc: meas.gps.utc.timestamp() as f32,
            dth1: bp.x[0] as f32,
            dth2: bp.x[1] as f32,
            dth3: bp.x[2] as f32,
            dth4: bp.x[3] as f32,
            dth5: bp.x[4] as f32,
            dth6: bp.x[5] as f32,
            dth7: bp.x[6] as f32,
            dth8: bp.x[7] as f32,
            dth9: bp.x[8] as f32,
            th1: bp.x[9] as f32,
            th2: bp.x[10] as f32,
            th3: bp.x[11] as f32,
            th4: bp.x[12] as f32,
            th5: bp.x[13] as f32,
            th6: bp.x[14] as f32,
            th7: bp.x[15] as f32,
            th8: bp.x[16] as f32,
            th9: bp.x[17] as f32,
            h1: bp.x[18] as f32,
            h2: bp.x[19] as f32,
            h3: bp.x[20] as f32,
            //flex
            flx1: flex.eta[0] as f32,
            flx2: flex.eta[1] as f32,
            flx3: flex.eta[2] as f32,
            flx4: flex.eta[3] as f32,
            flx5: flex.eta[4] as f32,
            flx6: flex.eta[5] as f32,
            flx7: flex.eta[6] as f32,
            flx8: flex.eta[7] as f32,
            flx9: flex.eta[8] as f32,
            flx10: flex.eta[9] as f32,
            // flex measurements
            flx_yaw_r: flex.c_out[0] as f32,
            flx_yaw_p: flex.c_pos_out[0] as f32,
            flx_bow_r: flex.c_out[1] as f32,
            flx_bow_p: flex.c_pos_out[1] as f32,
            flx_stern_r: flex.c_out[2] as f32,
            flx_stern_p: flex.c_pos_out[2] as f32,
            flx_port_r: flex.c_out[3] as f32,
            flx_port_p: flex.c_pos_out[3] as f32,
            flx_sb_r: flex.c_out[4] as f32,
            flx_sb_p: flex.c_pos_out[4] as f32,
            // ctrl
            err_roll_comb: ctrl.error.err_comb_th[0] as f32,
            err_pitch_comb: ctrl.error.err_comb_th[1] as f32,
            err_yaw_comb: ctrl.error.err_comb_th[2] as f32,
            err_roll: ctrl.error.err_gmb_th.roll as f32,
            err_pitch: ctrl.error.err_gmb_th.pitch as f32,
            err_yaw: ctrl.error.err_gmb_th.yaw as f32,
            err_fine_yaw: ctrl.error.err_b_th.yaw as f32,
            err_fine_pitch: ctrl.error.err_b_th.pitch as f32,
            err_fine_roll: ctrl.error.err_b_th.roll as f32,
            err_rate_roll: ctrl.error.err_rate[0] as f32,
            err_rate_pitch: ctrl.error.err_rate[1] as f32,
            err_rate_yaw: ctrl.error.err_rate[2] as f32,
            err_rate_sum_roll: ctrl.error.err_rate_sum[0] as f32,
            err_rate_sum_pitch: ctrl.error.err_rate_sum[1] as f32,
            err_rate_sum_yaw: ctrl.error.err_rate_sum[2] as f32,
            rate_des_roll: ctrl.error.rate_des[0] as f32,
            rate_des_pitch: ctrl.error.rate_des[1] as f32,
            rate_des_yaw: ctrl.error.rate_des[2] as f32,
            omega_roll: ctrl.state.omega[0] as f32,
            omega_pitch: ctrl.state.omega[1] as f32,
            omega_yaw: ctrl.state.omega[2] as f32,
            omega_gmm_i_roll: ctrl.error._d_theta[0] as f32,
            omega_gmm_i_pitch: ctrl.error._d_theta[1] as f32,
            omega_gmm_i_yaw: ctrl.error._d_theta[2] as f32,
            az: sim_st.hor.az as f32,
            el: sim_st.hor.el as f32,
            ir: sim_st.hor.ir as f32,
            az_d: ctrl.state.hor_d.az as f32,
            el_d: ctrl.state.hor_d.el as f32,
            ir_d: ctrl.state.hor_d.ir as f32,
            ra: sim_st.eq_k.ra as f32,
            dec: sim_st.eq_k.dec as f32,
            fr: sim_st.eq_k.fr as f32,
            ra_d: ctrl.state.eq_d.ra as f32,
            dec_d: ctrl.state.eq_d.dec as f32,
            fr_d: ctrl.state.eq_d.fr as f32,
            //
            roll: meas.roll as f32,
            pitch: meas.pitch as f32,
            yaw: meas.yaw_p as f32,
            roll_d: ctrl.state.gmb_d.roll as f32,
            pitch_d: ctrl.state.gmb_d.pitch as f32,
            yaw_d: ctrl.state.gmb_d.yaw as f32,
            //
            tau_roll: ctrl.fmot_roll.tau_applied as f32,
            tau_pitch: ctrl.fmot_pitch.tau_applied as f32,
            tau_yaw: ctrl.rw.tau_applied as f32,
            tau_flex_yaw: flex_c.u[0] as f32,
            tau_flex_bow: flex_c.u[1] as f32,
            tau_flex_stern: flex_c.u[2] as f32,
            tau_flex_port: flex_c.u[3] as f32,
            tau_flex_sb: flex_c.u[4] as f32,
            // Estimation
            ra_hat: est.eq_hat_k.ra as f32,
            dec_hat: est.eq_hat_k.dec as f32,
            fr_hat: est.eq_hat_k.fr as f32,
            gyro_bs_bias_roll: est.gyros_bs.bias[0] as f32,
            gyro_bs_bias_pitch: est.gyros_bs.bias[1] as f32,
            gyro_bs_bias_yaw: est.gyros_bs.bias[2] as f32,
            gyro_bs_bias_roll_act: meas.gyros_bs.bias[0] as f32,
            gyro_bs_bias_pitch_act: meas.gyros_bs.bias[1] as f32,
            gyro_bs_bias_yaw_act: meas.gyros_bs.bias[2] as f32,
            // servo control
            yaw_z1: servo.yaw.z[0] as f32,
            yaw_z2: servo.yaw.z[1] as f32,
            roll_z1: servo.roll.z[0] as f32,
            roll_z2: servo.roll.z[1] as f32,
            pitch_z1: servo.pitch.z[0] as f32,
            pitch_z2: servo.pitch.z[1] as f32,
            yaw_torque: servo.tau_y_lim as f32,
            roll_torque: servo.tau_r_lim as f32,
            pitch_torque: servo.tau_p_lim as f32,
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
            _t = rec.t.clone() as f64;
            meas.gps.utc = chrono::DateTime::<Utc>::from_utc(NaiveDateTime::from_timestamp_opt(rec.utc as i64, 0).unwrap(), Utc);
            bp.x[0] = rec.dth1 as f64;
            bp.x[1] = rec.dth2 as f64;
            bp.x[2] = rec.dth3 as f64;
            bp.x[3] = rec.dth4 as f64;
            bp.x[4] = rec.dth5 as f64;
            bp.x[5] = rec.dth6 as f64;
            bp.x[6] = rec.dth7 as f64;
            bp.x[7] = rec.dth8 as f64;
            bp.x[8] = rec.dth9 as f64;
            bp.x[9] = rec.th1 as f64;
            bp.x[10] = rec.th2 as f64;
            bp.x[11] = rec.th3 as f64;
            bp.x[12] = rec.th4 as f64;
            bp.x[13] = rec.th5 as f64;
            bp.x[14] = rec.th6 as f64;
            bp.x[15] = rec.th7 as f64;
            bp.x[16] = rec.th8 as f64;
            bp.x[17] = rec.th9 as f64;
            bp.x[18] = rec.h1 as f64;
            bp.x[19] = rec.h2 as f64;
            bp.x[20] = rec.h3 as f64;
            // ctrl
            ctrl.error.err_comb_th[0] = rec.err_roll_comb as f64;
            ctrl.error.err_comb_th[1] = rec.err_pitch_comb as f64;
            ctrl.error.err_comb_th[2] = rec.err_yaw_comb as f64;
            ctrl.error.err_gmb_th.roll = rec.err_roll as f64;
            ctrl.error.err_gmb_th.pitch = rec.err_pitch as f64;
            ctrl.error.err_gmb_th.yaw = rec.err_yaw as f64;
            ctrl.error.err_b_th.yaw = rec.err_fine_yaw as f64;
            ctrl.error.err_b_th.pitch = rec.err_fine_pitch as f64;
            ctrl.error.err_b_th.roll = rec.err_fine_roll as f64;
            ctrl.error.err_rate[0] = rec.err_rate_roll as f64;
            ctrl.error.err_rate[1] = rec.err_rate_pitch as f64;
            ctrl.error.err_rate[2] = rec.err_rate_yaw as f64;
            ctrl.error.err_rate_sum[0] = rec.err_rate_sum_roll as f64;
            ctrl.error.err_rate_sum[1] = rec.err_rate_sum_pitch as f64;
            ctrl.error.err_rate_sum[2] = rec.err_rate_sum_yaw as f64;
            ctrl.error.rate_des[0] = rec.rate_des_roll as f64;
            ctrl.error.rate_des[1] = rec.rate_des_pitch as f64;
            ctrl.error.rate_des[2] = rec.rate_des_yaw as f64;
            ctrl.state.omega[0] = rec.omega_roll as f64;//bore
            ctrl.state.omega[1] = rec.omega_pitch as f64;//pitch
            ctrl.state.omega[2] = rec.omega_yaw as f64;//cross_pitch
            ctrl.error._d_theta[0] = rec.omega_gmm_i_roll as f64;
            ctrl.error._d_theta[1] = rec.omega_gmm_i_pitch as f64;
            ctrl.error._d_theta[2] = rec.omega_gmm_i_yaw as f64;
            sim_st.eq_k.ra = rec.ra as f64;
            sim_st.eq_k.dec = rec.dec as f64;
            sim_st.eq_k.fr = rec.fr as f64;
            ctrl.state.eq_d.ra = rec.ra_d as f64;
            ctrl.state.eq_d.dec = rec.dec_d as f64;
            ctrl.state.eq_d.fr = rec.fr_d as f64;
            meas.roll = rec.roll as f64;
            meas.pitch = rec.pitch as f64;
            meas.yaw_p = rec.yaw as f64;
            ctrl.state.gmb_d.roll = rec.roll_d as f64;
            ctrl.state.gmb_d.pitch = rec.pitch_d as f64;
            ctrl.state.gmb_d.yaw = rec.yaw_d as f64;
            ctrl.fmot_roll.tau_applied = rec.tau_roll as f64;
            ctrl.fmot_pitch.tau_applied = rec.tau_pitch as f64;
            ctrl.rw.tau_applied = rec.tau_yaw as f64;
            sim_st.hor.az = rec.az as f64;
            sim_st.hor.el = rec.el as f64;
            sim_st.hor.ir = rec.ir as f64;
            ctrl.state.hor_d.az = rec.az_d as f64;
            ctrl.state.hor_d.el = rec.el_d as f64;
            ctrl.state.hor_d.ir = rec.ir_d as f64;
            est.eq_hat_k.ra = rec.ra_hat as f64;
            est.eq_hat_k.dec = rec.dec_hat as f64;
            est.eq_hat_k.fr = rec.fr_hat as f64;
            est.gyros_bs.bias[0] = rec.gyro_bs_bias_roll as f64; 
            est.gyros_bs.bias[1] = rec.gyro_bs_bias_pitch as f64;
            est.gyros_bs.bias[2] = rec.gyro_bs_bias_yaw as f64;
            meas.gyros_bs.bias[0] = rec.gyro_bs_bias_roll_act as f64;
            meas.gyros_bs.bias[1] = rec.gyro_bs_bias_pitch_act as f64;
            meas.gyros_bs.bias[2] = rec.gyro_bs_bias_yaw_act as f64;
            info!(
                "Previous bp.x read in and being using in simulation {:?}",
                rec
            );
        } else {
            info!("CSV file not availabe/no state to read");
        }
    }
}
