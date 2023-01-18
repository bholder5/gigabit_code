#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

// include!("bindings.rs");

#[allow(dead_code)]
// #[allow(unused_attributes)]
mod bindings;
mod initialization;
mod json_handling;
mod logging;
mod sim_csv;
mod measurements;
mod disturbances;
mod flex_sim;

use adcs::{control::*, estimation::*, verification::*};
use initialization::*;
use json_handling as js;
use logging::*;
use sim_csv as sc;
use measurements as ms;
use disturbances as dist;
use flex_sim as f_sim;

use bindings::bit_one_step;

pub use log::LevelFilter;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
pub use log4rs::append::file::FileAppender;
pub use log4rs::config::{Appender, Config, Root};
pub use log4rs::encode::pattern::PatternEncoder;

extern crate csv;
extern crate serde;
extern crate time;
// #[macro_use]
extern crate serde_derive;

fn main() {

    init_log();
    // the verification function for all tihngs in the control module.
    // test_control();
    // test_estimation();
    // run initialization and have all relevant
    // bit parameters in the struct
    let mut bp = init_bit();
    let mut flex = f_sim::init_flex_model();
    bp.get_orientation_vec(&flex);
    

    let mut est = Estimator::new();
    let mut ctrl = Ctrl::new();
    let mut sim_state = state::State::new();
    let mut meas = ms::Meas::new();
    let mut late_meas = ms::Meas::new();
    let mut wind = dist::Wind_torque::new();

    // define Simulation parameters
    let mut step: u64 = 0; //integer for tracking rate of control
    let mut tau_applied: [f64; 9] = [0., 0., 0., 0., 0., 0., 0., 0., 0.];
    let mut t = 0.0;

    //********************************************************************
    //********************************************************************
    //                   SIMULATION BEGINS
    //********************************************************************

    let mut y_result: [f64; 21] = [0.; 21];

    js::read_gains(&mut ctrl, &mut bp, &mut est); // read in gains from json file (for tuning)
                               // ctrl.update_gain_matrices();
                               // ctrl.update_desired_orientation_matrix();
    // sc::read_last_state(t, &mut bp, &mut est, &mut ctrl, &mut meas, &mut sim_state);
                               // bp.get_orientation_vec();

    tau_applied[6] = ctrl.rw.tau_applied.clone(); //yaw
    tau_applied[7] = ctrl.fmot_roll.tau_applied.clone(); //roll
    tau_applied[8] = ctrl.fmot_pitch.tau_applied.clone(); //pitch

    sc::push_record(&t, &bp, &est, &ctrl, &meas, &sim_state, &flex).unwrap();

    let mut first_LIS = true;

    // Initialize with an achievable target
    ctrl.state.gmb_d.roll = 0.1;
    ctrl.state.gmb_d.pitch = -0.6;
    ctrl.state.gmb_d.yaw = 1.0;
    ctrl.state.gmb_d.calculate_rotation_matrix();
    ctrl.state.update_desired_eq_from_gmb();


    trace!("START");
    for _step in 0..10 as usize {
        
        ///////// beginning of the simulation loop
        /////////////////////////////////////////
        unsafe {
            trace!("bit_one_step start");
            bit_one_step(
                bp.x.as_ptr(),
                tau_applied.as_ptr(),
                bp.unlock.as_ptr(),
                ctrl.pivot.omega_request,
                false as u8,
                bp._dt,
                bp._num_steps,
                bp._tau_piv_max,
                bp.pitch_nom,
                y_result.as_mut_ptr(),
            );
            trace!("bit_one_step end");
        }
        
        if flex.flex_enable{
            flex.propogate_flex(&[tau_applied[6], tau_applied[7], tau_applied[8]], bp._dt, bp._num_steps);
        // flex.propogate_flex(&[1., 1., 1.], bp._dt, bp._num_steps);
        }

        // println!("Pivot request {:}", ctrl.pivot.omega_request);

        step = step + 1;
        t = t + 0.001;
        
        late_meas = meas.clone();
        meas.gps.utc = meas.gps.utc + chrono::Duration::milliseconds(1);
        sim_state.gps = meas.gps.clone();

        // reassign the state matrix to the last state matrix
        //send theta to bp for use in calculating the angular velocity vector
        bp.x = Vector21::from_row_slice(y_result.as_ref());
        bp.update_state();
        bp.get_omega_meas();
        bp.get_orientation_vec(&flex);
        bp.get_rw_speed();

        sim_state.update_current_equatorial_coordinates(&bp.phi_act);

        //Meas is the struct holding the actual gyro calibration terms
        meas.cbh = sim_state.hor.rot.clone();
        meas.read_measurements(&bp, &sim_state, &flex);

        if bp.latency {
            est.gyros_bs.read_gyros(late_meas.gyros_bs.omega_m, t.clone());
            ctrl.state.gps = late_meas.gps.clone();
        } else {
            est.gyros_bs.read_gyros(meas.gyros_bs.omega_m, t.clone());
            ctrl.state.gps = meas.gps.clone();
        }

        
        est.propogate();

        if (step % 1000) < 1 {
            est.corr.read_LIS(&sim_state.eq_k.rot);
            if {est.reset | first_LIS} {
                est.eq_hat_k.rot = sim_state.eq_k.rot.clone();
                est.gyros_bs.reset();
                est.reset = false;
                first_LIS = false;
            }
            est.correct_estimate();  
        }
        
        if (step % 100) < 1 {  
            if bp.ctrl_from_est {
                ctrl.state.eq_k = est.eq_hat_k.clone();
                ctrl.state.omega[0] = est.gyros_bs.omega_k[0].clone();
                ctrl.state.omega[1] = est.gyros_bs.omega_k[1].clone();
                ctrl.state.omega[2] = est.gyros_bs.omega_k[2].clone();
                if bp.latency {
                    ctrl.state
                        .gmb_k
                        .update_gimbal_coordinates(&[late_meas.roll, late_meas.pitch, late_meas.yaw_p]);
                } else {
                    ctrl.state
                        .gmb_k
                        .update_gimbal_coordinates(&[meas.roll, meas.pitch, meas.yaw_p]);
                }
            } else {
                ctrl.state.eq_k = sim_state.eq_k.clone();
                ctrl.state.omega[0] = bp.omega_m[0].clone();
                ctrl.state.omega[1] = bp.omega_m[1].clone();
                ctrl.state.omega[2] = bp.omega_m[2].clone();
                ctrl.state
                .gmb_k
                .update_gimbal_coordinates(&[bp.x[16], bp.x[17], -sim_state.hor.az.clone()]);
            }

            // grab_vec3(&mut ctrl.state.omega, &bp.omega_m);
            
            ctrl.rw.omega = bp.omega_rw;
            

            ctrl.update_ctrl();

            //Disturbance generation
            wind.generate();
            tau_applied[0] = wind.tau[0].clone();
            tau_applied[1] = wind.tau[1].clone();
            tau_applied[2] = wind.tau[2].clone();

            // update actual torque vector with applied torques
            tau_applied[6] = 1.0*ctrl.rw.tau_applied; //yaw
            tau_applied[7] = 1.0*ctrl.fmot_roll.tau_applied; //roll
            tau_applied[8] = 1.0*ctrl.fmot_pitch.tau_applied; //pitch

            // if (step % 1000) < 1{
            //     println!("step: {:}", step);
            // }
        }

        // record the data
        sc::push_record(&t, &bp, &est, &ctrl, &meas, &sim_state, &flex).unwrap();
        js::read_gains(&mut ctrl, &mut bp, &mut est); // read in gains from json file (for tuning)

    }
    trace!("END");
}
