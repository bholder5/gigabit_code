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

use adcs::{control::*, estimation::*, miscellaneous::*, verification::*};
use initialization::*;
use json_handling as js;
use logging::*;
use sim_csv as sc;

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
    test_control();
    // run initialization and have all relevant
    // bit parameters in the struct
    let mut bp = init_bit();
    bp.get_orientation_vec();
    let est = Est::init();
    let mut ctrl = Ctrl::new();

    // define Simulation parameters
    let mut step: u64 = 0; //integer for tracking rate of control
    let mut tau_applied: [f64; 9] = [0., 0., 0., 0., 0., 0., 0., 0., 0.];
    let mut t = 0.0;

    //********************************************************************
    //********************************************************************
    //                   SIMULATION BEGINS
    //********************************************************************

    let mut y_result: [f64; 21] = [0.; 21];

    js::read_gains(&mut ctrl); // read in gains from json file (for tuning)
                               // ctrl.update_gain_matrices();
                               // ctrl.update_desired_orientation_matrix();
                               // read_last_state(&mut bp, &mut est, &mut ctrl);
                               // bp.get_orientation_vec();

    tau_applied[6] = ctrl.rw.tau_applied.clone(); //yaw
    tau_applied[7] = ctrl.fmot_roll.tau_applied.clone(); //roll
    tau_applied[8] = ctrl.fmot_pitch.tau_applied.clone(); //pitch

    sc::push_record(&t, &bp, &est, &ctrl);

    trace!("START");
    for _step in 0..0 as usize {
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
        // println!("Pivot request {:}", ctrl.pivot.omega_request);

        step = step + 1;
        t = t + 0.001;
        bp.x = Vector21::from_row_slice(y_result.as_ref());

        // reassign the state matrix to the last state matrix
        //send theta to bp for use in calculating the angular velocity vector
        bp.update_state();
        bp.get_omega_meas();
        bp.get_orientation_vec();
        bp.get_rw_speed();

        // debug!("new state array is: {:}", &bp.x);
        // Calculate control terms based on current step
        // println!("check for control calcs: 20 % step = {:}, step: {:}", (step % 20), step);

        if (step % 1) < 1 {
            // js::read_gains(&mut ctrl); // read in gains from json file (for tuning)

            ctrl.state
                .update_current_equatorial_coordinates(&bp.phi_act);
            // grab_vec3(&mut ctrl.state.omega, &bp.omega_m);
            ctrl.state.omega[0] = bp.omega_m[0].clone();
            ctrl.state.omega[1] = bp.omega_m[1].clone();
            ctrl.state.omega[2] = bp.omega_m[2].clone();
            ctrl.rw.omega = bp.omega_rw;
            ctrl.state
                .gmb_k
                .update_gimbal_coordinates(&[bp.x[16], bp.x[17], bp.x[15]]);

            ctrl.update_ctrl();

            // update actual torque vector with applied torques
            tau_applied[6] = ctrl.rw.tau_applied; //yaw
            tau_applied[7] = ctrl.fmot_roll.tau_applied; //roll
            tau_applied[8] = ctrl.fmot_pitch.tau_applied; //pitch

            // if (step % 1000) < 1{
            //     println!("step: {:}", step);
            // }
        }

        // record the data
        sc::push_record(&t, &bp, &est, &ctrl);
    }
    trace!("END");
}
