//! The control crate is used for all control calculations for BIT.
//!
//! This crate implements a struct `Ctrl` complete with ans initialization and methods that are required for calculating the control torques for the BIT gondola (simulation).

pub mod error;
pub mod gains;
pub mod motors;
pub mod state;
pub mod flex_control;

use error as er;
use gains as gns;
use motors as mot;
use state as st;
use flex_control as fc;

// import crate for matrix math
extern crate nalgebra as na;

// use std::fmt;
pub use std::f64::consts::{E, PI};

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

// Defining a struct to hold the definitions of the system
// matrices
#[derive(Debug)]
pub struct Ctrl {
    /// Struct to represent the roll frameless motors
    pub fmot_roll: mot::TorqueMotor,
    /// Struct to represent the pitch frameless motors
    pub fmot_pitch: mot::TorqueMotor,
    /// Struct to represent the reaction wheel motor
    pub rw: mot::TorqueMotor,
    /// Struct to represent the pivot stepper motor
    pub pivot: mot::StepperMotor,
    /// Struct to represent the fine pointing controller gains
    pub fine_gains: gns::Gains,
    /// Struct to represent the coarse pointing controller gains
    pub slew_gains: gns::Gains,
    /// Struct to represent the error state of the system
    pub error: er::Error,
    /// Struct to represent the state of the system
    pub state: st::State,
    /// flag to signal if the gondola should be slewing
    pub slew_flag: bool,
}

impl Ctrl {
    /// Function that calculates the requested control torques based on the current error and integrated error
    ///
    /// # Detailed Explanation
    ///
    /// This function takes the current yaw pitch roll error, integrated error and rate error
    /// (gyro measurement) and calculates the torques to request on each axis using a PID controller
    /// $$\tau_{req} = K_p\phi_{err rate,k} + K_i \phi_{err rate,sum,k} $$,
    ///  
    ///  -these error terms are calculated in functions implemented on the `error` struct
    ///
    /// # Arguments
    ///  -`self.error.err_th:Vector3<f64>`: $\phi_{err,k}$ the current error in euler angles
    ///
    ///  -`self.error.err_sum:Vector3<f64>`: $\phi_{err,sum,k}$ current integrated error (with a decay term)
    ///
    ///  -`self.error.err_om_k:Vector3<f64>`: $\omega_k$ the current rotational velocity error
    ///
    ///  -`gains.kp:Matrix3<f64>`: $K_p$ the proportional gain constants
    ///
    ///  -`gains.kp:Matrix3<f64>`: $K_i$ the integral gain constants
    ///
    ///  -`gains.kp:Matrix3<f64>`: $K_d$ the derivitive gain constants
    ///
    /// # Results
    /// - `self.fmot_roll.tau_applied: f64`: the applied torque in roll
    /// - `self.fmot_pitch.tau_applied: f64`: the applied torque in pitch
    /// - `self.fmot_rw.tau_applied: f64`: the applied torque in yaw
    pub fn calculate_applied_torque(&mut self) {
        trace!("calculate_applied_torque start");

        let mut temp2 = na::Vector3::<f64>::zeros();
        let mut temp3 = na::Vector3::<f64>::zeros();
        let mut tau_requested = na::Vector3::<f64>::zeros();

        let gains = self.fine_gains.clone();

        gains.kp.ad_mul_to(&self.error.err_rate, &mut tau_requested);
        // println!("{:} x {:} = {:}", &gains.kp, &self.error.err_rate, &tau_requested);

        gains.ki.mul_to(&self.error.err_rate_sum, &mut temp2);

        gains.kip.mul_to(&self.error.err_fine_sum, &mut temp3);
        // println!("{:} x {:} = {:}", &gains.ki, &self.error.err_rate_sum, &temp2);
        tau_requested = tau_requested + temp2 + temp3;
        // println!("tau_requested {}", &tau_requested);

        self.fmot_roll.tau_request = 1.0*tau_requested[0];
        self.fmot_pitch.tau_request = 1.0*tau_requested[1];
        // reverse rw torque because opposite acts on gondola (for simulation)
        self.rw.tau_request = -1.0*tau_requested[2];

        self.fmot_roll.bound_requested_torque();
        self.fmot_pitch.bound_requested_torque();
        self.rw.bound_requested_torque();
        // println!(
        //     "tau applied roll {:} \n pitch: {:} \n rw: {:}",
        //     self.fmot_roll.tau_applied, self.fmot_pitch.tau_applied, self.rw.tau_applied
        // );
        trace!("calculate_applied_torque end");
    }
    /// Function that ultimately updates the applied torque
    ///
    /// # Detailed Explanation
    ///
    /// This function completes the following steps to calculate the applied torque and applied pivot speed:
    /// - `self.update_error_terms()`
    /// - `_phi = self.calculate_coupling_matrix();`
    /// - `self.calculate_applied_torque(_phi));`
    /// - `self.pivot.calculate_pivot_speed(self.rw);`
    ///
    pub fn update_ctrl(&mut self) {
        trace!("update_ctrl start");
        self.state.eq_d.calculate_rotation_matrix();
        // update Chat_k based on yaw pitch roll in xhat_k
        self.state.update_desired_gimbal_rpy();
        let phi = self.state.calculate_coupling_matrix();
        self.state.gmb_k.calculate_inverse_gimbal_mapping_matrix();
        self.slew_flag = self.error
            .update_pointing_positional_error(&self.state, phi, self.slew_flag);
        self.error
            .update_pointing_velocity_error_terms(&mut self.state, &self.slew_flag);
        println!("slew flag {}", &self.slew_flag);

        self.calculate_applied_torque();
        self.pivot.calculate_pivot_speed(&self.rw);
        trace!("update_ctrl end");
    }

    /// Function to initialize the control matrix with default values
    ///
    pub fn new() -> Ctrl {
        // DEFINE SYSTEMATIC TERMS
        // ---------------------------------
        // define params struct
        let ctrl: Ctrl = Ctrl {
            fmot_roll: mot::TorqueMotor::fmot_new(),
            fmot_pitch: mot::TorqueMotor::fmot_new(),
            rw: mot::TorqueMotor::rw_new(),
            pivot: mot::StepperMotor::pivot_new(),
            fine_gains: gns::Gains::new(),
            slew_gains: gns::Gains::new(),
            error: er::Error::new(),
            state: st::State::new(),
            slew_flag: true,
        };
        info!("Ctrl struct initialized");
        debug!("Control initialized: \n {:?}", ctrl);
        return ctrl;
    }
}
