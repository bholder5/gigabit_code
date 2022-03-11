//! The control crate is used for all control calculations for BIT.
//!
//! This crate implements a struct `Ctrl` complete with ans initialization and methods that are required for calculating the control torques for the BIT gondola (simulation).

pub mod gains;
pub mod motors;
pub mod state;
pub mod error;

use gains as gns;
use motors as mot;
use state as st;
use error as er;

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
    /// Fine pointing motors max torque
    pub fmot_roll: mot::TorqueMotor,
    pub fmot_pitch: mot::TorqueMotor,
    pub rw: mot::TorqueMotor,
    pub pivot: mot::StepperMotor,
    pub fine_gains: gns::Gains,
    pub slew_gains: gns::Gains,
    pub error: er::Error,
    pub state: st::State,
    pub slew_flag: bool,
}

impl Ctrl {
    /// Function that calculates the requested control torques based on the current error and integrated error
    /// # Detailed Explanation
    /// This function takes the current yaw pitch roll error, integrated error and rate error (gyro measurement) and calculates the torques to request on each axis using a PID controller
    /// $$\tau_{req} = K_p\phi_{err,k} + K_d \omega_{err,k} + K_i \phi_{err,sum,k} $$, where $\omega_{err,k} = \omega_k$ for fine pointing
    ///  
    ///  -these error terms are calculated in `update_error_terms`
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
    fn calculate_applied_torque(
        &mut self,
    ) {
        trace!("calculate_applied_torque start");
        let mut temp = na::Vector3::<f64>::zeros();
        let mut temp2 = na::Vector3::<f64>::zeros();
        let mut tau_requested = na::Vector3::<f64>::zeros();

        let mut gains = self.fine_gains.clone();
        let phi = na::Matrix3::<f64>::identity();

        phi.mul_to(&self.error.err_rate, &mut temp);
        info!("{:} x {:} = {:}", &phi, &self.error.err_rate, &temp);
        gains.kp.ad_mul_to(&temp, &mut tau_requested);
        info!("{:} x {:} = {:}", &gains.kp, &temp, &tau_requested);

        // phi.mul_to(&(self.error.err_om), &mut temp);
        // info!("{:} x {:} = {:}", &phi, &self.error.err_om, &temp);
        // gains.kd.mul_to(&temp, &mut temp2);
        // info!("{:} x {:} = {:}", &gains.kd, &temp, &temp2);
        // tau_requested = tau_requested + temp2;

        phi.mul_to(&self.error.err_rate_sum, &mut temp);
        info!("{:} x {:} = {:}", &phi, &self.error.err_rate_sum, &temp);
        gains.ki.mul_to(&temp, &mut temp2);
        info!("{:} x {:} = {:}", &gains.ki, &temp, &temp2);
        tau_requested = tau_requested + temp2;

        self.fmot_roll.tau_request = tau_requested[0];
        self.fmot_pitch.tau_request = tau_requested[1];
        // reverse rw torque because opposite acts on gondola
        self.rw.tau_request = -tau_requested[2];

        self.fmot_roll.bound_requested_torque();
        self.fmot_pitch.bound_requested_torque();
        self.rw.bound_requested_torque();
        debug!(
            "tau applied roll {:} \n pitch: {:} \n rw: {:}",
            self.fmot_roll.tau_applied, self.fmot_pitch.tau_applied, self.rw.tau_applied
        );
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
        self.error.update_pointing_positional_error(&self.state, phi);
        self.error.update_pointing_velocity_error_terms(&self.state, &self.slew_flag);

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