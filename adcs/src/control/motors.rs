//! The Torque Motor submodule details the structs and implementations for the torque and stepper motors calculations and commanding.
//!
//! This crate defines structs `TorqueMotor` and `StepperMotor` complete with implementations necessary to 
//! safely update and command the motor
use std::f64::consts::{PI};
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
pub struct TorqueMotor {
    pub _tau_max: f64,
    pub _tau_thresh: f64,
    pub tau_applied: f64,
    pub tau_request: f64,
    pub _ts: f64,
    pub omega: f64,
}

impl TorqueMotor {
    /// Function takes in the requested motor torque and bounds it to within the motor limits
    ///
    ///  # Detailed Explanation
    ///
    /// This function is invoked once the requested torque of the motor is updated and sets the applied
    ///  torque to within the upper/lower limits of the motors torque capabilities.
    ///
    /// # Arguments
    ///
    /// - `self.tau_request` - requested torque output from the motor controller
    /// - `self.tau_max` - maximum torque the motor can output
    ///
    /// # Results
    ///
    /// - `self.tau_applied` - the applied torque sent to the motor
    pub fn bound_requested_torque(&mut self) {
        trace!("bound_requested_torque start");
        let mag = self.tau_request.abs();
        if mag > self._tau_max {
            self.tau_applied = self.tau_request.signum() * self._tau_max;
        } else {
            self.tau_applied = self.tau_request;
        }
        trace!("bound_requested_torque end");
    }

    pub fn fmot_new() -> TorqueMotor {
        let fmot: TorqueMotor = TorqueMotor {
            omega: 0.0,
            // define elevation/pitch fine motor limitations
            _tau_max: 11.14,              // max fine torque [Nm]
            _ts: 21.965 / 7.718 / 1000.0, // fine motor time constant [s]
            _tau_thresh: 0.001,           // threshhold torque to activate motors
            tau_applied: 0.0,
            tau_request: 0.0,
        };
        fmot
    }
    
    pub fn rw_new() -> TorqueMotor {
        let rw: TorqueMotor = TorqueMotor {
            omega: 0.0,
            _tau_max: 80.0,               // max fine torque [Nm]
            _ts: 21.965 / 7.718 / 1000.0, // fine motor time constant [s]
            _tau_thresh: 0.001,           // threshhold torque to activate motors
            tau_applied: 0.0,
            tau_request: 0.0,
        };
        rw
    }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////                      STEPPER MOTOR
//////////////////////////////////////////////////////////////////////////
#[derive(Debug)]
pub struct StepperMotor {
    pub omega: f64,
    pub omega_max: f64,
    pub omega_request: f64,
    pub omega_rw_nom: f64,
    pub gain1: f64,
    pub gain2: f64,
    pub _ts: f64,
}

impl StepperMotor {
    pub fn calculate_pivot_speed(&mut self, rw: &TorqueMotor) {
        trace!("calculate_pivot_speed start");
        let temp1: f64 = -self.gain1 * (rw.omega - self.omega_rw_nom);
        // println!("copntribution from temp1 {:}", temp1);
        let temp2: f64 = -self.gain2 * rw.tau_request;
        // println!("copntribution from temp2 {:}", temp2);
        self.omega_request = (temp1 + temp2);
        // println!("requested pivot speed is: {:}", self.omega_request);
        trace!("calculate_pivot_speed end");
    }

    pub fn pivot_new() -> StepperMotor{
        let gain1: f64 = 0.030;
        let k_flight_train: f64 = 0.001745329251994;
        let i_rw: f64 = 4.5;
        let gain2: f64 = 0.8 * (2.0 * ((gain1 / (i_rw * k_flight_train)).sqrt()));

        let pivot: StepperMotor = StepperMotor {
            omega: 0.0,
            omega_max: 0.0,
            omega_request: 0.0,
            omega_rw_nom: PI,
            gain1,
            gain2,
            _ts: 0.001,
        };
        pivot
    }
}