//! The Motor submodule details the structs and implementations for the torque and stepper motors calculations and commanding.
//!
//! This crate defines structs `TorqueMotor` and `StepperMotor` complete with implementations necessary to
//! safely update and command the motor
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use std::f64::consts::PI;

#[derive(Debug, Clone)]
/// a stuct to describe necessary parameters of a torque motor
pub struct TorqueMotor {
    /// maximum commandable torque
    pub _tau_max: f64,
    /// maximum sustainable torque
    pub _tau_thresh: f64,
    /// torque applied by the motor
    pub tau_applied: f64,
    /// torque requested from the controller
    pub tau_request: f64,
    /// time constant of the motor
    pub _ts: f64,
    /// current rotational velocity of the motor
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

    /// Function to instantiate a new TorqueMotor struct for a frameless motor with default values
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

    /// Function to instantiate a new TorqueMotor struct for a reaction wheel motor with default values
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
/// a struct to describe necessary parameters of a stepper motor
pub struct StepperMotor {
    /// current velocity of the stepper motor
    pub omega: f64,
    /// maximum velocity of the stepper motor
    pub omega_max: f64,
    /// requested velocity of the stepper motor from the controller
    pub omega_request: f64,
    /// nominal reaction wheel velocity
    pub omega_rw_nom: f64,
    /// first gain for pivot controller
    pub gain1: f64,
    /// second gain for pivot controller
    pub gain2: f64,
    /// stepper motor time constant
    pub _ts: f64,
}

impl StepperMotor {
    /// Function to calculate the desired pivot speed from the reaction commanded torque and current speed
    ///
    /// # Detailed Explanation
    ///
    /// This function is invoked once the requested torque of the reaction wheel motor is updated and determines the command for pivot speed according to the momentum dumping scheme outlined in J.Romualdez's thesis which is as follows
    /// $$\omega_{pivot} = -g_1 \left( \omega_{rw} - \omega_{rw,nom}\right) - g_2 \tau_{rw}$$
    ///
    /// # Arguments
    ///
    /// - `rw : &TorqueMotor` - the torque motor struct for the reaction wheel
    ///
    /// # Results
    ///
    /// - `omega_request` - the desired pivot speed
    pub fn calculate_pivot_speed(&mut self, rw: &TorqueMotor) {
        trace!("calculate_pivot_speed start");
        let temp1: f64 = -self.gain1 * (self.omega_rw_nom - rw.omega);
        // println!("copntribution from temp1 {:}", temp1);
        let temp2: f64 = -self.gain2 * rw.tau_applied;
        // println!("copntribution from temp2 {:}", temp2);
        let req = -(temp1 + temp2);

        if req.abs() > 0.2 {
            self.omega_request = req.signum() * 0.2;
        } else {
            self.omega_request = req;
        }

        // println!("requested pivot speed is: {:}", self.omega_request);
        trace!("calculate_pivot_speed end");
    }

    /// Function to instantiate a new StepperMotor struct for the pivot motor with default values
    pub fn pivot_new() -> StepperMotor {
        let gain1: f64 = 0.000005;
        let k_flight_train: f64 = 0.001745329251994;
        let i_rw: f64 = 4.5;
        let gain2: f64 = 0.05 * (2.0 * ((gain1 / (i_rw * k_flight_train)).sqrt()));

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
