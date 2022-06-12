//! The Gains submodule details the structs and implementations for the control gains.
//!
//! This crate implements a struct `Gains` complete with implementations necessary to
//! read in new gains from an input source and re-build the gain matrices for PID control

// import crate for matrix math
extern crate nalgebra as na;

#[derive(Debug, Clone)]
/// a struct to hold the gains for the PID control
pub struct Gains {
    /// proportional gain diagonal matrix
    pub kp: na::Matrix3<f64>,
    /// Positional Integral gain diagonal matrix
    pub kip: na::Matrix3<f64>,
    /// integral gain diagonal matrix
    pub ki: na::Matrix3<f64>,
}

impl Gains {
    /// Function that updates the gain matricies
    ///
    /// # Detailed Explanation
    ///
    /// This function is called when new gains are read in. The gains are read into a vector
    /// and stored in the ctrl struct. This function builds the matricies and places them back
    /// into the struct ctrl.
    ///
    /// # Arguments
    ///
    /// -`self.kp_vec: [f64;3]` : a vector of the proportional control gains for each axis
    /// -`self.ki_vec: [f64;3]` : a vecton of the integral control gains for each axis
    /// -`self.kd_vec: [f64;3]` : a vector of the derivative control gains for each axis
    ///
    /// # Results
    ///
    /// -`self.kp_f: Matrix3<f64>`: A diagonal matrix of the proportional control gains from `self.kp_vec`
    /// -`self.ki_f: Matrix3<f64>`: A diagonal matrix of the integral control gains from `self.ki_vec`
    /// -`self.kd_f: Matrix3<f64>`: A diagonal matrix of the derivative control gains from `self.kd_vec`
    pub fn update_gain_matrices(
        &mut self,
        kp_vec: &[f64; 3],
        kip_vec: &[f64; 3],
        ki_vec: &[f64; 3],
    ) -> () {
        self.kp = na::Matrix3::<f64>::from_partial_diagonal(kp_vec);
        self.ki = na::Matrix3::<f64>::from_partial_diagonal(ki_vec);
        self.kip = na::Matrix3::<f64>::from_partial_diagonal(kip_vec);
    }

    /// Function to instantiate a new Gains struct with default values
    pub fn new() -> Gains {
        let kp_vec = [0.10, 0.10, 3.0];
        let ki_vec = [0.00001, 0.000, 0.0001];
        let kip_vec = [0.010, 0.0010, 0.01];
        let kp = na::Matrix3::<f64>::from_partial_diagonal(&kp_vec);
        let kip = na::Matrix3::<f64>::from_partial_diagonal(&kip_vec);
        let ki = na::Matrix3::<f64>::from_partial_diagonal(&ki_vec);

        let gains: Gains = Gains { kp, kip, ki };

        gains
    }
}
