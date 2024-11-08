// import crate for matrix math
extern crate nalgebra as na;

pub mod correction;
pub mod gyros;
pub mod propogation;

// use std::fmt;
use crate::control::state::equitorial as eq;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
//define types
pub type Vector9 = na::SMatrix<f64, 9, 1>;
pub type Vector18 = na::SMatrix<f64, 18, 1>;

// Defining a struct to hold the definitions of the system
// matrices
#[derive(Debug)]
pub struct Estimator {
    // state terms
    /// the boresight gyros and accompanying calibration terms
    pub gyros_bs: gyros::GyroBs,
    pub eq_hat_k: eq::Equatorial,
    /// Progation struct
    pub prop: propogation::Propogation,
    //Correction terms
    pub corr: correction::Correction,
    pub reset: bool,
    /// piv estimation terms for lqr
    pub piv_angle: f64,
    pub piv_speed: f64,
    
}
impl Estimator {
    /// Function to instantiate a new Est struct
    ///
    /// # Detailed Explanation
    ///
    /// This function instantiates a new Est struct with default values
    pub fn new() -> Estimator {
        // angular velocity measurement (gyros)
        let gyros_bs = gyros::GyroBs::new();
        let eq_hat_k = eq::Equatorial::new();
        let prop = propogation::Propogation::new();
        let corr = correction::Correction::new();
        let reset = true;
        let piv_angle = 0.0;
        let piv_speed: f64 = 0.0;
        // DEFINE SYSTEMATIC TERMS
        // ----------------------------------

        // define params struct
        let est = Estimator {
            gyros_bs,
            eq_hat_k,
            prop,
            corr,
            reset,
            piv_angle,
            piv_speed,
        };
        info!("Estimation initialized: \n {:?}", est);
        return est;
    }
    pub fn propogate(&mut self) {
        self.prop.propogate(&mut self.eq_hat_k, &self.gyros_bs);
        self.piv_angle = 0.999999999*self.piv_angle + 0.001*self.piv_speed;
        // print!("piv angle {:.3} piv speed {:.5} \n", self.piv_angle, self.piv_speed);
    }
    pub fn correct_estimate(&mut self) {
        self.corr.correct_estimate(&mut self.eq_hat_k, &mut self.prop, &mut self.gyros_bs);
    }
}
