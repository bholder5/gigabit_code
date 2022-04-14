//! Struct and finctuionality to correct the propogated
//!  estimated state using absolute measurements

use nalgebra as na;
use crate::estimation::propogation as prp;
use crate::control::state::equitorial as eq;

type Matrix3x18 = na::SMatrix<f64, 3, 18>;

///Struct holding the necessary terms for the kalman filter correction step
#[derive(Debug)]
pub struct Correction{
    pub _r_k: na::Matrix3<f64>,
    pub _m_k: na::Matrix3<f64>,
    pub _h_k: Matrix3x18,
    pub eq_act: eq::Equatorial,  
}

impl Correction{
    /// Function to instantiate a new Correction struct
    /// 
    /// # Detailed Explanation
    /// 
    /// This function instantiates a new Correction struct with default values
    pub fn new() -> Correction{
        Correction{
            _r_k: na::Matrix3::identity() * 1e-8,
            _m_k: na::Matrix3::identity(),
            _h_k: Matrix3x18::identity(),
            eq_act: eq::Equatorial::new(),
        }
    }

    /// Function to correct the estimated state using absolute measurements
    /// 
    /// # Detailed Explanation
    /// 
    /// This function corrects the estimated state using absolute measurements
    /// the majority of which depends on star camera sensors. THe initial 
    /// implementation will assume a valid star camera solution is ready at 
    /// a set cadence before introducing the combination of star cameras and
    ///  the irregular cadence. The correction step follows the kalman filter correction.
    /// 
    /// # Arguments
    /// 
    /// `eq_hat: eq::Equatorial` - Estimated equatorial coordinates
    /// `prop: prp::Propogation` - Propogation struct
    /// 
    /// # Returns
    /// 
    /// `eq::Equatorial` - Corrected estimated equatorial coordinates
    /// `prop: prp::Propogation` - Propogation struct with correctiong covariance
    pub fn correct_estimate(eq_hat: &mut eq::Equatorial, prop: &mut prp::Propogation){

    }

    }

}