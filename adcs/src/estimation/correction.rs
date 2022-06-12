//! Struct and finctuionality to correct the propogated
//!  estimated state using absolute measurements

use nalgebra as na;
use crate::estimation::propogation as prp;
use crate::estimation::gyros as gy;
use crate::control::state::equitorial as eq;
use std::time::{Duration, Instant};

type Matrix3x18 = na::SMatrix<f64, 3, 18>;
type Vector9 = na::SVector<f64, 9>;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

///Struct holding the necessary terms for the kalman filter correction step
#[derive(Debug)]
pub struct Correction{
    pub _r_k: na::Matrix3<f64>,
    pub _m_k: na::Matrix3<f64>,
    pub _h_k: Matrix3x18,
    pub _mrmi: na::Matrix3<f64>,
    pub eq_act: eq::Equatorial,  
}

impl Correction{
    /// Function to instantiate a new Correction struct
    /// 
    /// # Detailed Explanation
    /// 
    /// This function instantiates a new Correction struct with default values
    pub fn new() -> Correction{
        let mut corr = Correction{
            _r_k: na::Matrix3::identity() * 1e-8,
            _m_k: na::Matrix3::identity(),
            _h_k: Matrix3x18::identity(),
            _mrmi: na::Matrix3::identity(),
            eq_act: eq::Equatorial::new(),
        };
        corr._mrmi = corr._m_k * corr._r_k * corr._m_k.transpose();

        corr
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
    pub fn correct_estimate(&mut self, eq_hat: &mut eq::Equatorial, prop: &mut prp::Propogation, gyros_bs: &mut gy::GyroBs){

        // calculate the error between the estimation and the most recent absolut measurement
        let err_mat = self.eq_act.rot * eq_hat.rot.inverse();
        let err_eul = err_mat.inverse().euler_angles();
        let err = na::Vector3::<f64>::new(err_eul.0, err_eul.1, err_eul.2);
        // println!("Err mat {:.6} err euler angles {} {} {} err vector {:}", err_mat, err_eul.0, err_eul.1, err_eul.2, err);

        let wk: na::Matrix3::<f64> = (self._h_k * prop.p_hat * self._h_k.transpose()) + (self._mrmi);
        // println!("WK {:}", wk);
        let wk_d_i = match wk.cholesky() {
            Some(d) => d.inverse(),
            None => {println!("this is happening for some reason"); error!("Error: inverse failed"); return;},
        };
        //Faster but more rounding error/ more unstable?
        // let wk_d_i = match wk.try_inverse() {
        //     // this gives inverse of wk, not wk_d
        //     Some(x) => x,
        //     None => {error!("Error: inverse failed"); return;}
        //     };
        
        let kalman_gain = prop.p_hat * self._h_k.transpose() * wk_d_i;
        // println!("KG {:}", kalman_gain);
        let zeta = kalman_gain * err;
        // println!("zeta {:}", zeta);
        // Decompose Zeta into the correction factors for each various term
        let zeta_bek  = na::Vector3::<f64>::new(zeta[0], zeta[1], zeta[2]);
        let zeta_cbgk = na::Vector3::<f64>::new(zeta[3], zeta[4], zeta[5]);
        let zeta_abg = Vector9::from_row_slice(&[zeta[6], zeta[7], zeta[8], zeta[9], zeta[10], zeta[11], zeta[12], zeta[13], zeta[14]]);
        let zeta_bgk = na::Vector3::<f64>::new(zeta[15], zeta[16], zeta[17]);
        
        // Correct the orientation estimate
        let zeta_bek_rot = na::Rotation3::from_axis_angle(
            &na::Unit::new_normalize(zeta_bek.clone()),
            zeta_bek.norm().clone(),
        )
        .inverse();
        // println!("zeta_bek_rot {:.10}", zeta_bek_rot);
        eq_hat.rot = zeta_bek_rot * eq_hat.rot;
        eq_hat.extract_ra_dec_fr_from_rotmat();

        // Correct the gyro position estimate
        let zeta_cbgk_rot = na::Rotation3::from_axis_angle(
            &na::Unit::new_normalize(zeta_cbgk.clone()),
            zeta_cbgk.norm().clone(),
        )
        .inverse();
        // println!("zeta_cbgk_rot {:.15}", zeta_cbgk_rot);
        gyros_bs.c_bg = zeta_cbgk_rot * gyros_bs.c_bg;

        // Correct gyro internal alignment and scaling
        let zeta_abg_mat = na::Matrix3::from_row_slice(zeta_abg.as_slice());
        // println!("zeta_abg_mat {:+e}", zeta_abg_mat);
        gyros_bs.a_g = zeta_abg_mat + gyros_bs.a_g;

        // Correct the gyro bias estimate
        gyros_bs.bias = gyros_bs.bias + zeta_bgk;

        // Update the covariance matrix
        prop.p_hat = prop.p_hat - (kalman_gain * self._h_k * prop.p_hat) - (prop.p_hat * self._h_k.transpose() * kalman_gain.transpose()) + (kalman_gain * wk * kalman_gain.transpose());
        // println!("a_g {:}, c_g {:}, b_g {:}", gyros_bs.a_g, gyros_bs.c_bg, gyros_bs.bias);    
    }
        
    pub fn read_LIS(&mut self, lis_rot: &na::Rotation3::<f64>){
        self.eq_act.rot = lis_rot.clone();
        self.eq_act.extract_ra_dec_fr_from_rotmat();
    }
}