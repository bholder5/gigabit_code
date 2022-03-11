extern crate nalgebra as na;
use std::f64::consts::{E, PI};
extern crate statrs as stat;

use crate::control::state as st;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};


use crate::miscellaneous as misc;

#[derive(Debug)]
pub struct Error {
    pub err_gmb_th: na::Vector3<f64>,
    pub err_fine_th: na::Vector3<f64>,
    pub err_th: na::Vector3<f64>,
    pub err_rate: na::Vector3::<f64>,
    pub err_rate_sum: na::Vector3<f64>,
    pub _err_tc: f64,
    pub _err_decay: f64,
    pub _ctrl_dt: f64,
    pub rot_err: na::Rotation3<f64>,
    pub u_lower: f64,
    pub u_upper: f64,
}

impl Error {
    pub fn new() -> Error{
        let _err_tc: f64 = 10.0;
        let _ctrl_dt: f64 = 0.01;
        let _err_decay: f64 = E.powf(-_ctrl_dt / _err_tc);
        let err_gmb_th = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_fine_th = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_th = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_rate = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_rate_sum = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let rot_err = na::Rotation3::<f64>::identity();
        let u_lower = 0.00001;
        let u_upper = 0.001;

        let error: Error = Error {
            
            err_fine_th,
            err_gmb_th,
            err_th,
            err_rate,
            err_rate_sum,
            rot_err,
            _err_tc,
            _ctrl_dt,
            _err_decay,
            u_lower,
            u_upper,
        };

        error
        
    }

    /// Function that updates the error terms for control
    ///
    /// # Detailed Explanation
    ///
    /// This function takes $\hat x_k$ that is passed to the ctrl struct from the estimation algorithm
    /// and follows the following steps to calculate the necessary error terms for control:
    /// - $\hat C_k = C_x(\hat x_1)C_y(\hat x_2)C_z(\hat x_3)$: Current orientation of the telescope
    /// - $C_{err} = \hat C_k \times C_{des}^T$: current error (rotation matrix)
    /// - $\phi_{err} = unxmat(1-C_{err})$
    /// - $\phi_{err,sum,k} = \phi_{err,sum,k-1} e^{dt_{ctrl}/\tau_{err}} + \phi_{err}dt$ which provides a decay term for
    /// the error sum to prevent integral windup
    ///
    /// # Arguments
    /// `self.xhat_k:Vector3<f64>` the current orientation estimate as euler angles
    ///
    /// `self.Cdes:Rotation3<f64>`  the desired orientation in rotation matrix format
    ///
    /// `self.err_decay:f64` Decay term (rate) for error sum
    ///
    /// # Results
    /// `self.Chat_k:Rotation3<f64>` the current orientation in rotation matrix format
    ///
    /// `self.err_sum_k:Vector3<f64>` the current error sum
    pub fn update_pointing_positional_error(&mut self, state: &st::State, phi: na::Matrix3::<f64>) {
        trace!("update_positional_pointing_error start");
        // calculate the fine pointing coupled error
        println!("eq desired: {} eq k: {}", state.eq_d.rot, state.eq_k.rot);
        self.rot_err = state.eq_d.rot.rotation_to(&state.eq_k.rot);
        let _rot_err = state.eq_k.rot.rotation_to(&state.eq_d.rot);
        
        misc::unxmat(
            &(na::Rotation3::<f64>::identity().matrix() - self.rot_err.matrix()),
            &mut self.err_fine_th,
        );

        println!("linear error rot: {}, error: {}, rot to eul {} {} {}", (na::Rotation3::<f64>::identity().matrix() - self.rot_err.matrix()), self.err_fine_th, _rot_err.euler_angles().0, _rot_err.euler_angles().1, _rot_err.euler_angles().2);

        let norm_fine_err = self.err_fine_th.norm();
        // println!("norm fine error: {}", norm_fine_err);
        
        if norm_fine_err < self.u_lower{
            self.err_th = phi * self.err_fine_th.clone();
        }
        else{
            let mut _err_gmb = na::Vector3::<f64>::new(
                state.gmb_d.roll - state.gmb_k.roll, 
                state.gmb_d.pitch - state.gmb_k.pitch,
                (state.gmb_d.yaw - state.gmb_k.yaw)
            );
            println!("gmb err: {}", _err_gmb);

            for ind in 0..3 {
                if _err_gmb[ind].abs() > PI{
                    let over = _err_gmb[ind] > PI;
                    let under = _err_gmb[ind] < PI;
                    _err_gmb[ind] = (_err_gmb[ind] - (2.0*PI* (over as u8 as f64)) + (2.0*PI* (under as u8 as f64)));

                }
            }
            println!("wrapped gmb err: {}", _err_gmb);

            self.err_gmb_th = _err_gmb.clone();
        
            if norm_fine_err > self.u_upper {
                // calculate the gimbal error
                self.err_th = self.err_gmb_th.clone();
            }
            else {
                let err_weight: f64 = (norm_fine_err - self.u_lower) / (self.u_upper-self.u_lower);

                let comb_err = ((1.0 - err_weight) * (phi*self.err_fine_th)) 
                                                                                + (err_weight * self.err_gmb_th);

                self.err_th = comb_err;
                println!("error weight: {}, error term 1: {}, error term 2: {}, gmb err: {}, fine err: {}", err_weight, ((1.0 - err_weight) * (phi*self.err_fine_th)), (err_weight * self.err_gmb_th), self.err_fine_th, self.err_gmb_th);
            }
            trace!("update_positional_pointing_error end");
        }
        
    }

    pub fn update_pointing_velocity_error_terms(&mut self, state: &st::State, slew_flag: &bool) {
        // println!("start of pointing velocity error terms");
        trace!("update_pointing_velocity_error_terms start");
        let _d_theta: na::Vector3::<f64> = state.gmb_k.gmm.cholesky().unwrap().inverse() * state.omega;

        let mut _roll_rate_des: f64 = 0.0;
        let mut _pitch_rate_des: f64 = 0.0;
        let mut _yaw_rate_des: f64 = 0.0;

        if slew_flag == &false {

            let g_const: f64 = 1.0;
            let g_lin: f64 = 100.0;

            let temp1: f64 = g_lin.powi(2) * PI / (4.0 * g_const.powi(2));
            let temp2: f64 = temp1.sqrt();
            let temp3: f64 =  (2.0 * g_const) / (g_lin * PI);

            let _roll_rate_des: f64 = -self.err_th.x.signum() * (2.0 * g_const * (self.err_th.x.abs() * 
                stat::function::erf::erf(self.err_th.x.abs() * temp2)) + 
                (temp3 * ((-self.err_th.x.powi(2) * temp1).exp()-1.0))).sqrt();

            let _pitch_rate_des: f64 = -self.err_th.y.signum() * (2.0 * g_const * (self.err_th.y.abs() * 
                stat::function::erf::erf(self.err_th.y.abs() * temp2)) + 
                (temp3 * ((-self.err_th.y.powi(2) * temp1).exp()-1.0))).sqrt();

            let _yaw_rate_des: f64 = -self.err_th.z.signum() * (2.0 * g_const * (self.err_th.z.abs() * 
                stat::function::erf::erf(self.err_th.z.abs() * temp2)) + 
                (temp3 * ((-self.err_th.z.powi(2) * temp1).exp()-1.0))).sqrt();
        }


        let err_gmb_rate: [f64; 3] = [(_d_theta.x - _roll_rate_des),
                                      (_d_theta.y - _pitch_rate_des),
                                      (_d_theta.z - _yaw_rate_des)];


        self.err_rate = na::Vector3::<f64>::from_row_slice(&err_gmb_rate);
        trace!("update_pointing_velocity_error_terms end");
        self.err_rate_sum = (self.err_rate_sum * self._err_decay)
            + (self.err_rate * self._ctrl_dt)
    }

}