//! This module performs two main functions, first to calculate 
//! the positional error in both the gimbal frame and in the
//! equatorial frame. These errors are then combined in a smooth
//!  function (to maintain stability) to use each inparticular 
//! scenarios. The second function takes the positional error 
//! and generates a desired velocity profile for each gimbal axis.
//!  This is all achieved using an Error struct with implementations.
extern crate nalgebra as na;
use std::f64::consts::{E, PI};
extern crate statrs as stat;

use crate::control::state as st;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

use crate::miscellaneous as misc;

#[derive(Debug)]
/// The error struct contains all the various components related to 
/// error necessary for gimbal control. There is a combination of 
/// calculated error terms and tuneable parameters.
pub struct Error {
    /// Positional gimbal error (roll, pitch, yaw)
    pub err_gmb_th: st::gimbal::Gimbal,
    /// Positional equatorial error (fr, dec, ra)
    pub err_b_th: st::gimbal::Gimbal,
    /// Combined positional error as a function of gimbal + equatorial error
    pub err_comb_th: na::Vector3<f64>,
    /// The error in gimbal rates (roll, pitch, yaw)
    pub err_rate: na::Vector3<f64>,
    /// Integrated gimbal rate error (roll, pitch, yaw)
    pub err_rate_sum: na::Vector3<f64>,
    /// Integrated rate error time constant of decay
    pub _err_tc: f64,
    /// Decay constant, function of control rate and decay time constant
    pub _err_decay: f64,
    /// Control algorithm time step
    pub _ctrl_dt: f64,
    /// Positional error in rotation matrix form
    pub rot_err: na::Rotation3<f64>,
    /// Lower bound constant for combined positional error calculation
    pub u_lower: f64,
    /// Upper bound constant for combined positional error calculation
    pub u_upper: f64,
}

impl Error {
    /// This function generates a new instance of the Error struct with default values
    pub fn new() -> Error {
        let _err_tc: f64 = 10.0;
        let _ctrl_dt: f64 = 0.01;
        let _err_decay: f64 = E.powf(-_ctrl_dt / _err_tc);
        let err_gmb_th = st::gimbal::Gimbal::new();
        let err_b_th = st::gimbal::Gimbal::new();
        let err_comb_th = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_rate = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_rate_sum = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let rot_err = na::Rotation3::<f64>::identity();
        let u_lower = 0.0001;
        let u_upper = 0.1;

        let error: Error = Error {
            err_b_th,
            err_gmb_th,
            err_comb_th,
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

    /// Function that updates the positional error terms for control calculations
    ///
    /// # Detailed Explanation
    ///
    /// This function takes $\hat x_k$ that is passed to the ctrl struct from the estimation algorithm
    /// and follows the following steps to calculate the necessary error terms for control:
    /// - $\hat C_k = C_x(\hat x_1)C_y(\hat x_2)C_z(\hat x_3)$: Current orientation of the 
    /// telescope in Equatorial coordiantes
    /// - $C_{err} =  C_{des} \times \hat C_{k}^T$: current error (rotation matrix)
    /// - $\phi_{err} = Parameterization(C_{err})$: The error in the telescope frame in euler angle form
    /// - $\theta_{err} = \theta_{des} - \theta_{k}$: The error in the gimbal frame
    /// - $\psi_{err,comb} = \alpha \theta_{err} + (1-\alpha) \phi_{err}$: The combined error 
    /// in the gimbal frame where $\alpha$ is a weighting factor determined by tunable upper and lower 
    /// bounds on the norm of the telescope frame
    ///
    /// # Arguments
    /// - `state: state::State` The current state object with desired and current orientation updated, 
    /// and updated gimbal mapping matrix (gmb_k.gmm)
    /// - `phi: na::Matrix3<f64>`  the coupling matrix
    ///
    /// # Results
    /// - `self.err_b_th: gimbal::Gimbal` the error in gimbal coordinates mapped from telescope frame error
    /// - `self.err_gmb_th: gimbal::Gimbal` the error in gimbal coordinates directly calculated from
    ///  current and desired gimbal coordinates
    /// - `self.err_comb_th: na::Vector3<f64>` the combined error for use in velocity profile calculation
    pub fn update_pointing_positional_error(&mut self, state: &st::State, phi: na::Matrix3<f64>) {
        trace!("update_positional_pointing_error start");
        
        //rot err in equatorial space
        let _rot_err = state.eq_d.rot * state.eq_k.rot.inverse();

        // parameterize the rotation error into euler angles
        let _rot_err_eul = _rot_err.inverse().euler_angles();

        let err_vec =
            na::Vector3::<f64>::new(_rot_err_eul.0, _rot_err_eul.1, _rot_err_eul.2);

        let norm_fine_err = err_vec.norm();
        println!("norm fine error: {}", norm_fine_err);

        let mapped_err = state.gmb_k.gmm.try_inverse().unwrap() * err_vec;


        self.err_b_th.roll = mapped_err[0];
        self.err_b_th.pitch = mapped_err[1];
        self.err_b_th.yaw = mapped_err[2];

        if norm_fine_err < self.u_lower {
            self.err_comb_th = phi * err_vec;
            println!("less than u.lower");
        } else {
            let mut _err_gmb = na::Vector3::<f64>::new(
                state.gmb_d.roll - state.gmb_k.roll,
                state.gmb_d.pitch - state.gmb_k.pitch,
                state.gmb_d.yaw - state.gmb_k.yaw,
            );
            // println!("gmbd{:?} \n\n gmbk{:?}",state.gmb_d, state.gmb_k);
            // println!("gmb err: {}, fine_err: {}", _err_gmb, phi* b_vec);

            for ind in 0..3 {
                if _err_gmb[ind].abs() > PI {
                    let over = _err_gmb[ind] > PI;
                    let under = _err_gmb[ind] < PI;
                    _err_gmb[ind] = _err_gmb[ind] - (2.0 * PI * (over as u8 as f64))
                        + (2.0 * PI * (under as u8 as f64));
                }
            }
            // println!("wrapped gmb err: {}", _err_gmb);

            self.err_gmb_th
                .update_gimbal_coordinates(&[_err_gmb[0], _err_gmb[1], _err_gmb[2]]);

            if norm_fine_err > self.u_upper {
                // calculate the gimbal error
                self.err_comb_th = _err_gmb.clone();
                println!("greater than u.upper");
            } else {
                println!("between u.lower and u.upper");
                let err_weight: f64 =
                    (norm_fine_err - self.u_lower) / (self.u_upper - self.u_lower);

                let comb_err = ((1.0 - err_weight) * (phi * err_vec)) + (err_weight * _err_gmb);

                self.err_comb_th = comb_err;
                println!("error weight: {}, combined error: {}\n term1 {} term2 {}", err_weight, comb_err, (1.0 - err_weight) * (phi * err_vec), (err_weight * _err_gmb));
            }
            trace!("update_positional_pointing_error end");
        }
    }
    /// $\phi_{err,sum,k} = \phi_{err,sum,k-1} e^{dt_{ctrl}/\tau_{err}} + \phi_{err}dt$ which 
    /// provides a decay term for the error sum to prevent integral windup

    /// Function to take the positional pointing error and determine a desired velocity profile for each gimbal axis
    /// 
    /// # Detailed Explanation 
    /// 
    /// This function takes the positional error and determines the desired velocity profile for each gimbal axis. It defaults to desired gimbal rates of 0, unless the slew flag is true in which case a desired velocity is calculated using the following formula per axis: $$\dot\theta_i = -sign(\theta_{i,err}) \sqrt{ 2G_{const} \left( ||\theta_{i,err}|| \text{erf} \left( ||\theta_{i,err}|| \frac{g_{lin} \sqrt{\pi}}{2G_{const}}\right) + \frac{2G_{const}}{\pi g_{lin}}\left[ \text{exp}\left( -\theta_{i,err}^2\frac{g_{lin}^2 \pi}{4G_{const}^2}\right)-1\right]\right)}$$
    /// 
    /// # Arguments	
    /// 
    /// `state: state::State` - The current state object, imported for the current gimbal angles and gimbal mapping matrix
    /// `slew_flag: bool` - Flag to determine if the gimbal is slewing
    /// `g_const: f64` - The gimbal constant (tuning parameter for velocity profile)
    /// `g_lin: f64` - The gimbal linearity constant (tuning parameter for velocity profile)
    /// 
    /// # Results 
    /// 
    /// - `self.err_rate: na::Vector3<f64>` - The error in gimbal rates
    /// - `self.err_rate_sum: na::Vector3<f64>` - The error in gimbal rates summed over time
    pub fn update_pointing_velocity_error_terms(
        &mut self,
        state: &mut st::State,
        slew_flag: &bool,
    ) {
        // println!("start of pointing velocity error terms");
        trace!("update_pointing_velocity_error_terms start");
        state.gmb_k.calculate_gimbal_mapping_matrix();
        let _d_theta: na::Vector3<f64> =
            state.gmb_k.gmm.cholesky().unwrap().inverse() * state.omega;

        let mut _roll_rate_des: f64 = 0.0;
        let mut _pitch_rate_des: f64 = 0.0;
        let mut _yaw_rate_des: f64 = 0.0;

        if slew_flag == &true {
            let g_const: f64 = 1.0;
            let g_lin: f64 = 100.0;

            let temp1: f64 = g_lin.powi(2) * PI / (4.0 * g_const.powi(2));
            let temp2: f64 = temp1.sqrt();
            let temp3: f64 = (2.0 * g_const) / (g_lin * PI);

            _roll_rate_des = -self.err_comb_th.x.signum()
                * (2.0
                    * g_const
                    * (self.err_comb_th.x.abs()
                        * stat::function::erf::erf(self.err_comb_th.x.abs() * temp2))
                    + (temp3 * ((-self.err_comb_th.x.powi(2) * temp1).exp() - 1.0)))
                    .sqrt();

            _pitch_rate_des = -self.err_comb_th.y.signum()
                * (2.0
                    * g_const
                    * (self.err_comb_th.y.abs()
                        * stat::function::erf::erf(self.err_comb_th.y.abs() * temp2))
                    + (temp3 * ((-self.err_comb_th.y.powi(2) * temp1).exp() - 1.0)))
                    .sqrt();
            _yaw_rate_des = -self.err_comb_th.z.signum()
                * (2.0
                    * g_const
                    * (self.err_comb_th.z.abs()
                        * stat::function::erf::erf(self.err_comb_th.z.abs() * temp2))
                    + (temp3 * ((-self.err_comb_th.z.powi(2) * temp1).exp() - 1.0)))
                    .sqrt();
        }

        let err_gmb_rate: [f64; 3] = [
            (_roll_rate_des - _d_theta.x),
            (_pitch_rate_des - _d_theta.y),
            (_yaw_rate_des - _d_theta.z),
        ];

        self.err_rate = na::Vector3::<f64>::from_row_slice(&err_gmb_rate);
        trace!("update_pointing_velocity_error_terms end");
        // println!("err_rate_sum {}, err_decay {}, err_rate {}, _ctrl_dt {}", self.err_rate_sum, self._err_decay, self.err_rate, self._ctrl_dt);
        println!("err_rate_sum {}", self.err_rate_sum);

        self.err_rate_sum = (self.err_rate_sum * self._err_decay) + (self.err_rate * self._ctrl_dt)
    }
    // archived function call
    pub fn update_pointing_positional_error2(&mut self, state: &st::State, phi: na::Matrix3<f64>) {
        trace!("update_positional_pointing_error start");
        // calculate the fine pointing coupled error
        // println!("eq desired: {} eq k: {}", state.eq_d.rot, state.eq_k.rot);
        // self.rot_err = state.eq_d.rot.rotation_to(&state.eq_k.rot);
        //rot err in body space
        let _rot_err = state.eq_k.rot * state.eq_d.rot.inverse();
        println!("current {:?} \n\n desired{:?}\n\n", state.eq_k, state.eq_d);
        // misc::unxmat(
        //     &(na::Rotation3::<f64>::identity().matrix() - self.rot_err.matrix()),
        //     &mut self.err_eq_th,
        // );

        // our rotation matrices are the opposite of most crates (including this one);
        self.err_b_th.rot = _rot_err;
        self.err_b_th.extract_gimbal_rpy();

        let b_vec =
            na::Vector3::<f64>::new(self.err_b_th.roll, self.err_b_th.pitch, self.err_b_th.yaw);
        let norm_fine_err = b_vec.norm();
        // println!("norm fine error: {}", norm_fine_err);

        if norm_fine_err < self.u_lower {
            self.err_comb_th = phi * b_vec;
        } else {
            let mut _err_gmb = na::Vector3::<f64>::new(
                state.gmb_d.roll - state.gmb_k.roll,
                state.gmb_d.pitch - state.gmb_k.pitch,
                state.gmb_d.yaw - state.gmb_k.yaw,
            );
            println!("gmbd{:?} \n\n gmbk{:?}", state.gmb_d, state.gmb_k);
            println!("gmb err: {}, fine_err: {}", _err_gmb, phi * b_vec);

            for ind in 0..3 {
                if _err_gmb[ind].abs() > PI {
                    let over = _err_gmb[ind] > PI;
                    let under = _err_gmb[ind] < PI;
                    _err_gmb[ind] = _err_gmb[ind] - (2.0 * PI * (over as u8 as f64))
                        + (2.0 * PI * (under as u8 as f64));
                }
            }
            // println!("wrapped gmb err: {}", _err_gmb);

            self.err_gmb_th
                .update_gimbal_coordinates(&[_err_gmb[0], _err_gmb[1], _err_gmb[2]]);

            if norm_fine_err > self.u_upper {
                // calculate the gimbal error
                self.err_comb_th = _err_gmb.clone();
                // println!("yaw error {}", self.err_comb_th[2]);
            } else {
                let err_weight: f64 =
                    (norm_fine_err - self.u_lower) / (self.u_upper - self.u_lower);

                let comb_err = ((1.0 - err_weight) * (phi * _err_gmb)) + (err_weight * _err_gmb);

                self.err_comb_th = comb_err;
                // println!("error weight: {}, error term 1: {}, error term 2: {}, gmb err: {}, fine err: {}", err_weight, ((1.0 - err_weight) * (phi*self.err_eq_th)), (err_weight * self.err_gmb_th), self.err_eq_th, self.err_gmb_th);
            }
            trace!("update_positional_pointing_error end");
        }
    }
}