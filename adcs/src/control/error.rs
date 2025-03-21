//! This module performs two main functions, first to calculate
//! the positional error in both the gimbal frame and in the
//! equatorial frame. These errors are then combined in a smooth
//!  function (to maintain stability) to use each inparticular
//! scenarios. The second function takes the positional error
//! and generates a desired velocity profile for each gimbal axis.
//!  This is all achieved using an Error struct with implementations.
extern crate nalgebra as na;
use std::f64::consts::{E, PI};
use std::collections::VecDeque;
extern crate statrs as stat;

use crate::control::state as st;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use na::Rotation;

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
    /// err rate lqr (gondola gyros)
    pub err_rate_lqr: na::Vector3<f64>,
    /// Integrated gimbal rate error (roll, pitch, yaw)
    pub err_rate_sum: na::Vector3<f64>,
     /// Integrated gimbal fine error (roll, pitch, yaw)
     pub err_fine_sum: na::Vector3<f64>,
    /// desired gimbal rates (roll, pitch, yaw)
    pub rate_des: na::Vector3<f64>,
    /// Integrated rate error time constant of decay
    pub _err_tc_v: f64,
    /// Decay constant for rate error integral, function of control rate and decay time constant
    pub _err_decay_v: f64,
    /// Integrated fine error time constant of decay
    pub _err_tc_p: f64,
    /// Decay constant for position error integral, function of control rate and decay time constant
    pub _err_decay_p: f64,
    /// Control algorithm time step
    pub _ctrl_dt: f64,
    /// Positional error in rotation matrix form
    pub rot_err: na::Rotation3<f64>,
    /// Lower bound constant for combined positional error calculation
    pub u_lower: f64,
    /// Upper bound constant for combined positional error calculation
    pub u_upper: f64,
    /// the mapped gimbal rates from omega decomposed using the gmm_i
    pub _d_theta: na::Vector3<f64>,
    /// scale factor for tuning the speed profile
    pub scale: f64,
    /// blended velocity scale 
    pub scale_blend: f64,
    /// the blend scalar describing where in between fine pointing and slewing
    pub err_weight: f64,
    
    pub amplitude: f64,
    pub period: f64,
    pub max_err: f64,
    pub max_err_old: f64,
    pub current_time: f64,
    pub phase_shift: f64,
}

impl Error {
    /// This function generates a new instance of the Error struct with default values
    pub fn new() -> Error {
        let _err_tc_v: f64 = 1.0;
        let _err_tc_p: f64 = 80.0;
        let _ctrl_dt: f64 = 0.01;
        let _err_decay_v: f64 = E.powf(-_ctrl_dt / _err_tc_v);
        let _err_decay_p: f64 = E.powf(-_ctrl_dt / _err_tc_p);
        let err_gmb_th = st::gimbal::Gimbal::new();
        let err_b_th = st::gimbal::Gimbal::new();
        let err_comb_th = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_rate = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_rate_lqr = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_rate_sum = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let err_fine_sum = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let rate_des = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let rot_err = na::Rotation3::<f64>::identity();
        let u_lower = 0.05;
        let u_upper = 0.25;
        let _d_theta = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let scale = 2.25;
        let scale_blend = 1.0;
        let err_weight = 1.0;
        let amplitude = 1.0;
        let period = 30.0;
        let max_err = 0.0;
        let max_err_old = 0.0;
        let current_time = 0.0;
        let phase_shift = 0.0;

        let error: Error = Error {
            err_b_th,
            err_gmb_th,
            err_comb_th,
            err_rate,
            err_rate_lqr,
            err_rate_sum,
            err_fine_sum,
            rate_des,
            rot_err,
            _err_tc_v,
            _err_tc_p,
            _ctrl_dt,
            _err_decay_v,
            _err_decay_p,
            u_lower,
            u_upper,
            _d_theta,
            scale,
            scale_blend,
            err_weight,
            amplitude,
            period,
            max_err,
            max_err_old,
            current_time,
            phase_shift,

        };
        error
    }

    pub fn update_decay(&mut self, _err_tc_v: &f64, _err_tc_p: &f64){
        self._err_decay_v = E.powf(-self._ctrl_dt / _err_tc_v);
        self._err_decay_p = E.powf(-self._ctrl_dt / _err_tc_p);
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
    pub fn update_pointing_positional_error(&mut self, state: &st::State, phi: na::Matrix3<f64>, mut slew_flag:bool) -> (bool) {
        trace!("update_positional_pointing_error start");
        // let phi = phi.try_inverse().unwrap();
        //rot err in equatorial space
        let _rot_err = state.eq_d.rot * state.eq_k.rot.inverse();

        // parameterize the rotation error into euler angles
        let _rot_err_eul = _rot_err.inverse().euler_angles();
        
        let err_vec = na::Vector3::<f64>::new(_rot_err_eul.0, _rot_err_eul.1, _rot_err_eul.2);
        // println!("des ra {}, dec {}, fr {}", &state.eq_d.ra, &state.eq_d.dec, &state.eq_d.fr);
        // println!("cur ra {}, dec {}, fr {}", &state.eq_k.ra, &state.eq_k.dec, &state.eq_k.fr);
        
        // println!("eq_des: {:.7}, eq_k: {:.7}, eq_err: {:.7}", &state.eq_d.rot, &state.eq_k.rot, &_rot_err);
        // println!("gps time {}, gps lat {}, gps lon {}", &state.gps.utc, &state.gps.lat, &state.gps.lon);
        // println!("desired gimbal roll {}, pitch {}, yaw {}", &state.gmb_d.roll, &state.gmb_d.pitch, &state.gmb_d.yaw);
        // println!("current gimbal roll {}, pitch {}, yaw {}", &state.gmb_k.roll, &state.gmb_k.pitch, &state.gmb_k.yaw);

        let norm_fine_err = err_vec.norm();

        let _err_tc_v = (2.0 + (1.0/ (1.0+E.powf(60.0*&norm_fine_err.powf(2.0))))).powf(1.0);
        // println!("err_decay{}", &_err_tc_v);

        self._err_decay_p = E.powf(-self._ctrl_dt / _err_tc_v);
        self._err_decay_v = E.powf(-self._ctrl_dt / _err_tc_v);

        // println!("norm fine error: {}", norm_fine_err);

        let mapped_err = state.gmb_k.gmm_i * err_vec;
        // let mapped_err = err_vec;
        // let mut state_d = state.gmb_d.clone();
        // state_d.calculate_gimbal_mapping_matrix();
        // state_d.calculate_inverse_gimbal_mapping_matrix();
        // let mapped_err = state_k.gmm_i * err_vec;

        self.err_b_th.roll = mapped_err[0];
        self.err_b_th.pitch = mapped_err[1];
        self.err_b_th.yaw = mapped_err[2];

        // println!("mapped err {:?}", &mapped_err);

        let mut _err_gmb = na::Vector3::<f64>::new(
            state.gmb_d.roll - state.gmb_k.roll,
            state.gmb_d.pitch - state.gmb_k.pitch,
            state.gmb_d.yaw - state.gmb_k.yaw,
        );


        // println!("gmb err {:?}", &_err_gmb);
        for ind in 0..3 {
            if _err_gmb[ind].abs() > PI {
                let over = _err_gmb[ind] > PI;
                let under = _err_gmb[ind] < PI;
                _err_gmb[ind] = _err_gmb[ind] - (2.0 * PI * (over as u8 as f64))
                    + (2.0 * PI * (under as u8 as f64));
            }
        }
        let err_mat_phi = phi*err_vec;
        //decay the fine suym error regarless of error but only integrate it if error is in fine error bracket
        self.err_fine_sum = self.err_fine_sum * self._err_decay_p;

        // if norm_fine_err < self.u_lower {
        if norm_fine_err < self.u_lower {
            // self.err_comb_th = phi * err_vec;
            self.err_comb_th = phi * mapped_err;
            self.err_fine_sum = self.err_fine_sum + (self.err_comb_th * self._ctrl_dt);
            slew_flag = false;
            self.err_weight = 0.0;
            self.scale_blend = 1.0;
            // println!("{}\n {} \n{}\n", _err_gmb.x, _err_gmb.y,_err_gmb.z );
            
        } else {
            
            self.err_gmb_th
                .update_gimbal_coordinates(&[_err_gmb[0], _err_gmb[1], _err_gmb[2]]);

            if norm_fine_err > self.u_upper {
                // calculate the gimbal error
                if _err_gmb.z.abs() > 0.015{
                    _err_gmb.x = 0.0;
                    _err_gmb.y = 0.0;
                }
                self.err_comb_th = _err_gmb.clone();
                slew_flag = true;
                self.err_weight = 1.0;
                self.scale_blend = 1.0;

                
                // println!("greater than u.upper {}", &norm_fine_err);
            } else {
                
                // println!("between u.lower and u.upper {}", &norm_fine_err);
                let err_weight: f64 =
                    (norm_fine_err - self.u_lower) / (self.u_upper - self.u_lower);

                self.err_weight = err_weight;

                // let comb_err = ((1.0 - err_weight) * (phi * err_vec)) + (err_weight * _err_gmb);
                let comb_err = ((1.0 - err_weight) * (phi * mapped_err)) + (err_weight * _err_gmb);
                self.err_fine_sum = self.err_fine_sum + (((1.0 - err_weight) * (phi * err_vec)) * self._ctrl_dt);

                // self.scale_blend = ((1.0 - err_weight) * (self.scale)) + (err_weight * 0.05);
                self.scale_blend = 1.0;

                self.err_comb_th = comb_err;
                // println!("error weight: {}, combined error: {}\n term1 {} term2 {}", err_weight, comb_err, (1.0 - err_weight) * (phi * err_vec), (err_weight * _err_gmb));
            }
            // self.err_comb_th = _err_gmb.clone();
            trace!("update_positional_pointing_error end");
            // println!("yaw gmb {:.4} vec {:.4} comb {:.4}",_err_gmb, err_vec, self.err_comb_th);
        }
        return slew_flag
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

        let ceh = state.ceh;
        let w_di = na::Vector3::<f64>::new(0.0, 0.0, 15.04108/206265.0);
        let w_di_h = ceh.transpose() * w_di;

        state.gmb_k.calculate_rotation_matrix();
        let c7h = state.gmb_k.calculate_rotation_matrix_7h();
        let c8h = state.gmb_k.calculate_rotation_matrix_8h();

        let wdi_7 = c7h * w_di_h;
        let wdi_8 = c8h * w_di_h;
        let wdi_9 = state.gmb_k.rot * w_di_h;



        // println!("7: {} 8: {} 9: {}", wdi_7, wdi_8, wdi_9);
        // println!("c8h: {}, roll {}, c7h {}, yaw {}, cbh {} pitch {}", c8h, state.gmb_k.roll, c7h, state.gmb_k.yaw, state.gmb_k.rot, state.gmb_k.pitch);

        // println!("{:?}", state.gmb_k);
        let _d_theta: na::Vector3<f64> = state.gmb_k.gmm_i * state.omega;

        self._d_theta = _d_theta;

        let mut _roll_rate_des: f64 = 0.0;
        let mut _pitch_rate_des: f64 = 0.0;
        let mut _yaw_rate_des: f64 = 0.0;     

        let max_v = 0.005;
        let slope = 100.0;

        _roll_rate_des = max_v*self.err_comb_th.x.signum()
                    * stat::function::erf::erf(self.err_comb_th.x.abs() * slope) - 0.0*wdi_8[0];      

        let max_v = 0.005;
        let slope = 100.0*1.0;

        _pitch_rate_des = max_v*self.err_comb_th.y.signum()
                * stat::function::erf::erf(self.err_comb_th.y.abs() * slope) + 0.0*wdi_9[1];

        let max_v = 0.008/2.0;
        let slope = 15000.0/0.5;


        _yaw_rate_des = max_v*self.err_comb_th.z.signum()
                * stat::function::erf::erf(self.err_comb_th.z.abs() * slope) - wdi_7[2];

        self.rate_des = na::Vector3::new(_roll_rate_des, _pitch_rate_des, _yaw_rate_des);

        // Based on boresight gyros
        let mut err_gmb_rate: [f64; 3] = [
            (_roll_rate_des - _d_theta.x),
            (_pitch_rate_des - _d_theta.y),
            (_yaw_rate_des - _d_theta.z),
        ];
        
        // Based on the gimbal gyros
        // have to subtract the yaw motion that is read by roll
        let omega_yfp = state.omega_gond[0]*state.gmb_k.roll.sin();

        // add in diurnal motion
        
        let err_gmb_rate2: [f64; 3] = [
            (_roll_rate_des - state.omega_gond[1]),
            (_pitch_rate_des - (state.omega_gond[3] - omega_yfp)),
            (_yaw_rate_des - state.omega_gond[0]),
        ];

        self.err_rate_lqr =na::Vector3::<f64>::from_row_slice(&err_gmb_rate2);


        self.err_rate = na::Vector3::<f64>::from_row_slice(&err_gmb_rate);
        trace!("update_pointing_velocity_error_terms end");

        self.err_rate_sum = (self.err_rate_sum * self._err_decay_v) + (self.err_rate * self._ctrl_dt)
    }
}
