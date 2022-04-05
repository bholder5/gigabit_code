//! This module contains all of the initialization of BITs
//! parameters for use in the simulation
//!
//! provides the main struct formulation that carries all of
//! the parameters, defines matrix and vector types necessary,
//! then a function which is called to then initialize all of
//! the values. Changing BIT's parameters that have to do with the
//! bit simulation parameters happens here
//!
//! Examples would be:
//! -generating sensors measurements from the simulation state
//!

use crate::bindings::{compute_angular_velocity_C, compute_rotation_mat_C};
use adcs::control::Ctrl;
use adcs::estimation::Est;
// import crate for matrix math

extern crate nalgebra as na;
pub use na::{Matrix3, Matrix6, MatrixMN, Rotation3, Vector3, Vector6};
// use std::fmt;
use std::error::Error;
// use std::io;
use csv::WriterBuilder;
pub use std::env;
use std::f64::consts::PI;
pub use std::ffi::OsString;
pub use std::fs::{File, OpenOptions};
use std::io::BufReader;
pub use std::path::Path;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use serde::Deserialize;
use serde::Serialize;

// Create a matrix type
/// A stack-allocated, column-major, 9 dimensional vector
pub type Vector9 = na::MatrixMN<f64, na::U9, na::U1>;
/// A stack-allocated, column-major, 21 dimensional vector
pub type Vector18 = na::MatrixMN<f64, na::U18, na::U1>;
pub type Vector21 = na::MatrixMN<f64, na::U21, na::U1>;

/// number of degrees of freedom in the simulation (from balloon to telescope)
pub const N_DOF: u8 = 9;

// Defining a struct to hold the definitions of the system
// matrices
#[derive(Debug)]
pub struct Params {
    // //
    pub theta: Vector9,
    pub d_theta_dt: Vector9,
    pub dthet_thet: Vector18,
    pub x: Vector21,
    pub _zn: [[f64; 3]; 9],
    pub omega_rw: f64,
    pub _i_z_rw: f64,
    pub omega_m: Vector3<f64>,
    pub _dt: f64,
    pub _num_steps: u16,
    pub phi_act: [f64; 3],
    pub unlock: Vector9,
    pub piv_flag: bool,
    pub _tau_piv_max: f64,
    _roll_theta_max: f64,
    _pitch_theta_max: f64,
    _pitch_theta_min: f64,
    pub pitch_nom: f64,
    // //
    pub gps: Gps,
}

impl Params {
    /// Function to update the cureent and immediate past orientation rotation matricies
    ///
    /// # Details
    ///
    /// This function will update C0 to be C1, the same as saying $C_{k-1} = C_{k}$ in order o increase the index k on each step. Next the new/updated $C_k$
    /// is updated by calling the `compute_rotation_mat_C` from the bindings.rs crate (generated from bindgen using the output of the matlab code library functions)
    ///
    /// # Arguments
    /// `self.theta` - current angle of the 9 joints of the superbit simulation vector
    /// `self.C1` - the orientation rottion matrix of the telescope at the last time step
    ///
    /// # Results
    ///
    ///`self.C1` - the new current orientation rotation matrix
    pub fn get_omega_meas(&mut self) -> () {
        trace!("get_omega_meas start");
        let mut om: [f64; 3] = [0.; 3];
        unsafe {
            trace!("compute_angular_velocity_C start");
            compute_angular_velocity_C(
                self.dthet_thet.as_ptr(),
                self._zn.as_mut_ptr(),
                om.as_mut_ptr(),
            );
            trace!("compute_angular_velocity_C end");
        }
        self.omega_m = Vector3::<f64>::from_row_slice(om.as_ref());
        debug!("omega_m: {:}, om(pointer): {:?}", self.omega_m, om.as_ref());
        trace!("get_omega_meas end");
    }

    pub fn update_state(&mut self) -> () {
        trace!("update_state start");
        self.theta = Vector9::from_row_slice(&self.x.as_slice()[9..18]);
        self.d_theta_dt = Vector9::from_row_slice(&self.x.as_slice()[0..9]);
        self.bound_gondola_position();
        self.dthet_thet = Vector18::from_row_slice(&self.x.as_slice()[0..18]);
        debug!(
            "theta: {:}, state_array {:}, dtheta {:}, dthet_thet: {:}",
            self.theta, self.x, self.d_theta_dt, self.dthet_thet
        );
        trace!("update_state end");
    }
    /// get phi actual in equatorial frame, thus the declination is -ve
    pub fn get_orientation_vec(&mut self) -> () {
        // /////////////////////////////////////
        // this function is verified vs matlab
        // /////////////////////////////////////

        trace!("get_orientation_vector start");
        // let mut c_vec :[f64; 9] = [0.0; 9];
        let mut c_vec: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        unsafe {
            compute_rotation_mat_C(
                self._zn.as_mut_ptr(),
                self.theta.as_ptr(),
                c_vec.as_mut_ptr(),
            )
        }
        let c = c_vec.as_ref();
        // println!("{:?}", [c[0], c[1], c[2]].concat());
        let rotmat = Rotation3::<f64>::from_matrix(&Matrix3::<f64>::from_row_slice(
            &[c[0], c[1], c[2]].concat(),
        ));
        debug!("rotmat {:} from thet: {:}", rotmat, self.theta);
        let eul = rotmat.euler_angles();
        self.phi_act = [eul.0, -eul.1, -eul.2];
        debug!("x_act {:?}", self.phi_act);
        trace!("get_orientation_vector end");
    }

    pub fn get_rw_speed(&mut self) {
        trace!("get_rw_speed start");
        self.omega_rw = self.x[20] / self._i_z_rw;
        trace!("get_rw_speed end");
        // println!("Curent RW speed: {:}", self.omega_rw);
    }

    fn bound_gondola_position(&mut self) {
        trace!("bound_gondola_position start");

        //bound roll to +/- roll max
        if self.theta[7].abs() >= self._roll_theta_max {
            self.theta[7] = self.theta[7].signum() * self._roll_theta_max;
            self.x[16] = self.theta[7];
            if self.d_theta_dt[7].signum() == self.theta[7].signum() {
                self.d_theta_dt[7] = 0.0;
                self.x[7] = 0.0;
            }
        }

        // //bound pitch to pitch max
        // if self.theta[8] >= self._pitch_theta_max{
        //     self.theta[8] = self._pitch_theta_max;
        //     self.x[17] = self._pitch_theta_max;
        //     if self.d_theta_dt[8].signum() > 0.0{
        //         self.d_theta_dt[8] = 0.0;
        //         self.x[8] = 0.0;
        //     }
        // }
        // //bound pitch to pitch min
        // if self.theta[8] <= self._pitch_theta_min{
        //     self.theta[8] = self._pitch_theta_min;
        //     self.x[17] = self._pitch_theta_min;
        //     if self.d_theta_dt[8].signum() < 0.0{
        //         self.d_theta_dt[8] = 0.0;
        //         self.x[8] = 0.0;
        //     }
        // }
        trace!("bound_gondola_position end");
    }
}
// impl fmt::Display for Params {
//     // This trait requires `fmt` with this exact signature.
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         write!(f, "g0: {}, num dof: {}, Offset Vectors: {}, Rot Axes: {}, P_n: {}, Telescope Offset: {}, Joint Masses: {}, Joint COM: {}, MOI Matrices: {:?}", self._g0, self._n_dof, self._r_n1_n, self._z_n, self._p_n, self._tel_offset, self._m_n, self._c_n, self._i_n)
//     }
// }
///

pub fn init_bit() -> Params {
    //Set the inertial frame gravity vector at float

    // DEFINE INITIAL CONDITIONS
    // ----------------------------------
    let theta = Vector9::from_row_slice(&[
        0.,
        0.0 * 0.4 * PI / 180.0,
        0.0 * 0.4 * PI / 180.0,
        0.0 * 0.1 * PI / 180.0,
        0.0 * 0.1 * PI / 180.0,
        0.0 * 0.1 * PI / 180.0,
        0.0,
        0.0,
        1.0 * -40.0 * PI / 180.0,
    ]); // joint angles [rad]
    let d_theta_dt = Vector9::from_row_slice(&[0., 0., 0., 0., 0., 0., 0., 0., 0.]); // joint angle rates [rad/s]
    let dthet_thet = Vector18::from_row_slice(&[d_theta_dt.as_slice(), theta.as_slice()].concat());

    let omega_rw: f64 = 0.0;
    let _i_z_rw: f64 = 4.5;

    let mut x = Vector21::zeros();
    for _step in 0..N_DOF as usize {
        x[_step] = d_theta_dt[_step];
        x[_step + 9] = theta[_step];
    }
    // x[20] = PI * _i_z_rw;
    // let x0 = Vector21::from_row_slice(&d_theta_dt_0, &theta_0, &h_rw_0);

    let _zn = [
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];

    // [[0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
    // [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    // [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]]);

    let omega_m = Vector3::<f64>::zeros();
    let _dt = 0.0001;
    let _num_steps: u16 = 10;
    let phi_act = [0.0; 3];
    let unlock = Vector9::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let _tau_piv_max: f64 = 20.0;
    let _roll_theta_max: f64 = 6.0 * PI / 180.0;
    let _pitch_theta_max: f64 = -20.0 * PI / 180.0;
    let _pitch_theta_min: f64 = -60.0 * PI / 180.0;
    let pitch_nom: f64 = -40.0 * PI / 180.0;
    let gps = Gps::new();
    // define params struct_
    let _params = Params {
        theta,
        d_theta_dt,
        dthet_thet,
        x,
        _zn,
        omega_rw,
        _i_z_rw,
        omega_m,
        _dt,
        _num_steps,
        phi_act,
        unlock,
        piv_flag: true,
        _tau_piv_max,
        _roll_theta_max,
        _pitch_theta_max,
        _pitch_theta_min,
        pitch_nom,
        gps,
    };
    info!("Bit initialized: \n {:?}", _params);
    return _params;
}

#[derive(Debug, Deserialize, Default)]
pub struct Gps {
    lat: f64,
    lon: f64,
    utc: f64,
}

impl Gps {
    pub fn new() -> Gps {
        let lat = 47.002282304922659;
        let lon = -82.210336016074365;
        let utc = 946684800.0;

        let gps = Gps { lat, lon, utc };
        gps
    }
}