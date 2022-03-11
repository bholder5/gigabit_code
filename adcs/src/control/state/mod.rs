// import crate for matrix math
extern crate nalgebra as na;

pub mod coords;
pub mod horizontal;
pub mod equitorial;
pub mod gps;
pub mod gimbal;

use {equitorial as eq, horizontal as hor, gimbal as gb};
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
pub struct State {
    pub omega: na::Vector3<f64>,
    // //
    // superbit pitch, roll, yaw
    /// bit_thet is in same order as simulation yaw, roll, pitch
    pub dtheta: na::Vector3<f64>,
    gps: gps::Gps,
    pub eq_k: eq::Equatorial,
    pub eq_d: eq::Equatorial,
    hor: hor::Horizontal,
    pub gmb_k: gb::Gimbal,
    pub gmb_d: gb::Gimbal,
    ceh: na::Rotation3<f64>,
    b2h_offset: na::Rotation3<f64>,
}

impl State {
    /// Function to update the orientation vector and rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This functions takes in a slice containing the euler angles in 123 sequence, and
    /// updates the orientation vector and the orientation matrix
    ///
    /// # Arguments
    ///
    /// `x: &[f64;3]` - Desired orientation vector
    ///
    /// # Results
    ///
    /// - `self.eul: na::Vector3<f64>` - the orientation vector
    /// - `self.rot: na::Matrix3<f64>` - the orientation matrix

    pub fn new() -> State{
        let omega = na::Vector3::<f64>::zeros();
        let dtheta = na::Vector3::<f64>::zeros();
        let gps = gps::Gps::new();
        let eq_k = eq::Equatorial::new();
        let eq_d = eq::Equatorial::new();
        let hor = hor::Horizontal::new();
        let gmb_k = gb::Gimbal::new();
        let gmb_d = gb::Gimbal::new();
        let ceh = na::Rotation3::<f64>::identity();
        let b2h_offset = na::Rotation3::<f64>::identity();

        let state: State = State {
            omega,
            dtheta,
            gps, 
            eq_k,
            eq_d,
            hor,
            gmb_k,
            gmb_d,
            ceh,
            b2h_offset,
        };
        state
    }

    /// Calculates the "coupling" matrix for an euler angle system
    ///
    /// # Detailed explanation
    /// This function calculates the coupling matrix which essentially describes how each euler angle interacts with each other. More technically, this matrix is calculated as
    /// $$\Phi_{\theta} = (S_{\theta}^TS_{\theta})^{-1}$$, where $S_{\theta}$ is a matrix that maps the nominal gimbal (frame) angular rates to the telescope angular velocity It turns out, $$\Phi_{\theta} = \left[ \begin{matrix} 1 & 0 & 0\\ 0 & 1 & \sin(\theta_{roll}) \\ 0 & \sin(\theta_{roll} & 1) \end{matrix}\right]$$
    ///  # Arguments
    /// `self.bit_thet[1]` - the roll angle of the bit gondola, passed from the estimation algorithm
    ///
    /// # Result
    ///
    /// `_phi` - the coupling matrix for use in the control algorithm.
    pub fn calculate_coupling_matrix(
        &mut self,
    ) -> na::Matrix<f64, na::U3, na::U3, na::ArrayStorage<f64, na::U3, na::U3>> {
        trace!("calculate_coupling_matrix start");
        let st2: f64 = self.gmb_k.roll.sin();
        let _phi = na::Matrix3::<f64>::from_row_slice(&[1., 0., 0., 0., 1., st2, 0., st2, 1.]);
        trace!("calculate_coupling_matrix end");
        _phi
    }
    pub fn update_hor_to_eq_conversion(&mut self){
        trace!("update_hor_to_eq_conversion start");
        self.gps.get_greenwhich_apparent_sidereal_time();

        let last = self.gps.gast + self.gps.lon;
        let arg_z = (last + 180.0).to_radians();

        let arg_y = (self.gps.lat - 90.0).to_radians();
        
        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), arg_y);
        let rotz = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::z_axis(), arg_z);

        // println!("arg y {:}, argz {:}", arg_y, arg_z);
        // println!("rot y {:}, rotz {:}, ceh {:}", roty, rotz, rotz*roty);
        self.ceh = rotz*roty;
        trace!("update_hor_to_eq_conversion end");
    }

    pub fn update_horizontal_coordinates(&mut self){
        trace!("update_horizontal_coordinates start");
        self.hor.rot = self.ceh * self.eq_k.rot;
        self.hor.extract_az_el_ir_from_rotmat();
        trace!("update_horizontal_coordinates end");
    }

    pub fn update_current_equatorial_coordinates(&mut self, vec: &[f64; 3]){
        trace!("update_current_equatorial_coordinates start");
        self.eq_k.update_equatorial_coordinates(vec);
        self.update_horizontal_coordinates();
        trace!("update_current_equatorial_coordinates end");
    }

    pub fn update_desired_gimbal_rpy(&mut self){
        trace!("update_desired_gimbal_rpy start");
        self.update_hor_to_eq_conversion();
        self.gmb_d.rot = self.ceh * self.eq_d.rot;
        self.gmb_d.extract_gimbal_rpy();
        trace!("update_desired_gimbal_rpy end");

    }
}