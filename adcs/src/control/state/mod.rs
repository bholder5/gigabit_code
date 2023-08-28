//! The State module implements a struct of other modules
//! usefule for cleanly organizing all necessary information
//! and functions required in the ADCS system that has to do
//! with the current orientation of the instrument

// import crate for matrix math
extern crate nalgebra as na;

// Declare the sub-modules for the state module for organizational purposes
pub mod equitorial;
pub mod gimbal;
pub mod gps;
pub mod horizontal;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use {equitorial as eq, gimbal as gb, horizontal as hor};

#[derive(Debug, Clone)]
pub struct State {
    /// The angular velocity of the telescope boresight
    pub omega: na::Vector3<f64>,
    /// The angular velocity of the telescope boresight
    pub omega_gond: na::Vector5<f64>,
    /// The joint velocity of the gimbal yaw pitch roll in 321 order
    pub dtheta: na::Vector3<f64>,
    /// The longitude, latitude and sidereal time
    pub gps: gps::Gps,
    /// Current equatorial coordinates
    pub eq_k: eq::Equatorial,
    /// Desired equatorial coordinates
    pub eq_d: eq::Equatorial,
    /// Current horizontal coordinates (for telemetry purposes only)
    pub hor: hor::Horizontal,
    /// Desired Horizontal coordinates
    pub hor_d: hor::Horizontal,
    /// Current gimbal coordinates
    pub gmb_k: gb::Gimbal,
    /// Desired gimbal coordinates
    pub gmb_d: gb::Gimbal,
    /// eqautorial to horizontal mapping rotation matrix
    pub ceh: na::Rotation3<f64>,
    /// offset from imbalance between horizontal frame and body fized frame
    _b2h_offset: na::Rotation3<f64>,
    pub piv_angle: f64,
    pub piv_speed: f64,
    pub rw_hs_nom: f64,
    pub rw_hs: f64,
}

impl State {
    /// Function to instantiate a new state type
    ///
    /// # Detailed Explanation
    ///
    /// Initializes each sub struct/item to zeros or the
    /// described new implementation of that struct
    ///
    /// # Results
    ///
    /// - `returns state::state::State` -
    pub fn new() -> State {
        let omega = na::Vector3::<f64>::zeros();
        let omega_gond = na::Vector5::<f64>::zeros();
        let dtheta = na::Vector3::<f64>::zeros();
        let gps = gps::Gps::new();
        let eq_k = eq::Equatorial::new();
        let eq_d = eq::Equatorial::new();
        let hor = hor::Horizontal::new();
        let hor_d = hor::Horizontal::new();
        let gmb_k = gb::Gimbal::new();
        let gmb_d = gb::Gimbal::new();
        let ceh = na::Rotation3::<f64>::identity();
        let _b2h_offset = na::Rotation3::<f64>::identity();
        let piv_angle = 0.0;
        let piv_speed: f64 = 0.0;
        let rw_hs_nom = 0.0;
        let rw_hs: f64 = 0.0;

        let state: State = State {
            omega,
            omega_gond,
            dtheta,
            gps,
            eq_k,
            eq_d,
            hor,
            hor_d,
            gmb_k,
            gmb_d,
            ceh,
            _b2h_offset,
            piv_angle,
            piv_speed,
            rw_hs_nom,
            rw_hs,
        };
        state
    }

    /// Calculates the "coupling" matrix for an euler angle system
    ///
    /// # Detailed explanation
    /// This function calculates the coupling matrix which essentially
    ///  describes how each euler angle interacts with each other.
    /// More technically, this matrix is calculated as
    /// $$\Phi_{\theta} = (S_{\theta}^TS_{\theta})^{-1}$$, where
    /// $S_{\theta}$ is a matrix that maps the nominal gimbal (frame)
    /// angular rates to the telescope angular velocity It turns out,
    /// $$\Phi_{\theta} = \left[ \begin{matrix} 1 & 0 & 0\\ 0 & 1 &
    /// \sin(\theta_{roll}) \\ 0 & \sin(\theta_{roll} & 1) \end{matrix}\right]$$
    ///
    ///  # Arguments
    ///
    /// `self.gmb_k.roll` - the roll angle of the bit gondola, passed from the estimation algorithm
    ///
    /// # Result
    ///
    /// `_phi` - the coupling matrix for use in the control algorithm.
    pub fn calculate_coupling_matrix(&mut self) -> na::Matrix3<f64> {
        trace!("calculate_coupling_matrix start");
        let st2: f64 = self.gmb_k.roll.sin();
        let _phi = na::Matrix3::<f64>::from_row_slice(&[1., 0., 0., 0., 1., st2, 0., st2, 1.]);
        trace!("calculate_coupling_matrix end");
        _phi
    }
    /// Function to update the horizontal to equatorial conversion rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This function uses the local apparent sidereal time of the payload and
    /// the latitude of the payload to calculate the  equatorial wrt horizontal
    ///  rotation matrix as follows $$ C_{y}(\frac{\pi}{2} - \phi_{lat}) C_{z}(\pi + LAST) $$
    ///
    /// # Arguments
    ///
    /// `gps.gast: f64` - The greenwhich apparent sidereal time in degrees
    /// `gps.long: f64` - Payload longitude in degrees
    /// `gps.lat: f64` - Payload latitude in degrees
    ///
    /// # Results
    ///
    /// - `self.ceh` - horizontal wrt equatorial rotation matrix
    pub fn update_hor_to_eq_conversion(&mut self) {
        trace!("update_hor_to_eq_conversion start");
        self.gps.get_greenwhich_apparent_sidereal_time();

        let last = self.gps.gast + self.gps.lon;
        let arg_z = (last + 180.0).to_radians();

        let arg_y = (self.gps.lat - 90.0).to_radians();

        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), arg_y);
        let rotz = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::z_axis(), arg_z);

        // println!("arg y {:}, argz {:}", arg_y, arg_z);
        // println!("rot y {:}, rotz {:}, ceh {:}", roty, rotz, rotz*roty);
        self.ceh = rotz * roty;
        // println!("ceh: {}", self.ceh);
        trace!("update_hor_to_eq_conversion end");
    }
    /// Function to update the current equatorial coordinates
    /// from the current horizontal coordinates
    ///
    /// # Detailed Explanation
    ///
    /// This function takes in the current horizontal coordinates,
    /// rotates them and update the current equitorial coordinates
    ///  according to: $C_{b,E} = C_{b,H}C_{H,E}$ where $C_{b,E} is
    /// the orientation in the equatorial frame, and C_{b,H} is the
    ///  orientation in the horizontal frame.
    ///
    /// # Arguments
    ///
    /// `vec: &[f64;3]` - vector of the horizontal coordinates in [ir, el, az] order
    ///
    /// # Results
    ///
    /// - `eq_k: equitorial::Equatorial` - a fully updated current
    /// equatorial orientation
    pub fn update_current_equatorial_coordinates(&mut self, vec: &[f64; 3]) {
        trace!("update_current_equatorial_coordinates start");
        self.hor.ir = vec[0].clone();
        self.hor.el = vec[1].clone();
        self.hor.az = vec[2].clone();
        self.hor.calculate_rotation_matrix();

        self.update_hor_to_eq_conversion();
        self.eq_k.rot = self.hor.rot * self.ceh.inverse();
        // pre multiply seld.eq_k.rot by phi_

        self.eq_k.extract_ra_dec_fr_from_rotmat();
        trace!("update_current_equatorial_coordinates end");
    }

    /// Function to update the desired gimbal coordinates from the desired
    /// equatorial coordinates
    ///
    /// # Detailed Explanation
    ///
    /// This function converts the desired equatorial coordinates into the
    /// gimbal frame yaw pitch and roll through rotating the desired
    /// equatorial frame rotation matrix into the gimbal frame. NOTE:
    /// CURRENTLY NO OFFSET BETWEEN GIMBAL FRAME AND HORIZONTAL, NEED
    /// TO ADD. The conversion is done following
    ///  according to: $C_{b,G} = C_{b,E}C_{E,H}$ where $C_{b,E}$ is
    /// the orientation in the equatorial frame, and $C_{b,G}$ is the
    ///  orientation in the gimbal frame, $C_{E,H}$ is the conversion matrix
    /// using sidereal time and gps coordinates.
    ///
    /// # Arguments
    ///
    /// `eq_d.rot na::Rotation3<f64>` - Desired orientation rotation matrix
    /// in the equatorial frame
    ///
    /// # Results
    ///
    /// - `gmb_d: gimbal::Gimbal` - a fully updated desired gimbal rotation
    ///  matrix and corresponding euler angles.
    pub fn update_desired_gimbal_rpy(&mut self) {
        trace!("update_desired_gimbal_rpy start");
        self.update_hor_to_eq_conversion();
        self.gmb_d.rot = self.eq_d.rot * self.ceh;
        self.gmb_d.extract_gimbal_rpy();
        self.hor_d.rot = self.gmb_d.rot.clone();
        self.hor_d.extract_az_el_ir_from_rotmat();
        trace!("update_desired_gimbal_rpy end");
    }
    /// Function to update the desired equatorial coordinates from the desired
    /// gimbal coordinates
    ///
    /// # Detailed Explanation
    ///
    /// This function converts the desired gimbal coordinates into the
    /// equatorial frame through rotating the desired
    /// gimbal frame rotation matrix into the equatorial frame. NOTE:
    /// CURRENTLY NO OFFSET BETWEEN GIMBAL FRAME AND HORIZONTAL, NEED
    /// TO ADD. The conversion is done
    ///  according to: $C_{b,E} = C_{b,G}C_{H,E}$ where $C_{b,E}$ is
    /// the orientation in the equatorial frame, and $C_{b,G}$ is the
    ///  orientation in the gimbal frame, $C_{H,E}$ is the conversion matrix
    /// using sidereal time and gps coordinates.
    ///
    /// # Arguments
    ///
    /// `gmb_d gimbal::Gimbal` - Desired orientation in gimbal form
    ///
    /// # Results
    ///
    /// - `eq_d: equitorial::Equatorial` - a fully updated desired equatorial rotation
    ///  matrix and corresponding euler angles.
    pub fn update_desired_eq_from_gmb(&mut self) {
        trace!("update_desired_eq_rdf start");
        self.update_hor_to_eq_conversion();
        self.gmb_d.calculate_rotation_matrix();
        self.eq_d.rot = self.gmb_d.rot * self.ceh.inverse();
        self.eq_d.extract_ra_dec_fr_from_rotmat();
        // println!("desired: ra {} dec {} fr {}", self.eq_d.ra, self.eq_d.dec, self.eq_d.fr);
        trace!("update_desired_eq_rdf end");
    }

    /// Temp fucntion for testing error functions
    ///
    pub fn update_current_eq_from_gmb(&mut self) {
        trace!("update_current_eq_rdf start");
        self.update_hor_to_eq_conversion();
        self.gmb_k.calculate_rotation_matrix();
        self.eq_k.rot = self.gmb_k.rot * self.ceh;
        self.eq_k.extract_ra_dec_fr_from_rotmat();
        // println!("desired: ra {} dec {} fr {}", self.eq_d.ra, self.eq_d.dec, self.eq_d.fr);
        trace!("update_current_eq_rdf end");
    }
}
