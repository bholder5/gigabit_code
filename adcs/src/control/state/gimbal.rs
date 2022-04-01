//! Outlines a struct and implements functions associated 
//! specifically with handling gimbal coordinates for use
//!  in the ADCS Pointing system.

// import crate for matrix math
extern crate nalgebra as na;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
pub struct Gimbal {
    /// Roll, rotation of the middle frame in a 312 euler sequence
    pub roll: f64,
    /// Pitch, rotation of the inner frame in the 312 euler sequence
    pub pitch: f64,
    /// Yaw, rotation of the outer frame in the 312 euler sequence
    pub yaw: f64,
    /// Rotation matrix associated with the yaw roll pitch (312) euler sequence
    pub rot: na::Rotation3<f64>,
    /// the gimbal mapping matrix, maps from boresight telescope frame to 
    /// gimbal coordinates ($\dot\theta = \Phi_{gmm}\omega$)
    pub gmm: na::Matrix3<f64>,
}

impl Gimbal {
    /// Function to update the gimbal orientation rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the rotation matrix based on 
    /// the gimbal 312 euler parameters currently in the 
    /// gimbal struct according to: 
    /// $$ C_{b,G} = C_{y}(pitch)C_{x}(roll)C_{z}(yaw)$$
    /// 
    /// # Arguments
    ///
    /// - `roll: f64` - Roll
    /// - `pitch: f64` - Pitch
    /// - `yaw: f64` - Yaw
    ///
    /// # Results
    ///
    /// - `self.rot: na::Matrix3<f64>` - $C_{b,G}$ the orientation rotation matrix
    pub fn calculate_rotation_matrix(&mut self){
        trace!("calculate_rotation_matrix start");
        let rotx = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::x_axis(), self.roll).inverse();
        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), self.pitch).inverse();
        let rotz = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::z_axis(), self.yaw).inverse();

        // calculating the CbH using the gimbal euler 312 rotation.
        self.rot = roty * rotx * rotz;
        // println!("gimbal rot: {}", self.rot);
        trace!("calculate_rotation_matrix end");
    }
    /// Function to update the gimbal euler angles directly from a vector
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the yaw pitch and roll and subsequently 
    /// updates the rotation matrix based on the new parameters. 
    /// The vector contains the euler angles in the following order 
    /// - $$ Vec = [roll, pitch, yaw]^T$$
    ///
    /// # Arguments
    /// - `vec: &[f64; 3]` - the new equatorial euler angles
    /// 
    /// # Results
    /// - `roll: f64` - Roll
    /// - `pitch: f64` - Pitch
    /// - `yaw: f64` - Yaw
    /// - `self.rot: na::Rotation3<f64>` - orientation rotation matrix
    pub fn update_gimbal_coordinates(&mut self, vec: &[f64; 3]){
        trace!("update_gimbal_coordinates start");
        self.roll  = vec[0];
        self.pitch = vec[1];
        self.yaw  = vec[2];

        self.calculate_rotation_matrix();
        trace!("update_gimbal_coordinates end");
    }
    /// Function to parameterize the Yaw Pitch and Roll from a rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the ra dec and field rotation from
    /// a rotation matrix ($C_G$) according to 
    /// - $$ Pitch = \atan2(-C_G(1,3), C_G(3,3))$$
    /// - $$ Roll = \asin(C_G(2,3))$$
    /// - $$ Yaw = \atan2(-C_G(2,1), C_E(2,2))$$
    ///
    /// # Arguments
    /// - `vec: &[f64; 3]` - the new equatorial euler angles
    /// 
    /// # Results
    /// - `roll: f64` - Roll
    /// - `pitch: f64` - Pitch
    /// - `yaw: f64` - Yaw
    pub fn extract_gimbal_rpy(&mut self){
        trace!("extract_gimbal_rpy start");
        self.roll = self.rot[(1,2)].asin();
        self.pitch = (-self.rot[(0,2)]).atan2(self.rot[(2,2)]);
        self.yaw = (-self.rot[(1,0)]).atan2(self.rot[(1,1)]);
        // println!("mat: {} roll: {}, pitch {}, yaw {}", self.rot,self.roll, self.pitch, self.yaw);
        trace!("extract_gimbal_rpy end");
    }
    /// Function to calculate the gimbal mapping matrix for mapping joint rates to angular velocity
    /// 
    /// # Detailed Explanation
    /// 
    /// This function calculates the gimbal mapping matrix for use in the control system. This 
    /// matrix is necessary for determining joint torques necessary to correct for deviations in telescope
    /// angular velocity. The gimbal mapping matrix $S(\theta)$ can be calculated (for the 312 euler 
    /// sequence of the bit gondola) according to
    ///  $$ S(\theta) = C_{x}()$$
    pub fn calculate_gimbal_mapping_matrix(&mut self){
        trace!("calculate_gimbal_mapping_matrix start");
        let rotx = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::x_axis(), self.roll).inverse();
        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), self.pitch).inverse();

        let col1 = roty * na::Vector3::<f64>::new(1.0, 0.0, 0.0);
        let col2 = na::Vector3::<f64>::new(0.0, 1.0, 0.0);
        let col3 = roty * rotx * na::Vector3::<f64>::new(0.0, 0.0, 1.0);
        let gmm = na::Matrix3::<f64>::from_columns(&[col1, col2, col3]);
        self.gmm = gmm;
        trace!("calculate_gimbal_mapping_matrix end");
    }

    /// This function creates a new instance of the Gimbal struct
    /// initialized with selected initial values and then calculates
    /// the corresponding rotation matrix
    pub fn new() -> Gimbal {
        let roll: f64 = -0.0259154632;
        let pitch: f64 = -0.750145461;
        let yaw: f64 = -0.0278;
        let rot = na::Rotation3::<f64>::identity();
        let gmm = na::Matrix3::<f64>::identity();
        let mut gimbal = Gimbal{
            roll,
            pitch,
            yaw,
            rot,
            gmm,
        };
        gimbal.calculate_rotation_matrix();
        gimbal
    }
}