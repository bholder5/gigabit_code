//! Contains a struct and functions for dealing with 
//! horizontal coordinates
// import crate for matrix math
extern crate nalgebra as na;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
/// Struct to hold the horizontal coordinates
pub struct Horizontal {
    /// Azimuth angle radians
    pub az: f64,
    /// Elevation angle radians
    pub el: f64,
    /// Image rotation angle radians
    pub ir: f64,
    /// Horizontal frame rotation matrix
    pub rot: na::Rotation3<f64>,
}

impl Horizontal {
    /// Function to update the horizontal rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the rotation matrix based on 
    /// the horizontal euler parameters currently in the 
    /// horizontal struct according to: 
    /// $$ C_{b,H} = C_{x}(ir)C_{y}(-dec)C_{z}(-az)$$
    ///
    /// # Arguments
    ///
    /// - `ir: f64` - Image rotation
    /// - `el: f64` - Elevation
    /// - `az: f64` - Azimuth
    ///
    /// # Results
    ///
    /// - `self.rot: na::Matrix3<f64>` - $C_{b,H}$ the orientation rotation matrix
    pub fn calculate_rotation_matrix(&mut self){
        self.rot = na::Rotation3::from_euler_angles(self.ir, -self.el, -self.az).inverse();
    }

    /// Function to update the horizontal euler angles directly from a vector
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the az el and ir and subsequently 
    /// updates the rotation matrix based on the new parameters. 
    /// The vector contains the euler angles in the following order 
    /// - $$ Vec = [ir, el, az]^T$$
    ///
    /// # Arguments
    /// - `vec: &[f64; 3]` - the new horizontal euler angles
    /// 
    /// # Results
    /// - `ir: f64` - Image rotation
    /// - `el: f64` - Elevation
    /// - `az: f64` - Azimuth
    pub fn update_horizontal_coordinates(&mut self, vec: &[f64; 3]){
        self.ir = vec[0];
        self.el = vec[1];
        self.az = vec[2];

        self.calculate_rotation_matrix();
    }

    /// Function to parameterize the Az El Ir from a rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the Az El and ir from
    /// a rotation matrix ($C_H$) according to 
    /// - $$ Ir = \atan2(C_H(2,3), C_H(3,3))$$
    /// - $$ El = \asin(C_H(1,3))$$
    /// - $$ Az = \atan2(-C_H(1,2), C_H(1,1))$$
    ///
    /// # Arguments
    /// - `rot: na::Rotation3<f64?` - The horizontal to body rotation matrix
    /// 
    /// # Results
    /// - `ir: f64` - Image rotation
    /// - `el: f64` - Elevation
    /// - `az: f64` - Azimuth
    pub fn extract_az_el_ir_from_rotmat(&mut self){
        debug!("cbh{:}", self.rot);

        let az = (-self.rot[(0,1)]).atan2(self.rot[(0,0)]);
        let el = self.rot[(0,2)].asin();
        let ir = self.rot[(1,2)].atan2(self.rot[(2,2)]);

        self.az = az.clone();
        self.ir = ir.clone();
        self.el = el.clone();
        debug!("az {}, el {}, ir {}", az.to_degrees(),el.to_degrees(),ir.to_degrees());
    }

    /// This function creates a new instance of the Horizontal struct
    /// initialized to all $0$s which has the corresponding identity
    /// for a rotation matrix
    pub fn new() -> Horizontal {
        let az = 0.0;
        let el = 0.0;
        let ir = 0.0;
        let rot = na::Rotation3::<f64>::identity();
        let hor = Horizontal{
            az,
            el,
            ir,
            rot,
        };
        hor
    }
}