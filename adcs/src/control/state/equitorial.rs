//! Outlines a struct and implements functions associated
//! specifically with handling equatorial coordinates for
//! use in the ADCS Pointing system.

// import crate for matrix math
extern crate nalgebra as na;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
/// A Struct containing an equatorial orientation in euler form and rotation matrix form
pub struct Equatorial {
    /// Right Ascension
    pub ra: f64,
    /// Declination
    pub dec: f64,
    /// Field Rotation
    pub fr: f64,
    /// Rotation matrix corresponding to the ra, dec and fr euler parameters
    pub rot: na::Rotation3<f64>,
}

impl Equatorial {
    /// Function to update the equatorial rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the rotation matrix based on
    /// the equatorial euler parameters currently in the
    /// equatorial struct according to:
    /// $$ C_{b,E} = C_{x}(fr)C_{y}(-dec)C_{z}(ra)$$
    ///
    /// # Arguments
    ///
    /// - `ra: f64` - Right ascension
    /// - `dec: f64` - Declination
    /// - `fr: f64` - Field rotation
    ///
    /// # Results
    ///
    /// - `self.rot: na::Matrix3<f64>` - $C_{b,E}$ the orientation rotation matrix
    pub fn calculate_rotation_matrix(&mut self) {
        trace!("calculate_rotation_matric_eq start");
        self.rot = na::Rotation3::from_euler_angles(self.fr, -self.dec, self.ra).inverse();
        trace!("calculate_rotation_matric_eq end");
    }

    /// Function to update the equatorial euler angles directly from a vector
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the ra dec and field rotation and subsequently
    /// updates the rotation matrix based on the new parameters.
    /// The vector contains the euler angles in the following order
    /// - $$ Vec = [fr, dec, ra]^T$$
    ///
    /// # Arguments
    /// - `vec: &[f64; 3]` - the new equatorial euler angles
    ///
    /// # Results
    /// - `fr: f64` - Field rotation
    /// - `dec: f64` - Declination
    /// - `ra: f64` - Right ascension
    pub fn update_equatorial_coordinates(&mut self, vec: &[f64; 3]) {
        trace!("update_equatorial_coordinates start");
        self.fr = vec[0];
        self.dec = vec[1];
        self.ra = vec[2];

        self.calculate_rotation_matrix();
        trace!("update_equatorial_coordinates end");
    }

    /// Function to parameterize the RA DEC and FR from a rotation matrix
    ///
    /// # Detailed Explanation
    ///
    /// This functions updates the ra dec and field rotation from
    /// a rotation matrix ($C_E$) according to
    /// - $$ FR = \atan2(C_E(2,3), C_E(3,3))$$
    /// - $$ DEC = \asin(C_E(1,3))$$
    /// - $$ RA = \atan2(C_E(1,2), C_E(1,1))$$
    ///
    /// # Arguments
    /// - `rot: na::Rotation3<f64>` - The Equatorial to body rotation matrix
    ///
    /// # Results
    /// - `fr: f64` - Field rotation
    /// - `dec: f64` - Declination
    /// - `ra: f64` - Right ascension
    pub fn extract_ra_dec_fr_from_rotmat(&mut self) {
        trace!("extract_ra_dec_fr_from_rotmat start");
        let ra = (self.rot[(0, 1)]).atan2(self.rot[(0, 0)]);
        let dec = self.rot[(0, 2)].asin();
        let fr = self.rot[(1, 2)].atan2(self.rot[(2, 2)]);

        self.ra = ra.clone();
        self.dec = dec.clone();
        self.fr = fr.clone();
        trace!("extract_ra_dec_fr_from_rotmat end");
    }

    /// This function creates a new instance of the Equatorial struct
    /// initialized to all $0$s which has the corresponding identity
    /// for a rotation matrix
    pub fn new() -> Equatorial {
        let ra = 0.0;
        let dec = 0.0;
        let fr = 0.0;
        let rot = na::Rotation3::<f64>::identity();
        let eq = Equatorial { ra, dec, fr, rot };
        eq
    }
}
