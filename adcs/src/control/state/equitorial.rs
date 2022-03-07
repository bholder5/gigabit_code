// import crate for matrix math
extern crate nalgebra as na;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
pub struct Equatorial {
    pub ra: f64,
    pub dec: f64,
    pub fr: f64,
    pub rot: na::Rotation3<f64>,
}

impl Equatorial {
    pub fn calculate_rotation_matrix(&mut self){
        trace!("calculate_rotation_matric_eq start");
        self.rot = na::Rotation3::from_euler_angles(self.fr, -self.dec, self.ra).inverse();
        trace!("calculate_rotation_matric_eq end");
    }

    pub fn update_equatorial_coordinates(&mut self, vec: &[f64; 3]){
        trace!("update_equatorial_coordinates start");
        self.fr  = vec[0];
        self.dec = vec[1];
        self.ra  = vec[2];

        self.calculate_rotation_matrix();
        trace!("update_equatorial_coordinates end");
    }

    pub fn new() -> Equatorial {
        let ra = 0.0;
        let dec = 0.0;
        let fr = 0.0;
        let rot = na::Rotation3::<f64>::identity();
        let eq = Equatorial{
            ra,
            dec,
            fr,
            rot,
        };
        eq
    }

}