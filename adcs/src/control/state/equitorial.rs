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
        // println!("{}",self.rot);
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

    pub fn extract_ra_dec_fr_from_rotmat(&mut self){
        debug!("c {:}", self.rot);

        let ra = (self.rot[(0,1)]).atan2(self.rot[(0,0)]);
        let dec = self.rot[(0,2)].asin();
        let fr = self.rot[(1,2)].atan2(self.rot[(2,2)]);
        // println!("ra {}, dec {}, fr {}", ra,dec,fr);

        self.ra = ra.clone();
        self.dec = dec.clone();
        self.fr = fr.clone();
    }

    pub fn new() -> Equatorial {
        let ra = 0.0;
        let dec = 0.0;
        let fr = 0.0;
        let rot = na::Rotation3::<f64>::identity();
        let mut eq = Equatorial{
            ra,
            dec,
            fr,
            rot,
        };
        eq.calculate_rotation_matrix();
        eq.extract_ra_dec_fr_from_rotmat();

        eq
    }

}