// import crate for matrix math
extern crate nalgebra as na;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
pub struct Horizontal {
    az: f64,
    el: f64,
    ir: f64,
    pub rot: na::Rotation3<f64>,
}

impl Horizontal {
    pub fn calculate_rotation_matrix(&mut self){
        self.rot = na::Rotation3::from_euler_angles(self.ir, -self.el, -self.az).inverse();
    }

    pub fn update_horizontal_coordinates(&mut self, vec: &[f64; 3]){
        self.ir = vec[0];
        self.el = vec[1];
        self.az = vec[2];

        self.calculate_rotation_matrix();
    }

    pub fn extract_az_el_ir_from_rotmat(&mut self){
        debug!("cbh{:}", self.rot);

        let az = (-self.rot[(0,1)]).atan2(self.rot[(0,0)]);
        let el = self.rot[(0,2)].asin();
        let ir = self.rot[(1,2)].atan2(self.rot[(2,2)]);
        debug!("az {}, el {}, ir {}", az.to_degrees(),el.to_degrees(),ir.to_degrees());
    }

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