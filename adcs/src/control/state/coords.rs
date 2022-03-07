// import crate for matrix math
extern crate nalgebra as na;
use crate::control::state::{gps, equitorial as eq, horizontal as hor, gimbal as gb};

#[derive(Debug, Clone)]
pub struct Coords {
    gps: gps::Gps,
    pub eq_k: eq::Equatorial,
    pub eq_d: eq::Equatorial,
    hor: hor::Horizontal,
    pub gmb_k: gb::Gimbal,
    pub gmb_d: gb::Gimbal,
    ceh: na::Rotation3<f64>,
    b2h_offset: na::Rotation3<f64>,
}

impl Coords{
    pub fn new() -> Coords {
        let gps = gps::Gps::new();
        let eq_k = eq::Equatorial::new();
        let eq_d = eq::Equatorial::new();
        let hor = hor::Horizontal::new();
        let gmb_k = gb::Gimbal::new();
        let gmb_d = gb::Gimbal::new();
        let ceh = na::Rotation3::<f64>::identity();
        let b2h_offset = na::Rotation3::<f64>::identity();

        let mut coords = Coords{
            gps, 
            eq_k,
            eq_d,
            hor,
            gmb_k,
            gmb_d,
            ceh,
            b2h_offset,
        };

        coords.update_hor_to_eq_conversion();
        coords.update_horizontal_coordinates();
        coords
    }

    pub fn update_hor_to_eq_conversion(&mut self){
        let last = self.gps.gast + self.gps.lon;
        let arg_z = (last + 180.0).to_radians();

        let arg_y = (self.gps.lat - 90.0).to_radians();
        
        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), arg_y);
        let rotz = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::z_axis(), arg_z);

        println!("arg y {:}, argz {:}", arg_y, arg_z);
        println!("rot y {:}, rotz {:}, ceh {:}", roty, rotz, rotz*roty);
        self.ceh = rotz*roty;
    }

    pub fn update_horizontal_coordinates(&mut self){
        self.hor.rot = self.ceh * self.eq_k.rot;
        self.hor.extract_az_el_ir_from_rotmat();

    }
}