// import crate for matrix math
extern crate nalgebra as na;
use crate::control::state::{gps, equitorial as eq, horizontal as hor, gimbal as gb};

#[derive(Debug, Clone)]
pub struct Coords {
    pub eq_k: eq::Equatorial,
    pub eq_d: eq::Equatorial,
    pub gmb_k: gb::Gimbal,
    pub gmb_d: gb::Gimbal,
}

impl Coords{
    pub fn new() -> Coords {
        let eq_k = eq::Equatorial::new();
        let eq_d = eq::Equatorial::new();
        let gmb_k = gb::Gimbal::new();
        let gmb_d = gb::Gimbal::new();

        let mut coords = Coords{
            eq_k,
            eq_d,
            gmb_k,
            gmb_d,
        };
        coords
    }
}