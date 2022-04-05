//! Outlining the Coords (Coordinates) structure and implementing the associated functions for it for use in the state sub-module of the control sub-module of the ADCS crate
// import crate for matrix math
extern crate nalgebra as na;
use crate::control::state::{equitorial as eq, gimbal as gb};

#[derive(Debug, Clone)]
/// An umbrella structure to hold the the various current and desired orientations in the gimbal and equatorial coordinate frames
pub struct Coords {
    /// The current orientation in equatorial coordinates
    pub eq_k: eq::Equatorial,
    /// The desired orientation in equatorial coordinates
    pub eq_d: eq::Equatorial,
    /// The current state in gimbal coordinates
    pub gmb_k: gb::Gimbal,
    /// The desired state in gimbal coordinates
    pub gmb_d: gb::Gimbal,
}

impl Coords {
    pub fn new() -> Coords {
        let eq_k = eq::Equatorial::new();
        let eq_d = eq::Equatorial::new();
        let gmb_k = gb::Gimbal::new();
        let gmb_d = gb::Gimbal::new();

        let coords = Coords {
            eq_k,
            eq_d,
            gmb_k,
            gmb_d,
        };
        coords
    }
}
