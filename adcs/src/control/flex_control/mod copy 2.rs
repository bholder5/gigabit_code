extern crate nalgebra as na;
use crate::control::gains::{Gains};

mod amat;
mod bmat;
mod pc;

pub use pc::*;

#[derive(Clone)]
pub struct DampingControl{
    pc_slew: PassiveControl,
    pc_fine: PassiveControl,
    pub err_weight: f64,
    pub u: na::DVector<f64>, //the requested torque
    pub enable: bool,
    pub update: bool,
    pub r_sc_slew: f64,
    pub qr_sc_slew: f64,
    pub ql_sc_slew: f64,
    pub r_sc_fine: f64,
    pub qr_sc_fine: f64,
    pub ql_sc_fine: f64,
}

impl DampingControl{
    pub fn propogate(&mut self, gyro:  &[f64], dt: f64, num_steps: u16, gyro_des: &[f64]){
        self.pc_slew.propogate_control_state(gyro, dt, num_steps, &gyro_des);
        self.pc_fine.propogate_control_state(gyro, dt, num_steps, &gyro_des);

        let u = (1.0-self.err_weight) * self.pc_fine.u.clone() + (self.err_weight) * self.pc_slew.u.clone();
        self.u = u.clone()
    }

    pub fn reset(&mut self){
        self.pc_slew.reset();
        self.pc_fine.reset();
    }

    pub fn init_pc(gains: &Gains) -> DampingControl {
        let r_sc_slew = 0.5;
        let qr_sc_slew = 0.01;
        let ql_sc_slew = 100.0;
        let pc_slew = PassiveControl::init_pc(gains, r_sc_slew, qr_sc_slew, ql_sc_slew);
        
        let r_sc_fine = 15.0;
        let qr_sc_fine = 2.0;
        let ql_sc_fine = 100.0;
        let pc_fine = PassiveControl::init_pc(gains, r_sc_fine, qr_sc_fine, ql_sc_fine);
        let u = na::DVector::<f64>::zeros(5);
        let err_weight = 0.0;
        let enable = true;
        let update = false;

        let damping_control = DampingControl {
            pc_slew,
            pc_fine,
            err_weight,
            u,
            enable,
            update,
            r_sc_slew,
            qr_sc_slew,
            ql_sc_slew,
            r_sc_fine,
            qr_sc_fine,
            ql_sc_fine,

        };
        return damping_control;

    }
}