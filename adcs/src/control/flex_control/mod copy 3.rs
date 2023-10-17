extern crate nalgebra as na;
use crate::control::gains::{Gains};

mod amat;
pub mod bmat;
mod pc;
mod pc1;

pub use pc::*;
#[derive(Clone)]
pub struct DampingControl{
    pc_p: PassiveControl,
    pc_n: PassiveControl,
    pub err_weight: f64,
    pub u: na::DVector<f64>, //the requested torque
    pub enable: bool,
    pub update: bool,
    pub roll: f64,
    roll_torque: f64,
    pub pitch: f64,
    pitch_k: f64,
}

impl DampingControl{
    pub fn propogate(&mut self, gyro:  &[f64], dt: f64, num_steps: u16, gyro_des: &[f64]){
        self.pc_p.propogate_control_state(gyro, dt, num_steps, &gyro_des);
        self.pc_n.propogate_control_state(gyro, dt, num_steps, &gyro_des);
    }

    pub fn get_control_input(&mut self){
        let vec_pos = na::DVector::<f64>::from_row_slice(&[0.0,self.roll_torque,0.0, 0.0]);
        let vec_pitch = 1.0*na::DVector::<f64>::from_row_slice(&[0.0,0.0,self.pitch_k * (self.pitch+0.6981317007977318), 0.0]);
        let cont_pos = self.pc_p.u.clone() + vec_pos.clone();
        let cont_neg = self.pc_n.u.clone() - vec_pos.clone();
        
        let weight = (self.roll + 0.1)/0.2;
        let control_in: na::DVector<f64> = (weight * cont_pos) + ((1.0-weight)*cont_neg);
        self.u = control_in + vec_pitch.clone();
        println!("tau: {} \n pitch {} \n vec_pos {} \n weight {}",self.u, self.pitch, vec_pos+vec_pitch, weight);
    }

    pub fn reset(&mut self){
        self.pc_p.reset();
        self.pc_n.reset();
    }

    pub fn init_pc(gains: &Gains) -> DampingControl {
        let pc_p = PassiveControl::init_pc_p();
        let pc_n = PassiveControl::init_pc_n();
        let u = na::DVector::<f64>::zeros(5);
        let err_weight = 0.0;
        let enable = true;
        let update = false;
        let roll = 0.0;
        let roll_torque = 6.260767;
        let pitch = 0.0;
        let pitch_k = 62.60767;

        let damping_control = DampingControl {
            pc_p,
            pc_n,
            err_weight,
            u,
            enable,
            update,
            roll,
            roll_torque,
            pitch,
            pitch_k
        };
        return damping_control;

    }
}