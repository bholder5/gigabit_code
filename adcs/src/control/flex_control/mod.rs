extern crate nalgebra as na;
use crate::control::gains::{Gains};

mod amat;
pub mod bmat;
// mod pc;
mod pc_matlab;

pub use pc_matlab::*;
#[derive(Clone)]
pub struct DampingControl{
    pc_p: PassiveControl,
    pc_n: PassiveControl,
    pub err_weight: f64,
    pub u: na::DVector<f64>, //the requested torque
    pub piv_speed: f64,
    pub i_yaw: f64,
    pub enable: bool,
    pub update: bool,
    pub roll: f64,
    roll_torque: f64,
    roll_torque_max: f64,
    pub pitch: f64,
    pitch_k: f64,
    pub ff_r: f64,
    pub ff_p: f64,
    pub fine_point: bool,
    pub ff_flag: bool,
    pub pitch_nom: f64,
    pub pitch_nom_d: f64,
    pub move_pitch_nom: bool,
}

impl DampingControl{
    pub fn propogate(&mut self, gyro:  &[f64], dt: f64, num_steps: u16, gyro_des: &[f64]){
        let sp = na::Matrix3::<f64>::from_row_slice(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.1_f64.sin(), 0.0, 1.0]);
        let sn = na::Matrix3::<f64>::from_row_slice(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -0.1_f64.sin(), 0.0, 1.0]);

        self.pc_p.propogate_control_state(gyro, dt, num_steps, &gyro_des, sp);
        self.pc_n.propogate_control_state(gyro, dt, num_steps, &gyro_des, sn);
    }

    pub fn get_control_input(&mut self){
        if self.roll_torque < self.roll_torque_max{
            self.roll_torque = self.roll_torque + 0.005;
        } else {
            self.ff_flag = true;
        }
        
        let ff_p_d = 0.99* self.pitch * self.pitch_k;

        // if ff_p_d.abs() > 4.0 {
        //     // self.move_pitch_nom = true;
        //     self.pitch_nom_d = self.pitch_nom + ff_p_d.signum() * 0.01;
        // } else if self.ff_p.abs() < ff_p_d.abs(){
        //     self.ff_p = self.ff_p + ff_p_d.signum() * 0.005;
        //     self.move_pitch_nom = false;
        // }

        if ff_p_d.abs() > self.roll_torque.abs(){
            self.ff_p = ff_p_d.signum() * self.roll_torque.abs();
            // println!("ff_p_d: {} roll torque {} ff_p {} pitch: {}", &ff_p_d, &self.roll_torque, &self.ff_p, &self.pitch);
        } else {
            self.ff_p = ff_p_d;
            // println!("ELSE ff_p_d: {} roll torque {} ff_p {}", &ff_p_d, &self.roll_torque, &self.ff_p);
        }

        // println!("pitch ffd {}", self.ff_p);

        // let vec_pitch = 1.0*na::DVector::<f64>::from_row_slice(&[0.0,0.0, 0.0, self.pitch_k * (self.pitch+0.6981317007977318),self.pitch_k * (self.pitch+0.6981317007977318)]);
        
        let sp = na::DMatrix::<f64>::from_row_slice(3,3,&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.1_f64.sin(), 0.0, 1.0]);
        let sn = na::DMatrix::<f64>::from_row_slice(3,3,&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -0.1_f64.sin(), 0.0, 1.0]);
        
        let sts_p = sp.transpose() * sp;
        let stsi_p = sts_p.try_inverse().unwrap();

        let sts_n = sn.transpose() * sn;
        let stsi_n = sts_n.try_inverse().unwrap();
        
        // let cont_pos = (stsi_p*self.pc_p.u.clone());
        // let cont_neg = (stsi_n * self.pc_n.u.clone());

        let cont_pos = 1.0*self.pc_p.u.clone();
        let cont_neg = 1.0*self.pc_n.u.clone();
        

        let pos_1 = &stsi_p*na::DVector::<f64>::from_row_slice(&[cont_pos[0], cont_pos[1], cont_pos[3]]);
        let pos_2 = &stsi_p*na::DVector::<f64>::from_row_slice(&[cont_pos[0], cont_pos[2], cont_pos[4]]);
        let posi = na::DVector::<f64>::from_row_slice(&[pos_1[0], pos_1[1], pos_2[1], pos_1[2], pos_2[2]]);

        let neg_1 = &stsi_n*na::DVector::<f64>::from_row_slice(&[cont_neg[0], cont_neg[1], cont_neg[3]]);
        let neg_2 = &stsi_n*na::DVector::<f64>::from_row_slice(&[cont_neg[0], cont_neg[2], cont_neg[4]]);
        let negi = na::DVector::<f64>::from_row_slice(&[neg_1[0], neg_1[1], neg_2[1], neg_1[2], neg_2[2]]);
        
        
        let weight = (self.roll + 0.1)/0.2;
        let u: na::DVector<f64> = (weight * posi) + ((1.0-weight)*negi);
        self.u = u.clone();
        self.ff_r = 0.99* ((weight * self.roll_torque) + ((1.0-weight)*-1.0*self.roll_torque));
        println!("ffd roll {} ff_p {}", self.ff_r, self.ff_p);

        // println!("tau: {} \n pitch {} \n vec_pos {} \n weight {}",self.u, self.pitch, vec_pos+vec_pitch, weight);
    }

    pub fn reset(&mut self){
        self.pc_p.reset();
        self.pc_n.reset();
    }

    pub fn init_pc(gains: &Gains) -> DampingControl {
        let pc_p = PassiveControl::init_pc_p();
        let pc_n = PassiveControl::init_pc_n();
        let u = na::DVector::<f64>::zeros(5);
        let piv_speed = 0.0;
        let i_yaw = 979.0;
        let err_weight = 0.0;
        let enable = true;
        let update = false;
        let roll = 0.0;
        let roll_torque = 0.060767;
        let roll_torque_max = 6.260767;
        let pitch = 0.0;
        let pitch_k = 62.60767;
        let ff_r = 0.0;
        let ff_p = 0.0;
        let fine_point = false;
        let ff_flag = false;
        let pitch_nom = 0.0;
        let pitch_nom_d = 0.0;
        let move_pitch_nom = true;

        let damping_control = DampingControl {
            pc_p,
            pc_n,
            err_weight,
            u,
            piv_speed,
            i_yaw,
            enable,
            update,
            roll,
            roll_torque,
            roll_torque_max,
            pitch,
            pitch_k,
            ff_r,
            ff_p,
            fine_point,
            ff_flag,
            pitch_nom,
            pitch_nom_d,
            move_pitch_nom,
        };
        return damping_control;

    }
}