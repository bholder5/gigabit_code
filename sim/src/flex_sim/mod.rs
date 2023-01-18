mod amat;
mod bmat;
mod cmat;
mod kmat;

pub use amat::*;
pub use bmat::*;
pub use cmat::*;
pub use kmat::*;

extern crate nalgebra as na;
use std::time::{Duration, Instant};

pub type Eta = na::SMatrix<f64, 104, 1>;

pub struct Flex_model {
    // //
    a_mat: AMat,
    b_mat: BMat,
    k_mat: KMat,
    g1_mat: CMat,
    g2_mat: CMat,
    g1_pos_mat: CMat,
    c_pos_mat: CMat_gimb,
    pub eta: Eta,
    pub g1_out: na::Vector3<f64>,
    pub g2_out: na::Vector3<f64>,
    pub c_out: na::Vector5<f64>,
    pub g1_pos_out: na::Vector3<f64>,
    pub c_pos_out: na::Vector5<f64>,
    pub flex_enable: bool,
}
impl Flex_model{
    pub fn propogate_flex(&mut self, tau:  &[f64; 3], dt: f64, num_steps: u16){
        let tau = na::Vector5::from_row_slice(&[tau[0], tau[1]/2.0, tau[1]/2.0, tau[2]/2.0, tau[2]/2.0]);
        // println!("{}", self.eta);
        let now = Instant::now();
        for step in 0..num_steps{
            let k1 = ((self.a_mat*self.eta) + (self.b_mat * tau))*dt;            
            let k2 = ((self.a_mat*(self.eta + (k1/2.0))) + (self.b_mat * tau))*dt;
            let k3 = ((self.a_mat*(self.eta + (k2/2.0))) + (self.b_mat * tau))*dt;
            let k4 = ((self.a_mat*(self.eta + k3))     + (self.b_mat * tau))*dt;

            let eta_dot = (k1+(2.0*k2)+(2.0*k3)+k4)/6.0;
            self.eta = self.eta + eta_dot;

            self.update_gyro_outputs();
        }
        println!("elapsed: {}", now.elapsed().as_micros());
    }
    
    fn update_gyro_outputs(&mut self){
        self.g1_out = self.g1_mat * self.eta;
        self.g2_out = self.g2_mat * self.eta;
        self.c_out  = self.b_mat.transpose() * self.eta;
        self.g1_pos_out = self.g1_pos_mat * self.eta;
        self.c_pos_out = self.c_pos_mat * self.eta;
    }
} 

pub fn init_flex_model() -> Flex_model{
    let a_mat = init_amat();
    let b_mat = 1000.0 * init_bmat();
    let k_mat = init_kmat();
    let g1_mat = 1000.0 * init_g1mat();
    let g2_mat = 1000.0 * init_g2mat();
    let g1_pos_mat = 1000.0 * init_g1_pos_mat();
    let c_pos_mat = 1000.0 * init_c_pos_mat();
    let eta = Eta::zeros();
    let g1_out = na::Vector3::<f64>::zeros();
    let g2_out = na::Vector3::<f64>::zeros();
    let c_out = na::Vector5::<f64>::zeros();
    let g1_pos_out = na::Vector3::<f64>::zeros();
    let c_pos_out = na::Vector5::<f64>::zeros();
    let flex_enable: bool = true;

    let flex_model = Flex_model{
        a_mat,
        b_mat,
        k_mat,
        g1_mat,
        g2_mat,
        g1_pos_mat,
        c_pos_mat,
        eta,
        g1_out,
        g2_out,
        c_out,
        g1_pos_out,
        c_pos_out,
        flex_enable,
    };
    return flex_model
}