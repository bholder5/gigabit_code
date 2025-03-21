mod amat;
mod cmat;
mod kmat;

pub use amat::*;
// pub use bmat::*;
pub use cmat::*;
pub use kmat::*;

use adcs::control::flex_control::bmat::{init_bmat,BMat};

extern crate nalgebra as na;
use std::time::{Duration, Instant};
use adcs::control::flex_control::PassiveControl;

pub type Eta = na::SMatrix<f64, 104, 1>;
pub type HEta = na::SMatrix<f64, 52, 1>;
pub type Matrix1 = na::SMatrix<f64, 1, 1>;

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
    pub fn _propogate_flex(&mut self, tau:  &[f64; 3], dt: f64, num_steps: u16, fc: &PassiveControl){
        let tau = na::Vector5::from_row_slice(&[tau[0] + fc.u[0], (tau[1]/2.0) + fc.u[1], (tau[1]/2.0) + fc.u[2], (tau[2]/2.0) + fc.u[3], (tau[2]/2.0) + fc.u[4]]);
        // println!("{}", self.eta);
        // let now = Instant::now();
        let mut k1 = Eta::zeros();
        let mut k2 = Eta::zeros();
        let mut k3 = Eta::zeros();
        let mut k4 = Eta::zeros();

        let tau_contrib = self.b_mat * tau;

        for step in 0..num_steps{
            let eta0 = self.eta.clone();
            k1 = (self.calc_eta_dot(&eta0) + tau_contrib) * dt;
            // let k1 = ((self.a_mat*self.eta) + tau_contrib)*dt; 
            // k2 = ((self.a_mat*(self.eta + (k1/2.0))) + tau_contrib)*dt;
            k2 = (self.calc_eta_dot(&(&eta0 + (k1/2.0))) + tau_contrib)*dt;
            // println!("Test of new calc {}", k2-k2a);
            // k3 = ((self.a_mat*(self.eta + (k2/2.0))) + tau_contrib)*dt;
            k3 = (self.calc_eta_dot(&(&eta0 + (k2/2.0))) + tau_contrib)*dt;
            // k4 = ((self.a_mat*(self.eta + k3))     + tau_contrib)*dt;
            k4 = (self.calc_eta_dot(&(&eta0 + (k3))) + tau_contrib)*dt;

            let eta_dot = (k1+(2.0*k2)+(2.0*k3)+k4)/6.0;
            self.eta = self.eta + eta_dot;
        }
        self.update_gyro_outputs();
        // println!("elapsed: {}", now.elapsed().as_micros());
    }

    fn calc_eta_dot(&mut self, eta0: &Eta) -> Eta{
        let eta = HEta::from_row_slice(&eta0.as_slice()[0..52]);
        let d_eta = HEta::from_row_slice(&eta0.as_slice()[52..104]);

        let out2 = (-self.k_mat * eta);
        let out = Eta::from_row_slice(&[&d_eta.as_slice()[0..52], &out2.as_slice()[0..52]].concat());
        return out
    }
    
    fn update_gyro_outputs(&mut self){
        self.g1_out = self.g1_mat * self.eta;
        self.g2_out = self.g2_mat * self.eta;
        self.c_out  = self.b_mat.transpose() * self.eta;
        self.g1_pos_out = self.g1_pos_mat * self.eta;
        self.c_pos_out = self.c_pos_mat * self.eta;
    }

    pub fn update_eta(&mut self, eta1: &[f64]){
        self.eta = Eta::from_row_slice(eta1.as_ref());
        // println!("eta {}", &self.eta[0]);
        self.update_gyro_outputs();

    }
    /// Reorg the cmatrices to reflect the position then rate of a single mode following each other in the 
    /// state vector rather than all positions followed by all rates
    pub fn reorg_cmats(&mut self){

        // g1_mat
        let tmp = self.g1_mat.slice((0,0), (3,52)).clone();
        let tmp2 = self.g1_mat.slice((0,52), (3,52)).clone();

        let mut g1_mat = self.g1_mat.clone();

        let mut num_col = 0;

        // println!("{:.2}", &self.g1_mat.transpose());

        for col in tmp2.column_iter(){
            g1_mat.slice_mut((0,num_col*2), (3,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        let mut num_col = 1;
        for col in tmp.column_iter(){
            g1_mat.slice_mut((0,(num_col*2)-1), (3,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        // println!("{:.2}", &g1_mat.transpose());
        self.g1_mat = g1_mat.clone();

        // g2_mat
        let tmp = self.g2_mat.slice((0,0), (3,52)).clone();
        let tmp2 = self.g2_mat.slice((0,52), (3,52)).clone();

        let mut g2_mat = self.g2_mat.clone();

        let mut num_col = 0;

        // println!("{:.5}", &self.g2_mat.transpose());

        for col in tmp2.column_iter(){
            g2_mat.slice_mut((0,num_col*2), (3,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        let mut num_col = 1;
        for col in tmp.column_iter(){
            g2_mat.slice_mut((0,(num_col*2)-1), (3,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        // println!("{:.5}", &g2_mat.transpose());
        self.g2_mat = g2_mat.clone();

        // g1_pos_mat
        let tmp = self.g1_pos_mat.slice((0,0), (3,52)).clone();
        let tmp2 = self.g1_pos_mat.slice((0,52), (3,52)).clone();

        let mut g1_pos = self.g1_pos_mat.clone();

        let mut num_col = 0;

        // println!("{:.2}", &self.g1_pos_mat.transpose());

        for col in tmp2.column_iter(){
            g1_pos.slice_mut((0,num_col*2), (3,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        let mut num_col = 1;
        for col in tmp.column_iter(){
            g1_pos.slice_mut((0,(num_col*2)-1), (3,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        // println!("{:.2}", &g1_pos.transpose());
        self.g1_pos_mat = g1_pos.clone();

        // c_pos_mat
        let tmp = self.c_pos_mat.slice((0,0), (5,52)).clone();
        let tmp2 = self.c_pos_mat.slice((0,52), (5,52)).clone();

        let mut c_pos = self.c_pos_mat.clone();

        let mut num_col = 0;

        // println!("{:.2}", &self.c_pos_mat.transpose());

        for col in tmp.column_iter(){
            c_pos.slice_mut((0,num_col*2), (5,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        let mut num_col = 1;
        for col in tmp2.column_iter(){
            c_pos.slice_mut((0,(num_col*2)-1), (5,1)).copy_from(&col);
            num_col = num_col + 1;
        }
        // println!("{:.2}", &c_pos.transpose());
        self.c_pos_mat = c_pos.clone();
        // println!("cposmat {:.2}", &self.c_pos_mat.transpose());

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
    let mut eta = Eta::zeros();

    let val = 0.0;
    // eta.slice_mut((1,0), (1,1)).copy_from(&(Matrix1::from_row_slice(&[val])));
    // eta.slice_mut((3,0), (1,1)).copy_from(&(Matrix1::from_row_slice(&[val])));
    // eta.slice_mut((5,0), (1,1)).copy_from(&(Matrix1::from_row_slice(&[val])));
    // eta.slice_mut((7,0), (1,1)).copy_from(&(Matrix1::from_row_slice(&[val])));
    for i in (1..=104).step_by(2) {
        eta.slice_mut((i, 0), (1, 1)).copy_from(&(Matrix1::from_row_slice(&[val])));
    }
    

    let g1_out = na::Vector3::<f64>::zeros();
    let g2_out = na::Vector3::<f64>::zeros();
    let c_out = na::Vector5::<f64>::zeros();
    let g1_pos_out = na::Vector3::<f64>::zeros();
    let c_pos_out = na::Vector5::<f64>::zeros();
    let flex_enable: bool = false;

    // println!("bmat: {}, g1_mat {}, g2mat {}, g1pos_mat {}, c_pos_mat {}", &b_mat, &g1_mat, &g2_mat, &g1_pos_mat, &c_pos_mat);

    let mut flex_model = Flex_model{
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
    flex_model.reorg_cmats();
    // println!("bmat: {}, g1_mat {}, g2mat {}, g1pos_mat {}, c_pos_mat {}", &flex_model.b_mat, &flex_model.g1_mat.transpose(), &flex_model.g2_mat.transpose(), &flex_model.g1_pos_mat.transpose(), &flex_model.c_pos_mat.transpose());

    return flex_model;
}
