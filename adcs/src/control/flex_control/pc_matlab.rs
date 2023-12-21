extern crate nalgebra as na;
use crate::control::gains::{Gains};

use crate::control::flex_control as fc;
use crate::control::lqr;

type BRigid = na::SMatrix<f64, 9, 4>;
type BRigid2 = na::SMatrix<f64, 3, 4>;
type Matrix1 = na::SMatrix<f64,1,1>;
type Vector9 = na::SMatrix<f64, 9, 1>;
type Matrix49 = na::SMatrix<f64, 4, 9>;
type Matrix9 = na::SMatrix<f64, 9, 9>;


#[derive(Clone)]
pub struct PassiveControl {
    num_modes: usize,
    num_states: usize,
    a_flex: na::SMatrix<f64, 104, 104>,
    b_flex: na::SMatrix<f64, 104, 5>,
    a: na::SMatrix<f64,18,18>,
    b: na::SMatrix<f64, 18,5>,
    a_cl: na::DMatrix<f64>,
    b_c: na::DMatrix<f64>,
    c_c: na::DMatrix<f64>,
    state: na::DVector<f64>,
    tau_max: f64,
    pub u: na::DVector<f64>, //the requested torque
    pub enable: bool,
    pub update: bool,
}

impl PassiveControl{
    pub fn init_pc_p() -> PassiveControl{
        let a_flex = fc::amat::init_amat();
        let b_flex = 1000.0 *fc::bmat::init_bmat();
        let a = lqr::pos::init_AMat1_pos();
        let b = lqr::pos::init_AMat2_pos();

        let num_states = 18;
        let num_inputs = 5;
        let num_modes = 6;

        let mut a_cl = na::DMatrix::<f64>::zeros(num_states, num_states);
        a_cl.slice_mut((0,0), (num_states,num_states)).copy_from(&lqr::pos::init_AMat3_pos());
        let mut b_c = na::DMatrix::<f64>::zeros(num_states, num_inputs);
        b_c.slice_mut((0,0), (num_states,num_inputs)).copy_from(&lqr::pos::init_AMat4_pos());

        let mut c_c = na::DMatrix::<f64>::zeros(num_inputs, num_states);
        c_c.slice_mut((0,0), (num_inputs,num_states)).copy_from(&lqr::pos::init_AMat5_pos());


        let state = na::DVector::<f64>::zeros(num_states);
        let u = na::DVector::<f64>::zeros(num_inputs);
        let tau_max = 4.0;
        let enable = false;
        let update = false;  


        let passive_control = PassiveControl {
            a_flex,
            b_flex,
            a,
            b,
            num_modes,
            num_states,
            a_cl,
            b_c,
            c_c,
            state,
            tau_max,
            u,
            enable,
            update
        };

        return passive_control

    }

    pub fn init_pc_n() -> PassiveControl{
        let a_flex = fc::amat::init_amat();
        let b_flex = 1000.0 *fc::bmat::init_bmat();
        let a = lqr::neg::init_AMat1_pos();
        let b = lqr::neg::init_AMat2_pos();
        let num_states = 18;
        let num_inputs = 5;
        let num_modes = 6;

        let mut a_cl = na::DMatrix::<f64>::zeros(num_states, num_states);
        a_cl.slice_mut((0,0), (num_states, num_states)).copy_from(&lqr::neg::init_AMat3_pos());
        let mut b_c = na::DMatrix::<f64>::zeros(num_states, num_inputs);
        b_c.slice_mut((0,0), (num_states,num_inputs)).copy_from(&lqr::neg::init_AMat4_pos());
        let mut c_c = na::DMatrix::<f64>::zeros(num_inputs, num_states);
        c_c.slice_mut((0,0), (num_inputs,num_states)).copy_from(&lqr::neg::init_AMat5_pos());

        let state = na::DVector::<f64>::zeros(num_states);
        let u = na::DVector::<f64>::zeros(num_inputs);
        let tau_max = 4.0;
        let enable = false;
        let update = false;  


        let passive_control = PassiveControl {
            a_flex,
            b_flex,
            a,
            b,
            num_modes,
            num_states,
            a_cl,
            b_c,
            c_c,
            state,
            tau_max,
            u,
            enable,
            update
        };

        return passive_control

    }


    pub fn reset(&mut self){
        self.state = na::DVector::<f64>::zeros(self.num_states);
    }
    

    pub fn propogate_control_state(&mut self, gyro:  &[f64], dt: f64, num_steps: u16, gyro_des: &[f64], s: na::Matrix3::<f64>){
        let s_i = s.try_inverse().unwrap();
        // let mut gyro_v = s_t*na::DVector::from_row_slice(gyro);
        // let gyro_d = s_t*s*na::DVector::from_row_slice(gyro_des);

        // Create two 3-element vectors from the 5-element gyro vector
        let gyro_124 = [gyro[0], gyro[1], gyro[3]];
        let gyro_135 = [gyro[0], gyro[2], gyro[4]];

        // Perform the multiplication for each vector
        let gyro_v_124 = s_i * na::DVector::from_row_slice(&gyro_124);
        let gyro_v_125 = s_i * na::DVector::from_row_slice(&gyro_135);

        let gyro_d_i = na::DVector::from_row_slice(gyro_des);

        // Reassemble the results into a 5-element vector
        // Average the first elements
        let avg_first_element = (gyro_v_124[0] + gyro_v_125[0]) / 2.0;
        let gyro_v = na::DVector::from_row_slice(&[
            avg_first_element,
            gyro_v_124[1],
            gyro_v_125[1],
            gyro_v_124[2], // This is actually the 4th element of the original gyro
            gyro_v_125[2], // This is the 5th element of the original gyro
        ]);

        let gyro_d = na::DVector::from_row_slice(&[
            gyro_d_i[0],
            gyro_d_i[1],
            gyro_d_i[1],
            gyro_d_i[2],
            gyro_d_i[2],
        ]);

        // println!("Testing\n");
        // let mut gyro_v = &gyro;

        // if &gyro_v.norm() > &0.0001 {
        //     gyro_v = gyro_v * 0.0;
        // }
        // println!("{}", self.eta);
        // let now = Instant::now();
        let num_states = self.num_states;
        let mut k1 = na::DVector::<f64>::zeros(num_states);
        let mut k2 = na::DVector::<f64>::zeros(num_states);
        let mut k3 = na::DVector::<f64>::zeros(num_states);
        let mut k4 = na::DVector::<f64>::zeros(num_states);
        // println!("Testing\n");
        let tau_contrib = self.b_c.clone() * (&gyro_d - &gyro_v);
        // let tau_contrib = self.b_c.clone() * (&gyro_d - (&gyro_v));

        for step in 0..num_steps{
            let state0 = self.state.clone();
            let a_cl = self.a_cl.clone();
            k1 = ((&a_cl*&state0) + &tau_contrib) * dt;
            k2 = ((&a_cl*(&state0+(&k1/2.0))) + &tau_contrib) * dt;
            k3 = ((&a_cl*(&state0+(&k2/2.0))) + &tau_contrib) * dt;
            k4 = ((&a_cl*(&state0+&k3)) + &tau_contrib) * dt;
            

            let state_dot = (k1+(2.0*k2)+(2.0*k3)+k4)/6.0;
            // println!("state dot: {}", &state_dot);
            self.state = state0 + state_dot;
        }
        // let tau_des = -&self.c_c * &self.state;
        let tau_des = &self.c_c * &self.state;
        // println!("Testing\n");
        // println!("state {} tau_des {} gyros {} gyro_des {}", &self.state, &tau_des, &gyro_v, &gyro_d);
        ////// BOUND TORQUE
        let mut tau_bound = tau_des.clone();
        // println!("trque before: {}", &tau_des);
        for (loc) in (0..5){
            let torque = tau_des[loc].clone();
            let mag = torque.abs();
            if mag > self.tau_max {
                tau_bound[loc] = torque.signum() * self.tau_max;
            }
            if loc == 0 {
                tau_bound[loc] = -tau_bound[loc];
            }

        }
        // println!("Testing\n");
        // println!("trque after: {}", &tau_bound);
        ///////
        self.u = tau_bound.clone();
        // let u_r = na::DVector::from_row_slice(&[tau_des[0], ((tau_des[1] + tau_des[2])/2.0), ((tau_des[3] + tau_des[4])/2.0)]);
        // self.u_rigid = u_r;
        // println!("tau in fc: {} \n tau rigid: {}", tau_des.clone(), self.u_rigid.clone());
        // println!("dt: {} num_steps: {}", &dt, &num_steps);
    }
}