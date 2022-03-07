// Import modules (function and types)
mod types;
mod functions;
mod constants;

use types::*;
use functions::*;
use constants::*;

// import for writing csv
extern crate csv;
extern crate serde;
#[macro_use]
extern crate serde_derive;
use std::process;
use nalgebra::{DMatrix};

/// linear systems crates etc for automatica
extern crate automatica;
// use num_complex::Complex;
#[allow(unused_imports)]
use automatica::{signals::continuous, Seconds, Ss, TfMatrix};

#[allow(unused_assignments)]
#[allow(unused_variables)]
#[allow(dead_code)]
// for linear systems automatica toolbox examples
#[allow(clippy::many_single_char_names)]
#[allow(clippy::non_ascii_literal)]
fn main() {
    // initializing the model of the system
    let mut _model = SSmodel {
        _mass_mat: build_mass(MOI_IR, MOI_I1, MOI_I2),
        _stiff_mat: build_stiffness(K1, K2),
        _a_mat: DMatrix::<f64>::zeros(6,6),
        _b_mat: DMatrix::<f64>::zeros(6,1),
        _c_mat: DMatrix::<f64>::from_vec(1,6, vec![1.0,0.0,0.0,0.0,0.0,0.0]),
        _d_mat: DMatrix::<f64>::from_vec(1,1, vec![0.0]),
    };

    //assign identity to block 12 ode of theta-theta_dot
    _model._a_mat = mat_block_change(DMatrix::<f64>::identity(3,3), _model._a_mat, 2);
    // 

    //calculate block 21 using inverse.
    match _model._mass_mat.try_inverse() {
        Some(inv) => {
            let b_inv = inv.clone();
            _model._a_mat = mat_block_change(-inv * _model._stiff_mat, _model._a_mat, 3);
            // println!("{:.3}", _model._a_mat);
            //
            // build b matrix (botto involves inverse of mass)
            let mut _b_tdd = DMatrix::<f64>::from_vec(3,1, vec![1.0, 0.0, 0.0]);
            _b_tdd = b_inv * _b_tdd;
            _model._b_mat = mat_block_change(_b_tdd, _model._b_mat, 4);
            // println!("{:.3}", _model._b_mat);

        }
        None => {
            println!("mass matrix is not invertible!");
            process::exit(1);
        }
    }
// 
// 
//
    let a_mat = _model._a_mat.transpose();

// Build the statespace model of the yaw 
    let sys = Ss::new_from_slice(6, 1, 1, 
        a_mat.as_slice(), _model._b_mat.as_slice(),
        _model._c_mat.as_slice(), _model._d_mat.as_slice());

    // Initialize sim loop, looping sim at 10 times the frequency of the control.
    

    let ctrl_freq: f64 = 100.; // Hz
    let ctrl_step: f64 = 1.0/ctrl_freq;

    let loop_steps_per_ctrl: usize = 10;

    let loop_freq: f64 = ctrl_freq * (loop_steps_per_ctrl as f64);
    let loop_step: f64 = 1.0 / loop_freq;

    let mut t_cur: f64 = 0.0;
    let t_end: f64 = 200.0;
    let mut x_cur: Vec<f64> = vec![0., 0., 0., 0., 0., 0.];
    let mut x_prv: Vec<f64> = vec![0., 0., 0., 0., 0., 0.];

    // Controller constants and variables. 
    let xref:f64 = 0.1;
    let kp: f64 = 2.0;
    let ki: f64 = 0.0;
    let kd: f64 = 0.5;
    let mut e_cum: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut u: f64 = 0.0;
    let td_max: f64 = 3.;
    let mut td_1: f64 = 0.0;

    // diff eq implem of PDF
    let cont_a: f64 = 21450.0;
    let cont_b: f64 = -21449.9;
    let cont_c: f64 = 0.0006738;

    let mut e: f64 = 0.0;
    let mut ekm1: f64 = 0.0;
    let mut ukm1: f64 = 0.0;


    while t_cur <= t_end{
        // determine control input here
        //
        // Try PID control assuming positional feedback. 
        //
        e = xref - x_cur[0].clone();
        

        u = cont_a * e + cont_b * ekm1 + cont_c * ukm1;
        
        ukm1 = u.clone();
        ekm1 = e.clone();
        // inital PID control
        //
        //calulate error terms;
        // e = xref - x_cur[0];
        // e_cum += e;
        // if e_cum > 10.0{
        //     e_cum = 10.0;
        // } else if e_cum < -10.0{
        //     e_cum = -10.0;
        // }

        // if x_cur[3] < td_max{
        //     if x_cur[3] > -td_max{
        //         td_1 = x_cur[3];
        //     } else {
        //         td_1 = -td_max;
        //     }
        // } else {
        //     td_1 = td_max;
        // }
        // //
        // //apply control
        // u = (kp * e) + (ki * e_cum) + (kd * td_1);
        

        //
        //
        
        let input = continuous::step(u, 1);
        //
        // Simulate one ctrl loop frequencys worth of system
        let rk4 = sys.rk4(&input, &x_cur, Seconds(loop_step), loop_steps_per_ctrl);
        // let rk5 = rk4.clone();
        // println!("{:?}", rk5);
        // write data to csv file (including calculating proper time)
        // if env::args_os().len() > 1{
            let mut _count = 0;
            for i in rk4 {
                if _count < loop_steps_per_ctrl {
                    if let Err(err) = write_func_automatica(&i, u, t_cur) 
                    {
                        println!("{}", err);
                        process::exit(1);
                    }
                    // let last = i;
                    // x_cur[0] = last.state()[0];
                    // x_cur[1] = last.state()[1];
                    // x_cur[2] = last.state()[2];
                    // x_cur[3] = last.state()[3];
                    // x_cur[4] = last.state()[4];
                    // x_cur[5] = last.state()[5];
                    // let t= last.time().to_string().parse::<f64>().unwrap() + t_cur;
                    // println!("time: {:.3}, State: {:.3?}, u: {:.3}", t, x_cur, u);
                } else {
                    // Reset initial condition (x_cur) to last state, don't print to csv as
                    // as it will get printed on next round.
                    x_prv = x_cur.clone();

                    let last = i;
                    x_cur[0] = last.state()[0];
                    x_cur[1] = last.state()[1];
                    x_cur[2] = last.state()[2];
                    x_cur[3] = last.state()[3];
                    x_cur[4] = last.state()[4];
                    x_cur[5] = last.state()[5];
                    t_cur    = last.time().to_string().parse::<f64>().unwrap() + t_cur;
                    println!("time: {:.3}, State: {:.3?}, u: {:.3}, err: {:.3}, int: {:.3}", t_cur, x_cur, u, e, e_cum);
                }
                _count += 1;
            }
        // }
        // // reset initial condition to last state
        // t_cur = ((t_cur + loop_step) * 10.0).round() /10.0;
        // println!("time: {}", t_cur);
    }




    // Change to 'true' to print the result
    // if env::args_os().len() > 1{ 
        // for i in rk4 {
            // if let Err(err) = write_func_automatica(&i) {
            //     println!("{}", err);
            //     process::exit(1);
            // }
            // println!("{:?}", &i);
        // }
    // } else {
    //     println!("No arugment provided, csv not written");
    // }

}