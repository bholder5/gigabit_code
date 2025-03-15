#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

// include!("bindings.rs");

#[allow(dead_code)]
// #[allow(unused_attributes)]
mod bindings;
mod initialization;
mod json_handling;
mod logging;
mod sim_csvf64;
mod measurements;
mod disturbances;
mod flex_sim;

use adcs::{control::*, estimation::*, verification::*, dist_est};
use initialization::*;
use json_handling as js;
use logging::*;
use sim_csvf64 as sc;
use measurements as ms;
use disturbances as dist;
use flex_sim as f_sim;

use bindings::bit_one_step;

pub use log::LevelFilter;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
pub use log4rs::append::file::FileAppender;
pub use log4rs::config::{Appender, Config, Root};
pub use log4rs::encode::pattern::PatternEncoder;
use std::time::{Duration, Instant};

extern crate csv;
extern crate serde;
extern crate time;
extern crate nalgebra;

use nalgebra as na;
// #[macro_use]
extern crate serde_derive;

fn main() {

    init_log();
    // the verification function for all tihngs in the control module.
    // test_control();
    // test_estimation();
    // run initialization and have all relevant
    // bit parameters in the struct
    let mut bp = init_bit();
    let mut flex = f_sim::init_flex_model();
    bp.get_orientation_vec(&flex);
    

    let mut est = Estimator::new();
    let mut ctrl = Ctrl::new();
    let mut fc = flex_control::DampingControl::init_pc(&ctrl.fine_gains);
    let mut servo = servo::ServoControl::init();
    let mut sim_state = state::State::new();
    let mut meas = ms::Meas::new();
    let mut late_meas = ms::Meas::new();
    let mut wind = dist::Wind_torque::new();

    // define Simulation parameters
    let mut step: u64 = 0; //integer for tracking rate of control
    let mut tau_applied: [f64; 9] = [0., 0., 0., 0., 0., 0., 0., 0., 0.];
    let mut t = 0.0;

    //********************************************************************
    //********************************************************************
    //                   SIMULATION BEGINS
    //********************************************************************

    let mut y_result: [f64; 21] = [0.; 21];
    let mut flex_result: [f64; 104] = [0.; 104];

    js::read_gains(&mut ctrl, &mut bp, &mut est, &mut fc, &mut flex); // read in gains from json file (for tuning)
                               // ctrl.update_gain_matrices();
                               // ctrl.update_desired_orientation_matrix();
    // sc::read_last_state(t, &mut bp, &mut est, &mut ctrl, &mut meas, &mut sim_state);
                               // bp.get_orientation_vec();

    tau_applied[6] = ctrl.rw.tau_applied.clone(); //yaw
    tau_applied[7] = ctrl.fmot_roll.tau_applied.clone(); //roll
    tau_applied[8] = ctrl.fmot_pitch.tau_applied.clone(); //pitch

    
    sc::push_record(&t, &bp, &est, &ctrl, &meas, &sim_state, &flex, &fc, &servo).unwrap();

    let mut first_LIS = true;

    // Initialize with an achievable target
    ctrl.state.gmb_d.roll = 0.1;
    ctrl.state.gmb_d.pitch = -0.64;
    ctrl.state.gmb_d.yaw = -0.1;
    ctrl.state.gmb_d.calculate_rotation_matrix();
    ctrl.state.update_desired_eq_from_gmb();
    let mut servo_multiplier = 0.0;
    let mut tau_req_yaw = 0.0;
    let mut tau_req_roll = 0.0;
    let mut tau_req_pitch = 0.0;

    let mut total_duration = Duration::new(0,0);

    trace!("START");
    let now1 = Instant::now();

    for _step in 0..900000 as usize {
        
        ///////// beginning of the simulation loop
        /////////////////////////////////////////
        // fc.u = fc.u * 0.0;
        // println!("Tau applied: {} {} {}", tau_applied[6], tau_applied[7], tau_applied[8]);
        // println!("yaw:{} \nrol: {} \n pit: {}\n rw: {}\n bow: {}\nstn: {}\nprt: {}\n sb: {}\n t: {}\n", tau_applied[6], tau_applied[7], tau_applied[8], fc.u[0], fc.u[1], fc.u[2], fc.u[3], fc.u[4], &t);
        // tau_applied = [0., 0., 0., 0., 0., 0., 0., 0., 0.];

        // if _step == 30000{
        //     ctrl.state.gmb_d.roll = 0.01;
        //     ctrl.state.gmb_d.pitch = -0.69;
        //     ctrl.state.gmb_d.yaw = -0.05;
        //     ctrl.state.gmb_d.calculate_rotation_matrix();
        //     ctrl.state.update_desired_eq_from_gmb();
        // } //else if _step == 60000 {
        //     ctrl.state.gmb_d.roll = 0.0;
        //     ctrl.state.gmb_d.pitch = -0.69;
        //     ctrl.state.gmb_d.yaw = 0.01;
        //     ctrl.state.gmb_d.calculate_rotation_matrix();
        //     ctrl.state.update_desired_eq_from_gmb();
        // } else if _step == 150000 {
        //     ctrl.state.gmb_d.roll = 0.02;
        //     ctrl.state.gmb_d.pitch = -0.69;
        //     ctrl.state.gmb_d.yaw = -0.02;
        //     ctrl.state.gmb_d.calculate_rotation_matrix();
        //     ctrl.state.update_desired_eq_from_gmb();
        // }

        let start_time = Instant::now();
        unsafe {
            trace!("bit_one_step start");
            bit_one_step(
                bp.x.as_ptr(),
                tau_applied.as_mut_ptr(),
                bp.unlock.as_ptr(),
                1.0*ctrl.pivot.omega_request/1.0,
                true as u8,
                bp._dt,
                bp._num_steps,
                bp._tau_piv_max,
                bp.pitch_nom,
                flex.eta.as_ptr(),
                fc.u.as_ptr(),
                true as u8,
                false as u8,//True means SB
                y_result.as_mut_ptr(),
                flex_result.as_mut_ptr(),
            );
            // Assuming bp.x, flex.eta, and fc.u are Vec<f64> or similar
        // Convert to Vec<f32>
        // let bp_x_f32: Vec<f32> = bp.x.iter().map(|&x| x as f32).collect();
        // let flex_eta_f32: Vec<f32> = flex.eta.iter().map(|&x| x as f32).collect();
        // let fc_u_f32: Vec<f32> = fc.u.iter().map(|&x| x as f32).collect();

        // // Assuming tau_applied, y_result, and flex_result are mutable vectors of f64
        // let mut tau_applied_f32: Vec<f32> = tau_applied.iter().map(|&x| x as f32).collect();
        // let mut bp_unlock_f32: Vec<f32> = bp.unlock.iter().map(|&x| x as f32).collect();
        // let mut y_result_f32: Vec<f32> = y_result.iter().map(|&x| x as f32).collect();
        // let mut flex_result_f32: Vec<f32> = flex_result.iter().map(|&x| x as f32).collect();

        // // Now, call the function with these converted arrays
        // // and other parameters appropriately casted
        // bit_one_step(
        //     bp_x_f32.as_ptr(),
        //     tau_applied_f32.as_mut_ptr(),
        //     bp_unlock_f32.as_ptr(), // Assuming this is already f32 or does not need conversion
        //     (ctrl.pivot.omega_request / 1.0) as f32,
        //     true as u8,
        //     bp._dt as f32,
        //     bp._num_steps, // Assuming this is an integer type and does not need conversion
        //     bp._tau_piv_max as f32,
        //     bp.pitch_nom as f32,
        //     flex_eta_f32.as_ptr(),
        //     fc_u_f32.as_ptr(),
        //     true as u8,
        //     y_result_f32.as_mut_ptr(),
        //     flex_result_f32.as_mut_ptr(),
        // );

            trace!("bit_one_step end");
        }
        
        let duration = start_time.elapsed();
        total_duration += duration;   
        
        if flex.flex_enable{
            // flex.propogate_flex(&[tau_applied[6] + fc.u_rigid[0], tau_applied[7]+fc.u_rigid[1], tau_applied[8]+fc.u_rigid[2]], bp._dt, bp._num_steps, &fc);
        // flex.propogate_flex(&[1., 1., 1.], bp._dt, bp._num_steps);
            flex.update_eta(flex_result.as_ref());
        }
        
        // println!("flex {}", now2.elapsed().as_micros());

        // println!("Pivot request {:}", ctrl.pivot.omega_request);

        step = step + 1;
        t = t + 0.001;
        ctrl.error.current_time = t;
        
        late_meas = meas.clone();
        meas.gps.utc = meas.gps.utc + chrono::Duration::milliseconds(1);
        sim_state.gps = meas.gps.clone();

        // reassign the state matrix to the last state matrix
        //send theta to bp for use in calculating the angular velocity vector
        bp.x = Vector21::from_row_slice(y_result.as_ref());
        bp.update_state();
        bp.get_orientation_vec(&flex);
        bp.get_orientation_rots(&flex, &mut meas);
        bp.get_omega_meas();
        bp.get_rw_speed();

        // if &ctrl.slew_flag == &true {
        //     bp.pitch_nom = bp.theta[8].clone();
        //     println!("pitch nom {}", bp.pitch_nom);
        // }
        // println!("theta nom {}",&bp.pitch_nom);

        sim_state.update_current_equatorial_coordinates(&bp.phi_act);

        //Meas is the struct holding the actual gyro calibration terms
        meas.cbh = sim_state.hor.rot.clone();
        meas.read_measurements(&bp, &sim_state, &flex);

        if (step % 2) < 1{
            if bp.latency {
                est.gyros_bs.read_gyros(late_meas.gyros_bs.omega_m, t.clone());
                ctrl.state.gps = late_meas.gps.clone();
            } else {
                est.gyros_bs.read_gyros(meas.gyros_bs.omega_m, t.clone());
                ctrl.state.gps = meas.gps.clone();
            }

            est.piv_speed = ctrl.lqr.control_input[3];

            est.propogate();
        }

        ctrl.state.omega_gond = na::Vector5::<f64>::from_row_slice(&[
            late_meas.gyro_of.om_axial,
            late_meas.gyro_bow.om_axial,
            late_meas.gyro_stern.om_axial,
            late_meas.gyro_port.om_axial,
            late_meas.gyro_sb.om_axial
        ]);

        if fc.enable{
            // need to input the gyro data from rigid modes
            
            let mut gyro_in = [
                late_meas.gyro_of.om_axial,
                late_meas.gyro_bow.om_axial,
                late_meas.gyro_stern.om_axial,
                late_meas.gyro_port.om_axial,
                // fc.u[3],
                // bp.omega_rw * 4.5
                late_meas.gyro_sb.om_axial
            ];

            if bp.latency {
                 
            } else {
                gyro_in = [
                meas.gyro_of.om_axial,
                meas.gyro_bow.om_axial,
                meas.gyro_stern.om_axial,
                meas.gyro_port.om_axial,
                // fc.u[3],
                // bp.omega_rw * 4.5
                meas.gyro_sb.om_axial
            ]; 
            }

            // println!("{} {} {} {} {}", &gyro_in[0], &gyro_in[1], &gyro_in[2], &gyro_in[3], &gyro_in[4]);
            // if &ctrl.slew_flag==&true{
            //     fc.propogate_control_state(gyro_in.as_slice(), bp._dt, bp._num_steps, &[ctrl.error.rate_des.z.clone(), ctrl.error.rate_des.x.clone(),ctrl.error.rate_des.x.clone(),ctrl.error.rate_des.y.clone(),ctrl.error.rate_des.y.clone(),]);
            // } else {
                let sc: f64 = 1.0;
                fc.err_weight = ctrl.error.err_weight.clone();
                fc.pitch = (&ctrl.state.gmb_d.pitch - &bp.pitch_nom);
                fc.pitch_nom = bp.pitch_nom;
                // println!("gmb_d pitch {} pitch nom {} calc: {}",&ctrl.state.gmb_d.pitch, &bp.pitch_nom, (&ctrl.state.gmb_d.pitch - &bp.pitch_nom));
                // println!("err weight {}", fc.err_weight);
                // println!("fc call: t:{}", &t);
                // if fc.err_weight < 0.05{
                    fc.propogate(gyro_in.as_slice(), bp._dt, bp._num_steps, &[sc*ctrl.error.rate_des.z.clone(), sc*ctrl.error.rate_des.x.clone(),sc*ctrl.error.rate_des.y.clone()]);
                    
                // }
            // }
            
        }

        let centroid_on = ctrl.error.err_comb_th.norm() < 1e-5;

        if (step % 1000) < 1 || (step % 50) == 0 && centroid_on {

            // if (step % 100) == 0 && centroid_on{
            //     println!("Centroiding")
            // }
            // Your code here
        
            est.corr.read_LIS(&sim_state.eq_k.rot);
            if {est.reset | first_LIS} {
                est.eq_hat_k.rot = sim_state.eq_k.rot.clone();
                est.gyros_bs.reset();
                est.reset = false;
                first_LIS = false;
            }
            est.correct_estimate(); 
            // println!("Correction t: {} ",t);
        }        
        
        if bp.ctrl_from_est {
            ctrl.state.eq_k = est.eq_hat_k.clone();
            ctrl.state.update_current_horizontal_coordinates();           
            
            ctrl.state.omega[0] = est.gyros_bs.omega_k[0].clone();//omega_k in est is the corrected gyro.
            ctrl.state.omega[1] = est.gyros_bs.omega_k[1].clone();//omega_m in est is the measured gyro
            ctrl.state.omega[2] = est.gyros_bs.omega_k[2].clone();
            if bp.latency {
                ctrl.state
                    .gmb_k
                    .update_gimbal_coordinates(&[late_meas.roll, late_meas.pitch, late_meas.yaw_p]);
            } else {
                ctrl.state
                    .gmb_k
                    .update_gimbal_coordinates(&[meas.roll, meas.pitch, meas.yaw_p]);
            }
            ctrl.state.piv_angle = est.piv_angle;
            ctrl.state.piv_speed = ctrl.pivot.omega_request;
            ctrl.state.rw_hs = bp.x[20];
            ctrl.state.rw_hs_nom = ctrl.pivot.omega_rw_nom*ctrl.pivot._i_rw;
            fc.roll = ctrl.state.gmb_d.roll;
            fc.pitch = ctrl.state.gmb_d.pitch-bp.pitch_nom;
        } else {
            ctrl.state.eq_k = sim_state.eq_k.clone();
            ctrl.state.omega[0] = bp.omega_m[0].clone();
            ctrl.state.omega[1] = bp.omega_m[1].clone();
            ctrl.state.omega[2] = bp.omega_m[2].clone();
            ctrl.state
            .gmb_k
            .update_gimbal_coordinates(&[bp.x[16], bp.x[17], -sim_state.hor.az.clone()]);
            fc.roll = ctrl.state.gmb_d.roll;
            fc.pitch = ctrl.state.gmb_d.pitch-bp.pitch_nom;
        }

        // grab_vec3(&mut ctrl.state.omega, &bp.omega_m);
        
        ctrl.rw.omega = bp.omega_rw;
        

        // ctrl.update_ctrl();


        // Disturbance generation
        // wind.generate();
        // tau_applied[0] = wind.tau[0].clone() - bp.damp[0] * bp.d_theta_dt[0];
        // tau_applied[1] = wind.tau[1].clone() - bp.damp[1] * bp.d_theta_dt[1];
        // tau_applied[2] = wind.tau[2].clone() - bp.damp[2] * bp.d_theta_dt[2];


        // update actual torque vector with applied torques

        // if step < 1500{
        ctrl.update_ctrl();

        if servo.enable{
            if fc.ff_flag{
                servo.propogate(&ctrl.error.err_comb_th, &ctrl.error.err_rate_lqr);
                
            }
            
        }
// 
        if (step % 10) < 1 {  
                // ctrl.pivot.calculate_pivot_speed(&ctrl.rw, &fc.u[0]);
                // if ctrl.error.err_comb_th.norm() > 0.01{
                //     fc.reset();
                // }
                fc.get_control_input();

                if servo.enable{
                    servo.get_control_input(&ctrl.error.err_comb_th, &ctrl.error.err_rate);
                }

                
                if ctrl.error.err_comb_th.norm() < 0.001 && (servo.enable==false) {
                    // println!("IT HAPPENED");
                    servo.enable=true;
                    servo_multiplier = 1.0;

                } 
                servo.enable=true;
                    servo_multiplier = 1.0;
                // } else {
                //     println!("Norm: {}", ctrl.error.err_comb_th.norm() );
                // }
                // ctrl.pivot.calculate_pivot_speed(&ctrl.rw, &fc.u[0]);
                // tau_applied[6] = 1.0*ctrl.rw.tau_applied; //yaw
                // tau_applied[7] = 1.0*ctrl.fmot_roll.tau_applied; //roll
                // tau_applied[8] = 1.0*ctrl.fmot_pitch.tau_applied; //pitch
                tau_req_yaw = ctrl.rw.tau_applied - servo_multiplier*servo.tau_y_lim; //yaw
                tau_req_roll = ctrl.fmot_roll.tau_applied + fc.ff_r + servo_multiplier*servo.tau_r_lim; //roll
                tau_req_pitch = ctrl.fmot_pitch.tau_applied + (1.0*fc.ff_p) + 1.0*servo_multiplier*servo.tau_p_lim; //pitch
                // tau_req_yaw = 0.0; //yaw
                // tau_req_roll = fc.ff_r; //roll
                // tau_req_pitch = (1.0*fc.ff_p); //pitch
                // println!("pitch ff{}", &fc.ff_p);
                // fc.u = fc.u;
                // println!("yaw requested {} roll req: {} pitch req: {}", tau_req_yaw, tau_req_roll, tau_req_pitch);
                // let diff_yaw = tau_req_yaw - tau_applied[6];
                // let diff_roll = tau_req_roll - tau_applied[7];
                // let diff_pitch = tau_req_pitch - tau_applied[8];
                

                // if diff_yaw.abs() > 0.2 {
                //     tau_applied[6] = tau_applied[6] + 0.2*diff_yaw.signum();
                // } else {
                    tau_applied[6] = tau_req_yaw;
                // }

                // let thrsh = 0.5;

                // if diff_roll.abs() > thrsh {
                //     tau_applied[7] = tau_applied[7] + thrsh*diff_roll.signum();
                // } else {
                    tau_applied[7] = tau_req_roll;
                // }

                // if diff_pitch.abs() > thrsh {
                //     tau_applied[8] = tau_applied[8] + thrsh*diff_pitch.signum();
                // } else {
                    tau_applied[8] = tau_req_pitch;
                // }

                // if tau_applied[6].abs() > 10.0{
                //     tau_applied[6] = tau_applied[6].signum() * 10.0;
                // }
                // if tau_applied[7].abs() > 10.0{
                //     tau_applied[7] = tau_applied[7].signum() * 10.0;
                // }
                // if tau_applied[8].abs() > 10.0{
                //     tau_applied[8] = tau_applied[8].signum() * 10.0;
                // }

                let mut ctrl_rw_act = ctrl.rw.clone();
                ctrl_rw_act.tau_applied = tau_applied[6].clone();

                ctrl.pivot.calculate_pivot_speed(&ctrl_rw_act, &0.0);

                // tau_applied[6] = fc.u[0] + ctrl.rw.tau_applied; //yaw
                // tau_applied[7] = fc.u[1] + ctrl.fmot_roll.tau_applied; //roll
                // tau_applied[8] = fc.u[2] + ctrl.fmot_pitch.tau_applied; //pitch
                // println!("{} \n{} \n{}\n {}\n ", tau_applied[6], diff_yaw, fc.u[0], &ctrl.rw.tau_applied);

                // println!("yaw:{} \nrol: {} \n pit: {}\n rw: {}\n bow: {}\nstn: {}\nprt: {}\n sb: {}\n t: {}\n", tau_applied[6], tau_applied[7], tau_applied[8], fc.u[0], fc.u[1], fc.u[2], fc.u[3], fc.u[4], &t);
                
            // } else {
            //     ctrl.update_ctrl();
            //     ctrl.rw.tau_applied = 0.0;
            //     tau_applied[6] = 0.0*ctrl.rw.tau_applied; //yaw
            //     tau_applied[7] = 0.0*ctrl.fmot_roll.tau_applied; //roll
            //     tau_applied[8] = 0.0*ctrl.fmot_pitch.tau_applied; 
            // }
            

            // if (step % 1000) < 1{
            //     println!("step: {:}", step);
            // }
        }
        


        sc::push_record(&t, &bp, &est, &ctrl, &meas, &sim_state, &flex, &fc, &servo).unwrap();
        js::read_gains(&mut ctrl, &mut bp, &mut est, &mut fc, &mut flex); // read in gains from json file (for tuning)

        
        
        // if fc.update | (_step == 0){
        //     println!("Success");
        //     let mut new_fc = flex_control::DampingControl::init_pc(&ctrl.fine_gains);
        //     let mut fc = new_fc.clone();

        // }

        // println!("RW Applied: {} Request{}", tau_applied[6], tau_req_yaw);

    }
    println!("bit one step {}", now1.elapsed().as_secs_f64());

    let total_seconds = total_duration.as_secs_f64(); // Total time in seconds
    println!("Total time spent on the function: {:.2} seconds", total_seconds);
    trace!("END");
}
