/// This is a module for writing tests for the ADCS code
/// 
/// tests are performed against the live matlab script located in the flight code workspace. rust_control_verification.mlx
use crate::control::*;
extern crate nalgebra as na;

pub fn test_control(){
    // test the equitorial functionality
    let mut eq = state::equitorial::Equatorial::new();
    eq.fr = 0.45;
    eq.dec = -0.7;
    eq.ra = 1.5;

    eq.calculate_rotation_matrix();

    let res = na::Rotation3::<f64>::from_matrix(&na::Matrix3::<f64>::from_row_slice(&[0.0541, 0.7629, -0.6442, -0.8784, 0.3432, 0.3327, 0.4749, 0.5479,0.6887]));
    
    //the result should be identity
    println!("{}^T * {} = {}", res, eq.rot, res.inverse()*eq.rot);

    let mut eq2 = state::equitorial::Equatorial::new();
    eq2.rot = eq.rot.clone();

    eq2.extract_ra_dec_fr_from_rotmat();
    // these results should be 0
    println!("{} - {} = {}", eq.fr, eq2.fr, eq.fr-eq2.fr);
    println!("{} - {} = {}", eq.dec, eq2.dec, eq.dec-eq2.dec);
    println!("{} - {} = {}", eq.ra, eq2.ra, eq.ra-eq2.ra);
    // 
    //
    // test the gimbal functionality
    let mut gmb = state::gimbal::Gimbal::new();
    gmb.roll = -1.45;
    gmb.pitch = -0.14;
    gmb.yaw = 1.8216;

    gmb.calculate_rotation_matrix();

    let res = na::Rotation3::<f64>::from_matrix(&na::Matrix3::<f64>::from_row_slice(&[-0.3799, 0.9249, 0.0168, -0.1167, -0.0299, -0.9927, -0.9176, -0.3791, 0.1193]));
    
    //the result should be identity
    println!("{}^T * {} = {}", res, gmb.rot, res.inverse()*gmb.rot);

    let mut gmb2 = state::gimbal::Gimbal::new();
    gmb2.rot = gmb.rot.clone();

    gmb2.extract_gimbal_rpy();
    // these results should be 0
    println!("{} - {} = {}", gmb.roll, gmb2.roll, gmb.roll-gmb2.roll);
    println!("{} - {} = {}", gmb.pitch, gmb2.pitch, gmb.pitch-gmb2.pitch);
    println!("{} - {} = {}", gmb.yaw, gmb2.yaw, gmb.yaw-gmb2.yaw);

    //test gimbal mapping matrix
    gmb.calculate_gimbal_mapping_matrix();
    gmb2.gmm = na::Matrix3::<f64>::from_row_slice(&[
        0.990215996212637,
        0.0,
        0.016815331760778,
        0.0,
        1.000000000000000,
        -0.992712991037588,
        -0.139543114644236,
        0.0,
        0.119323769815489]);
    println!("gmm: {}, gmm_mat {}, gmm-gmm_mat {:.6}", gmb.gmm, gmb2.gmm, gmb.gmm-gmb2.gmm);

    //
    //
    //
    //
    // test the equitorial functionality
    let mut hor = state::horizontal::Horizontal::new();
    hor.ir = -1.5;
    hor.el = -1.470796326794896;
    hor.az = 0.75;

    hor.calculate_rotation_matrix();

    let res = na::Rotation3::<f64>::from_matrix(&na::Matrix3::<f64>::from_row_slice(&[
         0.073046999702127,  -0.068050326332037,  -0.995004165278026,
        -0.677992520184125,   0.728292044747507,  -0.099583332600765,
         0.731430296343333,   0.681879645277375,   0.007061936526523]));
         
    
    //the result should be identity
    println!("{}^T * {} = {}", res, hor.rot, res.inverse()*hor.rot);

    let mut hor2 = state::horizontal::Horizontal::new();
    hor2.rot = hor.rot.clone();

    hor2.extract_az_el_ir_from_rotmat();
    // these results should be 0
    println!("{} - {} = {}", hor.ir, hor2.ir, hor.ir-hor2.ir);
    println!("{} - {} = {}", hor.el, hor2.el, hor.el-hor2.el);
    println!("{} - {} = {}", hor.az, hor2.az, hor.az-hor2.az);
    //
    //
    //
    //
    // Test the state module
    //
    //
    //
    let mut st = state::State::new();
    st.gmb_k = gmb;

    let _phi = st.calculate_coupling_matrix();

    st.gmb_k.calculate_gimbal_mapping_matrix();
    let temp = st.gmb_k.gmm;

    let _phi2 = st.gmb_k.gmm.transpose() * temp;

    println!("coupling{}, gmm^T gmm{}",_phi, _phi2);

    //
    //
    //C_HE check
    // Note, have to run Matlab at same time as this to verify correctness
    // Also note, verified not accurately against maltab but accurately 
    // againt website 'http://neoprogrammics.com/sidereal_time_calculator/index.php'
    st.update_hor_to_eq_conversion();
    println!("GAST {} C_EH {}", st.gps.gast,st.ceh);

    //hor to eq conversion
    st.update_current_equatorial_coordinates(&[hor.ir, hor.el, hor.az]);
    println!("{} {} {}\n ra {} dec {} fr {}", st.hor.ir, st.hor.el, st.hor.az, st.eq_k.ra, st.eq_k.dec, st.eq_k.fr);

    // eq to gimbal conversion
    st.eq_d = eq.clone();
    st.update_desired_gimbal_rpy();
    println!("roll {} pitch {} yaw {}", st.gmb_d.roll, st.gmb_d.pitch, st.gmb_d.yaw);

    // gimbal to eq conversion
    let mut gmb3 = state::gimbal::Gimbal::new();
    gmb3.roll = -1.45;
    gmb3.pitch = -0.14;
    gmb3.yaw = 1.8216;
    st.gmb_d = gmb3.clone();

    st.update_desired_eq_from_gmb();
    println!("desired: ra {} dec {} fr {}", st.eq_d.ra, st.eq_d.dec, st.eq_d.fr);

}



