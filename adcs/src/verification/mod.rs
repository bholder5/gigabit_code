/// This is a module for writing tests for the ADCS code
///
/// tests are performed against the live matlab script located in the flight code workspace. rust_control_verification.mlx
use crate::control::{self, *};
use crate::estimation::{self};
extern crate nalgebra as na;

type Matrix18 = na::SMatrix<f64, 18, 18>;

pub fn test_control() {
    // test the equitorial functionality
    let mut eq = state::equitorial::Equatorial::new();
    eq.fr = 0.45;
    eq.dec = -0.7;
    eq.ra = 1.5;

    eq.calculate_rotation_matrix();

    let res = na::Rotation3::<f64>::from_matrix(&na::Matrix3::<f64>::from_row_slice(&[
        0.0541, 0.7629, -0.6442, -0.8784, 0.3432, 0.3327, 0.4749, 0.5479, 0.6887,
    ]));

    //the result should be identity
    println!("{}^T * {} = {}", res, eq.rot, res.inverse() * eq.rot);

    let mut eq2 = state::equitorial::Equatorial::new();
    eq2.rot = eq.rot.clone();

    eq2.extract_ra_dec_fr_from_rotmat();
    // these results should be 0
    println!("{} - {} = {}", eq.fr, eq2.fr, eq.fr - eq2.fr);
    println!("{} - {} = {}", eq.dec, eq2.dec, eq.dec - eq2.dec);
    println!("{} - {} = {}", eq.ra, eq2.ra, eq.ra - eq2.ra);
    //
    //
    // test the gimbal functionality
    let mut gmb = state::gimbal::Gimbal::new();
    gmb.roll = -1.45;
    gmb.pitch = -0.14;
    gmb.yaw = 1.8216;

    gmb.calculate_rotation_matrix();

    let res = na::Rotation3::<f64>::from_matrix(&na::Matrix3::<f64>::from_row_slice(&[
        -0.3799, 0.9249, 0.0168, -0.1167, -0.0299, -0.9927, -0.9176, -0.3791, 0.1193,
    ]));

    //the result should be identity
    println!("{}^T * {} = {}", res, gmb.rot, res.inverse() * gmb.rot);

    let mut gmb2 = state::gimbal::Gimbal::new();
    gmb2.rot = gmb.rot.clone();

    gmb2.extract_gimbal_rpy();
    // these results should be 0
    println!("{} - {} = {}", gmb.roll, gmb2.roll, gmb.roll - gmb2.roll);
    println!(
        "{} - {} = {}",
        gmb.pitch,
        gmb2.pitch,
        gmb.pitch - gmb2.pitch
    );
    println!("{} - {} = {}", gmb.yaw, gmb2.yaw, gmb.yaw - gmb2.yaw);

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
        0.119323769815489,
    ]);
    println!(
        "gmm: {}, gmm_mat {}, gmm-gmm_mat {:.6}",
        gmb.gmm,
        gmb2.gmm,
        gmb.gmm - gmb2.gmm
    );

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
        0.073046999702127,
        -0.068050326332037,
        -0.995004165278026,
        -0.677992520184125,
        0.728292044747507,
        -0.099583332600765,
        0.731430296343333,
        0.681879645277375,
        0.007061936526523,
    ]));

    //the result should be identity
    println!("{}^T * {} = {}", res, hor.rot, res.inverse() * hor.rot);

    let mut hor2 = state::horizontal::Horizontal::new();
    hor2.rot = hor.rot.clone();

    hor2.extract_az_el_ir_from_rotmat();
    // these results should be 0
    println!("{} - {} = {}", hor.ir, hor2.ir, hor.ir - hor2.ir);
    println!("{} - {} = {}", hor.el, hor2.el, hor.el - hor2.el);
    println!("{} - {} = {}", hor.az, hor2.az, hor.az - hor2.az);
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
    st.gmb_k.calculate_inverse_gimbal_mapping_matrix();
    let temp = st.gmb_k.gmm;

    let _phi2 = st.gmb_k.gmm.transpose() * temp;

    println!("coupling{}, gmm^T gmm{}", _phi, _phi2);
    println!(
        "gmm_i {} gmm^_1 {}",
        st.gmb_k.gmm_i,
        st.gmb_k.gmm.try_inverse().unwrap()
    );

    //
    //
    //C_HE check
    // Note, have to run Matlab at same time as this to verify correctness
    // Also note, verified not accurately against maltab but accurately
    // againt website 'http://neoprogrammics.com/sidereal_time_calculator/index.php'
    st.update_hor_to_eq_conversion();
    println!("GAST {} C_EH {}", st.gps.gast, st.ceh);

    //hor to eq conversion
    st.update_current_equatorial_coordinates(&[hor.ir, hor.el, hor.az]);
    println!(
        "{} {} {}\n ra {} dec {} fr {}",
        st.hor.ir, st.hor.el, st.hor.az, st.eq_k.ra, st.eq_k.dec, st.eq_k.fr
    );

    // eq to gimbal conversion
    st.eq_d = eq.clone();
    st.update_desired_gimbal_rpy();
    println!(
        "roll {} pitch {} yaw {}",
        st.gmb_d.roll, st.gmb_d.pitch, st.gmb_d.yaw
    );

    // gimbal to eq conversion
    let mut gmb3 = state::gimbal::Gimbal::new();
    gmb3.roll = -1.45;
    gmb3.pitch = -0.14;
    gmb3.yaw = 1.8216;
    st.gmb_d = gmb3.clone();

    st.update_desired_eq_from_gmb();
    println!(
        "desired: ra {} dec {} fr {}",
        st.eq_d.ra, st.eq_d.dec, st.eq_d.fr
    );
    //
    //
    //
    //
    //     Testing error calculations
    //
    //
    //
    //
    println!("\n\n\n TESTING POSITION ERRORS \n\n\n");
    let mut st = state::State::new();
    st.gmb_k.roll = -0.02;
    st.gmb_k.pitch = -0.35;
    st.gmb_k.yaw = 0.1;

    st.gmb_k.calculate_rotation_matrix();
    st.gmb_k.calculate_gimbal_mapping_matrix();

    st.gmb_d.roll = -0.02;
    st.gmb_d.pitch = -0.3501;
    st.gmb_d.yaw = 0.09;
    st.gmb_d.calculate_rotation_matrix();

    st.update_desired_eq_from_gmb();
    st.update_current_eq_from_gmb();

    let mut er = error::Error::new();
    er.update_pointing_positional_error(&st, na::Matrix3::identity());
    println!("{} {} {}", st.eq_d.ra, st.eq_d.dec, st.eq_d.fr);
    println!("eq desired: {:.5} eq k: {:.5}", st.eq_d.rot, st.eq_k.rot);
    println!(
        "mapped: {} {} {} \ngmb: {} {} {} \n combined {}",
        er.err_b_th.roll,
        er.err_b_th.pitch,
        er.err_b_th.yaw,
        er.err_gmb_th.roll,
        er.err_gmb_th.pitch,
        er.err_gmb_th.yaw,
        er.err_comb_th
    );

    er.err_comb_th = na::Vector3::new(-0.10, 2.0, 0.750);
    er.update_pointing_velocity_error_terms(&mut st, &true);
    println!("err_decay {}, crtl_dt {}", er._err_decay, er._ctrl_dt);
    println!("rate: {}", er.err_rate);
    println!("{}", er.err_rate_sum);
    er.update_pointing_velocity_error_terms(&mut st, &true);
    println!("{}", er.err_rate_sum);
    er.update_pointing_velocity_error_terms(&mut st, &true);
    println!("{}", er.err_rate_sum);

    //////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    ////
    ////
    ////
    //// Testing the ctrol update functions

    let mut ct = control::Ctrl::new();
    ct.error.err_rate = na::Vector3::new(-1000.0, -2000.0, -300.0);
    ct.error.err_rate_sum = na::Vector3::new(3.0, 4.0, 5.0);
    // have to make update_request_torque a public function for this test.
    ct.calculate_applied_torque();
}

pub fn test_estimation() {
    let mut gyro = estimation::gyros::GyroBs::new();
    let bias = na::Vector3::<f64>::new(0.000001, -0.000005, 0.0000075);
    // ax = [1; 1.5; -0.75];
    // ax = ax/norm(ax);
    // Abg_true= diag([1.005, 0.995, 1.0025]) + axis2rot(ax, 0.00025) - eye(3);
    let ax = na::Unit::<na::Vector3<f64>>::new_normalize(na::Vector3::<f64>::new(1.0, 1.5, -0.75));
    let Abg_true = na::Matrix3::<f64>::from_diagonal(&na::Vector3::new(1.005, 0.995, 1.0025))
        + na::Matrix3::<f64>::from(na::Rotation3::<f64>::from_axis_angle(&ax, 0.00025).transpose())
        - na::Matrix3::identity();
    // cbg_axis = [0.7; -0.2; 1.2]
    // cbg_axis = cbg_axis/norm(cbg_axis);
    // Cbg_true = axis2rot(cbg_axis, 0.2);
    let cbg_axis =
        na::Unit::<na::Vector3<f64>>::new_normalize(na::Vector3::<f64>::new(0.7, -0.2, 1.2));
    let cbg_true = na::Rotation3::<f64>::from_axis_angle(&cbg_axis, 0.2).transpose();

    let omega = na::Vector3::<f64>::new(0.1, -0.05, 0.150);

    gyro.bias = bias;
    gyro.a_g = Abg_true;
    gyro.c_bg = cbg_true;
    gyro.omega_m = omega;

    gyro.deconstruct_measurement();
    // println!("{} {} {} {} {}", gyro.omega_k, gyro.bias, gyro.c_bg, gyro.a_g, gyro.omega_m);
    println!("{}", gyro.omega_k);

    let omega2 = na::Vector3::<f64>::new(0.09, -0.06, 0.140);

    gyro.read_gyros(omega2, 0.001);
    println!("{} {} {} ", gyro.om_b, gyro.omega_m, gyro.t1);

    // Testing the propogation
    let mut est = estimation::Estimator::new();
    est.gyros_bs = gyro;

    est.eq_hat_k.rot = na::Rotation3::from_matrix(&na::Matrix3::<f64>::from_row_slice(&[
        0.764559789001659,
        0.065612475901131,
        -0.641204594531155,
        -0.001524258586071,
        0.994986631621922,
        0.099996397582231,
        0.644551010919618,
        -0.075475863027755,
        0.760826779512074,
    ]));

    println!(
        "bias {}, Ag {} Cbg {} omega {} t0 {} t1 {}",
        est.gyros_bs.bias,
        est.gyros_bs.a_g,
        est.gyros_bs.c_bg,
        est.gyros_bs.omega_m,
        est.gyros_bs.t0,
        est.gyros_bs.t1
    );

    est.gyros_bs.bias =
        na::Vector3::new(-0.105740540582954, 0.508595788854436, -0.033077338465802) * 1e-7;
    est.gyros_bs.a_g = na::Matrix3::<f64>::new(
        0.999999999888189,
        0.000000000537792,
        -0.000000000034976,
        0.000000000537792,
        0.999999997413301,
        0.000000000168230,
        -0.000000000034976,
        0.000000000168230,
        0.999999999989059,
    );
    est.gyros_bs.c_bg = na::Rotation3::identity();
    est.gyros_bs.read_gyros(
        na::Vector3::<f64>::new(-0.025344176657539,-0.037471380259213,-0.027217928882192),
        0.002,
    );

    println!(
        "bias {}, Ag {} Cbg {} omega {} t0 {} t1 {}",
        est.gyros_bs.bias,
        est.gyros_bs.a_g,
        est.gyros_bs.c_bg,
        est.gyros_bs.omega_k,
        est.gyros_bs.t0,
        est.gyros_bs.t1
    );

    est.prop.p_hat = Matrix18::from_row_slice(&[
        0.025340249735023,
        -0.000011500990957,
        -0.000012722943938,
        0.000011550010781,
        0.003578969071091,
        -0.003242043702472,
        -0.003562517893483,
        -0.003236063344672,
        -0.003582008195239,
        0.000005980363648,
        0.000008191490983,
        0.000006348460600,
        -0.000003039124526,
        -0.000005201550192,
        -0.000003352048815,
        0.158999054158019,
        -0.000391189617557,
        0.000217951453687,
        -0.000011500990957,
        0.025342425691719,
        -0.000011581288518,
        -0.003590447464977,
        0.000012630012504,
        0.003554302749814,
        -0.000005986271871,
        -0.000008208866989,
        -0.000006355858720,
        -0.003562511615569,
        -0.003236054382976,
        -0.003582001486395,
        0.000006274153779,
        0.000008445978218,
        0.000006642409336,
        0.000391863484855,
        0.158998589846070,
        -0.000407614105910,
        -0.000012722943938,
        -0.000011581288518,
        0.025340070339389,
        0.003229412993322,
        -0.003565853189450,
        0.000011453134522,
        0.000003025144901,
        0.000005173421594,
        0.000003336118591,
        -0.000006279712935,
        -0.000008462892108,
        -0.000006649438810,
        -0.003562517069641,
        -0.003236062437985,
        -0.003582007347514,
        -0.000216702205180,
        0.000408262996815,
        0.158999000587206,
        0.000011550010781,
        -0.003590447464977,
        0.003229412993322,
        0.999999997402519,
        -0.000000000537792,
        0.000000000034976,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000034976,
        0.000000000168230,
        -0.000000000010941,
        -0.000000000537792,
        0.000000002586699,
        -0.000000000168230,
        0.000000000000000,
        -0.000000003307734,
        -0.000000050859579,
        0.003578969071091,
        0.000012630012504,
        -0.003565853189450,
        -0.000000000537792,
        0.999999999877408,
        -0.000000000168230,
        0.000000000034976,
        -0.000000000168230,
        0.000000000010941,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000111811,
        0.000000000537792,
        -0.000000000034976,
        0.000000003307734,
        0.000000000000000,
        -0.000000010574054,
        -0.003242043702472,
        0.003554302749814,
        0.000011453134522,
        0.000000000034976,
        -0.000000000168230,
        0.999999997301650,
        0.000000000537792,
        -0.000000002586699,
        0.000000000168230,
        0.000000000111811,
        -0.000000000537792,
        0.000000000034976,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000050859579,
        0.000000010574054,
        -0.000000000000000,
        -0.003562517893483,
        -0.000005986271871,
        0.000003025144901,
        0.000000000000000,
        0.000000000034976,
        0.000000000537792,
        0.999999999888189,
        0.000000000537792,
        -0.000000000034976,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        -0.000000010574054,
        0.000000000000000,
        -0.000000000000000,
        -0.003236063344672,
        -0.000008208866989,
        0.000005173421594,
        -0.000000000000000,
        -0.000000000168230,
        -0.000000002586699,
        0.000000000537792,
        0.999999997413301,
        0.000000000168230,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        0.000000050859579,
        -0.000000000000000,
        0.000000000000000,
        -0.003582008195239,
        -0.000006355858720,
        0.000003336118591,
        0.000000000000000,
        0.000000000010941,
        0.000000000168230,
        -0.000000000034976,
        0.000000000168230,
        0.999999999989059,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        -0.000000003307734,
        0.000000000000000,
        -0.000000000000000,
        0.000005980363648,
        -0.003562511615569,
        -0.000006279712935,
        -0.000000000034976,
        0.000000000000000,
        0.000000000111811,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        0.999999999888189,
        0.000000000537792,
        -0.000000000034976,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        -0.000000010574054,
        0.000000000000000,
        0.000008191490983,
        -0.003236054382976,
        -0.000008462892108,
        0.000000000168230,
        -0.000000000000000,
        -0.000000000537792,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000537792,
        0.999999997413301,
        0.000000000168230,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        -0.000000000000000,
        0.000000050859579,
        -0.000000000000000,
        0.000006348460600,
        -0.003582001486395,
        -0.000006649438810,
        -0.000000000010941,
        0.000000000000000,
        0.000000000034976,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000034976,
        0.000000000168230,
        0.999999999989059,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        -0.000000003307734,
        0.000000000000000,
        -0.000003039124526,
        0.000006274153779,
        -0.003562517069641,
        -0.000000000537792,
        -0.000000000111811,
        -0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        0.999999999888189,
        0.000000000537792,
        -0.000000000034976,
        -0.000000000000000,
        0.000000000000000,
        -0.000000010574054,
        -0.000005201550192,
        0.000008445978218,
        -0.003236062437985,
        0.000000002586699,
        0.000000000537792,
        0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000537792,
        0.999999997413301,
        0.000000000168230,
        0.000000000000000,
        -0.000000000000000,
        0.000000050859579,
        -0.000003352048815,
        0.000006642409336,
        -0.003582007347514,
        -0.000000000168230,
        -0.000000000034976,
        -0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000034976,
        0.000000000168230,
        0.999999999989059,
        -0.000000000000000,
        0.000000000000000,
        -0.000000003307734,
        0.158999054158019,
        0.000391863484855,
        -0.000216702205180,
        0.000000000000000,
        0.000000003307734,
        0.000000050859579,
        -0.000000010574054,
        0.000000050859579,
        -0.000000003307734,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.999999000001175,
        0.000000000000001,
        -0.000000000000000,
        -0.000391189617557,
        0.158998589846070,
        0.000408262996815,
        -0.000000003307734,
        0.000000000000000,
        0.000000010574054,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000010574054,
        0.000000050859579,
        -0.000000003307734,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        0.000000000000001,
        0.999999000001173,
        0.000000000000000,
        0.000217951453687,
        -0.000407614105910,
        0.158999000587207,
        -0.000000050859579,
        -0.000000010574054,
        -0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000000000000,
        0.000000000000000,
        -0.000000010574054,
        0.000000050859579,
        -0.000000003307734,
        -0.000000000000000,
        0.000000000000000,
        0.999999000001175,
    ]);

    // println!("{:.4}", est.prop.p_hat);

    est.prop.propogate(&mut est.eq_hat_k, &est.gyros_bs);
    println!(
        "psi: {} psi_mat {:.10} fk{:.7} lk {:.15} p_hat {:.7}",
        est.prop.psi_hat, est.prop.psi_hat_mat, est.prop.fk, est.prop.lk, est.prop.p_hat
    );
}
