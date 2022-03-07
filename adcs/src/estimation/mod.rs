// import crate for matrix math
extern crate nalgebra as na;
pub use na::{
    Matrix3, Matrix6, MatrixMN, Vector3, Vector6, U1, U10, U18, U2, U21, U3, U36, U6, U9,
};
// use std::fmt;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

//define types
pub type Vector9 = MatrixMN<f64, U9, U1>;
pub type Vector18 = MatrixMN<f64, U18, U1>;

// Defining a struct to hold the definitions of the system
// matrices
#[derive(Debug)]
pub struct Est {
    // //
    pub _q_k: Matrix3<f64>,
    pub _r_k: Matrix3<f64>,
    // //
    pub omega_k: Vector3<f64>,
    pub xhat_k: Vector18,
    pub c_hat_k: Matrix3<f64>,
    pub xdes: Vector3<f64>,
    pub c_des: Matrix3<f64>,
    pub phi_hat_k: Vector3<f64>,
}
impl Est {
    pub fn init() -> Est {
        // angular velocity measurement (gyros)
        let omega_k = Vector3::<f64>::new(0.0, 0.0, 0.0);
        // current orientation estimate
        let phi_hat_k = Vector3::<f64>::new(0.0, 0.0, 0.);
        let xhat_k = Vector18::zeros();
        let c_hat_k = Matrix3::<f64>::identity();
        let xdes = Vector3::<f64>::new(0.0, 0.0, 0.);
        let c_des = Matrix3::<f64>::identity();

        let _q_k = Matrix3::<f64>::identity();
        let _r_k = Matrix3::<f64>::identity();
        // DEFINE SYSTEMATIC TERMS
        // ----------------------------------

        // define params struct
        let est = Est {
            omega_k,
            phi_hat_k,
            xhat_k,
            c_hat_k,
            xdes,
            c_des,
            _q_k,
            _r_k,
        };
        info!("Estimation initialized: \n {:?}", est);
        return est;
    }
}
