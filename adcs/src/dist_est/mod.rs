// import crate for matrix math
pub mod mats;

extern crate nalgebra as na;

type Matrix10 = na::SMatrix<f64, 10, 10>;
type Vector10 = na::SVector<f64, 10>;
type Matrix3x10 = na::SMatrix<f64, 3, 10>;
type Matrix6x10 = na::SMatrix<f64, 6, 10>;
type Vector3 = na::SVector<f64, 3>;

/// Struct for Kalman Filter
#[derive(Debug)]
pub struct KalmanFilter {
    pub a: Matrix10,          // State transition matrix
    pub c_vel: Matrix3x10,    // Measurement matrix for velocity
    pub c_pos: Matrix3x10,    // Measurement matrix for position
    pub c_combined: Matrix6x10,
    pub p: Matrix10,          // Estimated state covariance
    pub r_vel: na::Matrix3<f64>, // Measurement noise covariance for velocity
    pub r_pos: na::Matrix3<f64>, // Measurement noise covariance for position
    pub r_combined: na::Matrix6<f64>,
    pub q: Matrix10,          // Process noise covariance
    pub x_hat: Vector10,      // Estimated state
}

impl KalmanFilter {
    /// Function to instantiate a new Kalman Filter
    // pub fn new() -> KalmanFilter {
    //     let a = mats::init_AMat1_pos();
    //     let c_vel = mats::init_AMat2_pos();
    //     let c_pos = mats::init_AMat3_pos();
    //     let r_vel = na::Matrix3::<f64>::identity()* 1.0e-5;
    //     let r_pos = na::Matrix3::<f64>::identity()* 1.0e-5;
    //     let q = Matrix10::identity() * 1.0e-5;
    //     KalmanFilter {
    //         a,
    //         c_vel,
    //         c_pos,
    //         p: Matrix10::identity(), // Initial covariance
    //         r_vel,
    //         r_pos,
    //         q,
    //         x_hat: Vector10::zeros(), // Initial state estimate
    //     }
    // }
    pub fn new() -> KalmanFilter {
        let a = mats::init_AMat1_pos();
        let c_vel = mats::init_AMat2_pos();
        let c_pos = mats::init_AMat3_pos();
        
        // Initialize noise covariance matrices for velocity and position
        let r_vel = na::Matrix3::<f64>::identity() * 1.0e-2;
        let r_pos = na::Matrix3::<f64>::identity() * 1.0e-2;
    
        // Combine measurement matrices (c_vel and c_pos)
        let c_combined = Matrix6x10::from_rows(&[
            c_vel.row(0),
            c_vel.row(1),
            c_vel.row(2),
            c_pos.row(0),
            c_pos.row(1),
            c_pos.row(2),
        ]);
    
        // Combine noise covariance matrices (r_vel and r_pos)
        let r_combined = na::Matrix6::from_diagonal(&na::Vector6::new(
            r_vel[(0,0)], r_vel[(1,1)], r_vel[(2,2)],
            r_pos[(0,0)], r_pos[(1,1)], r_pos[(2,2)],
        ));
    
        let q = Matrix10::identity() * 1.0e-5;
    
        KalmanFilter {
            a,
            c_vel,
            c_pos,
            c_combined, // Combined measurement matrix
            p: Matrix10::identity(), // Initial covariance
            r_vel,
            r_pos,
            r_combined, // Combined noise covariance matrix
            q,
            x_hat: Vector10::zeros(), // Initial state estimate
        }
    }
    
    /// Propagation step of the Kalman Filter
    pub fn propagate(&mut self, dt: f64) {
        // Propagate the state forward in time
        self.x_hat = &self.x_hat + (&self.a * &self.x_hat) * dt;

        // Predict the state covariance
        self.p = &self.a * &self.p * self.a.transpose() * dt.powi(2) + &self.q;
    }


    // /// Correction step with velocity measurement
    // pub fn correct_with_velocity(&mut self, velocity_measurement: &[f64;3]) {

    //     let y_tilde = na::Vector3::<f64>::from_row_slice(velocity_measurement) - &self.c_vel * &self.x_hat;
    //     let s = &self.c_vel * &self.p * self.c_vel.transpose() + &self.r_vel;
    //     let k = &self.p * self.c_vel.transpose() * s.try_inverse().expect("S inverse failed");

    //     self.x_hat = &self.x_hat + k * y_tilde;
    //     let identity: Matrix10 = Matrix10::identity();
    //     self.p = (identity - k * self.c_vel) * self.p;
    // }

    // /// Correction step with positional measurement
    // pub fn correct_with_position(&mut self, position_measurement: &[f64;3]) {
    //     let y_tilde = na::Vector3::<f64>::from_row_slice(position_measurement) - &self.c_pos * &self.x_hat;
    //     let s = &self.c_pos * &self.p * self.c_pos.transpose() + &self.r_pos;
    //     let k = &self.p * self.c_pos.transpose() * s.try_inverse().expect("S inverse failed");

    //     self.x_hat = &self.x_hat + k * y_tilde;
    //     let identity: Matrix10 = Matrix10::identity();
    //     self.p = (identity - k * self.c_pos) * self.p;
    // }

    /// Combined correction step with both velocity and positional measurements
    pub fn correct_with_both_measurements(&mut self, velocity_measurement: &[f64;3], position_measurement: &[f64;3]) {
        // Augment measurement vector
        let y_vel = na::Vector3::<f64>::from_row_slice(velocity_measurement);
        let y_pos = na::Vector3::<f64>::from_row_slice(position_measurement);
        let y_tilde = na::Vector6::new(y_vel[0], y_vel[1], y_vel[2], y_pos[0], y_pos[1], y_pos[2]) 
                    - &self.c_combined * &self.x_hat;

        // Augment measurement matrix and noise covariance
        let s = &self.c_combined * &self.p * self.c_combined.transpose() + &self.r_combined;
        let k = &self.p * self.c_combined.transpose() * s.try_inverse().expect("S inverse failed");

        // Update state estimate and covariance
        self.x_hat = &self.x_hat + k * y_tilde;
        let identity: Matrix10 = Matrix10::identity();
        self.p = (identity - k * self.c_combined) * self.p;
    }
}