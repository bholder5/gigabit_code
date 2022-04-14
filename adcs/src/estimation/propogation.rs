//! Struct and functionality to propogate the estimated state

use nalgebra as na;

use crate::control::state::equitorial as eq;
use crate::estimation::gyros;
use crate::miscellaneous as misc;

type Matrix3x9 = na::SMatrix<f64, 3, 9>;
type Matrix18 = na::SMatrix<f64, 18, 18>;
type Matrix18x9 = na::SMatrix<f64, 18, 9>;
type Matrix9 = na::SMatrix<f64, 9, 9>;
/// Struct for propogating the state and covariance
#[derive(Debug)]
pub struct Propogation {
    /// expanded matrix of measured angular velocity + bias
    pub om_b: Matrix3x9,
    /// angular velocity scaled by dt
    pub psi_hat: na::Vector3<f64>,
    /// Rotation matrix representing the change in orientation from the gyro measurement
    pub psi_hat_mat: na::Rotation3<f64>,
    /// Estimation parameters transition matrix
    pub fk: Matrix18,
    /// Process noise input matrix
    pub lk: Matrix18x9,
    /// Estimated parameters covariance
    pub p_hat: Matrix18,
    /// Process noise covariance
    pub _q_k: Matrix9,
}

impl Propogation {
    /// Function to instantiate a new propogation struct
    ///
    /// # Detailed Explanation
    ///
    /// This function instantiates a new propogation struct with default values
    pub fn new() -> Propogation {
        Propogation {
            om_b: Matrix3x9::zeros(),
            psi_hat: na::Vector3::new(0.0, 0.0, 0.0),
            psi_hat_mat: na::Rotation3::identity(),
            fk: Matrix18::zeros(),
            lk: Matrix18x9::zeros(),
            p_hat: Matrix18::zeros(),
            _q_k: Matrix9::identity() * 1e-15,
        }
    }
    /// Function that propogates the estimation
    ///
    /// # Detailed Explanation
    ///
    /// This function propogates the estimation using the following equation
    /// $$ $$
    ///
    /// # Arguments
    ///
    /// `eq_hat: eq::Equatorial` - Estimated equatorial coordinates
    /// `gyro_bs: gyro::Gyro_bs` - Boresight gyroscope measurements
    ///
    /// # Results
    ///
    /// `eq_hat: eq::Equatorial` - Updated estimated equatorial coordinates
    /// `self.P_hat: Matrix18` - Updated estimated covariance
    pub fn propogate(&mut self, eq_hat: &mut eq::Equatorial, gyro_bs: &gyros::GyroBs) {
        // calculate the change in orientation from omega over the time step
        // and get it in the form of a rotation matrix
        let dt = gyro_bs.t1 - gyro_bs.t0;
        self.psi_hat = gyro_bs.omega_k * dt;
        self.psi_hat_mat = na::Rotation3::from_axis_angle(
            &na::Unit::new_normalize(self.psi_hat.clone()),
            self.psi_hat.norm().clone(),
        )
        .inverse();

        // Update current orientation
        eq_hat.rot = self.psi_hat_mat * eq_hat.rot;
        eq_hat.extract_ra_dec_fr_from_rotmat();

        // Build state transition matrix Fk
        // SECOND MATRIX IN Fk
        let mut psi_cross = na::Matrix3::<f64>::identity();
        misc::xmat(&self.psi_hat, &mut psi_cross);

        // THIRD MATRIX IN Fk
        self.om_b = self.calculate_om_b_expmat(&gyro_bs.om_b) * dt;
        let temp_mat_3 = gyro_bs.c_bg * self.om_b;

        // FOURTH MATRIX IN Fk
        let temp_mat_4 = gyro_bs.c_bg * gyro_bs.a_g * dt;

        // Fk = [Psi_hat, xmat(wb_k*dt_gyro), Cbg_hat * OMb_k, Cbg_hat*Abg_hat*dt_gyro; ...
        let mut fk = Matrix18::identity();
        fk.index_mut((0..3, 0..3))
            .copy_from(self.psi_hat_mat.matrix());
        fk.index_mut((0..3, 3..6)).copy_from(&psi_cross);
        fk.index_mut((0..3, 6..15)).copy_from(&temp_mat_3);
        fk.index_mut((0..3, 15..18)).copy_from(&temp_mat_4);

        self.fk = fk;

        // Build process noise input matrix Lk
        let mut lk = Matrix18x9::zeros();
        lk.index_mut((0..3, 0..3)).copy_from(&temp_mat_4);
        lk.index_mut((3..6, 3..6))
            .copy_from(&na::Matrix3::identity());
        lk.index_mut((15..18, 6..9))
            .copy_from(&na::Matrix3::<f64>::identity());

        self.lk = lk;

        // Update the estimated covariance
        self.p_hat = fk * self.p_hat * fk.transpose() + lk * self._q_k * lk.transpose();
    }

    /// Function to calculate the expanded matrix of measured angular
    /// velocity + bias for use in the scaling factor matrix estimation
    ///
    /// # Detailed Explanation
    ///
    /// This function calculates the expanded matrix of measured angular
    /// velocity + bias for use in the scaling factor matrix estimation according to
    /// om_b_expmat = [omb_k' zeros(1,6); zeros(1,3) omb_k' zeros(1,3); zeros(1,6)
    ///  omb_k']
    ///
    /// # Arguments
    ///
    /// `self.t1: f64` - Time stamp
    /// `self.t0: f64` - Previous time stamp
    /// `self.om_b: na::Vector3<f64>` - Measured angular velocity + bias
    ///
    /// # Results
    ///
    /// `om_b_expmat: na::Matrix3x9` - Expanded matrix of measured angular velocity + bias
    pub fn calculate_om_b_expmat(&mut self, om_b: &na::Vector3<f64>) -> Matrix3x9 {
        // calculate the expanded matrix of measured angular velocity + bias
        let om_b_expmat = Matrix3x9::from_row_slice(&[
            om_b[0], om_b[1], om_b[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, om_b[0],
            om_b[1], om_b[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, om_b[0], om_b[1],
            om_b[2],
        ]);

        // return the expanded matrix of measured angular velocity + bias
        om_b_expmat
    }
}
