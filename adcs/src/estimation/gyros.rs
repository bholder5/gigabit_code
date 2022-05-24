//! Sub module for generating realistic gyroscop measurements
//!
//! This includes injecting offsets, changes in scale and bias

use nalgebra as na;

/// struct for the boresight gyroscope measurements generation
#[derive(Debug)]
pub struct GyroBs {
    /// Actual boresight angular velocity
    pub omega_k: na::Vector3<f64>,
    /// Measured angular velocity (after offset and scale)
    pub omega_m: na::Vector3<f64>,
    /// Gyroscope offset (misalignment in mounting)
    pub c_bg: na::Rotation3<f64>,
    /// Gyroscope scale (gain and internal alignment offsets)
    pub a_g: na::Matrix3<f64>,
    /// Gyroscope bias
    pub bias: na::Vector3<f64>,
    /// Previous measurement time stamp
    pub t0: f64,
    /// Latest measurement time stamp
    pub t1: f64,
    /// measured angular velocity + bias for calibration estimation
    pub om_b: na::Vector3<f64>,
}

impl GyroBs {
    pub fn new() -> GyroBs {
        GyroBs {
            omega_k: na::Vector3::new(0.0, 0.0, 0.0),
            omega_m: na::Vector3::new(0.0, 0.0, 0.0),
            c_bg: na::Rotation3::identity(),
            a_g: na::Matrix3::identity(),
            bias: na::Vector3::new(0.0, 0.0, 0.0),
            t0: 0.0,
            t1: 0.0,
            om_b: na::Vector3::new(0.0, 0.0, 0.0),
        }
    }

    /// Deconstruct the gyroscope measurement to get the actual angular velocity
    ///
    /// # Detailed Explanation
    ///
    /// This function deconstructs a gyroscope measurement
    /// following the below equation
    ///
    /// omega_k = c_bg * a_g * (omega_m + bias)
    ///
    /// where omega_k is the actual angular velocity, c_bg is the
    /// gyroscope misalignment, a_g is the internal misalignment
    /// and scaling factor and bias is the gyroscope bias.
    ///
    /// # Arguments
    ///
    /// `self.omega_m: na::Vector3<f64>` - Measured angular velocity
    /// `self.c_bg: na::Rotation3<f64>` - Gyroscope misalignment
    /// `self.a_g: na::Matrix3<f64>` - Gyroscope scaling factor
    /// `self.bias: na::Vector3<f64>` - Gyroscope bias
    ///
    /// # Results
    ///
    /// `self.omega_k: na::Vector3<f64>` - Actual angular velocity
    pub fn deconstruct_measurement(&mut self) {
        // deconstruct the gyroscope measurement
        self.omega_k = self.c_bg * self.a_g * (self.omega_m + self.bias);
    }

    /// Function to read in a new measurement
    ///
    /// # Detailed Explanation
    ///
    /// This function takes in a new gyro reading along with the time stamp.
    /// The old time stamp is shifted from t1 to t0 so that the new dt can be
    /// calculated for the gyro progagation. The gyro measurement is then
    /// deconstructed in preparation for the propogation
    ///
    /// # Arguments
    ///
    /// `omega_m: na::Vector3<f64>` - Measured angular velocity
    /// `t1: f64` - Time stamp
    ///
    /// # Results
    ///
    /// `omega_k: na::Vector3<f64>` - actual angular velocity
    pub fn read_gyros(&mut self, omega_m: na::Vector3<f64>, t1: f64) {
        // shift the time stamps
        self.t0 = self.t1;
        self.t1 = t1;

        // update the gyroscope measurement
        self.omega_m = omega_m;

        // deconstruct the gyroscope measurement
        self.deconstruct_measurement();
        self.om_b = self.omega_m + self.bias;
    }
    pub fn reset(&mut self){
        self.c_bg = na::Rotation3::identity();
        self.a_g = na::Matrix3::identity();
        self.bias = na::Vector3::new(0.0, 0.0, 0.0);

    }
}
