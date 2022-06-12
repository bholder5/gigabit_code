//! Sub module for generating realistic gyroscop measurements
//! 
//! This includes injecting offsets, changes in scale and bias

use nalgebra as na;
#[derive(Clone)]
/// struct for the boresight gyroscope measurements generation
pub struct Gyro_bs {
    /// Actual boresight angular velocity
    pub omega_k: na::Vector3<f64>,
    /// Measured angular velocity (after offset and scale)
    pub omega_m: na::Vector3<f64>,
    /// Gyroscope offset (misalignment in mounting)
    pub c_bg_i: na::Rotation3<f64>,
    /// Gyroscope scale (gain and internal alignment offsets)
    pub a_g_inv: na::Matrix3<f64>,
    /// Gyroscope bias
    pub bias: na::Vector3<f64>,
    /// maximum bias drift magnitude
    pub bias_drift: f64,
}

impl Gyro_bs {
    pub fn new() -> Gyro_bs {
        Gyro_bs {
            omega_k: na::Vector3::new(0.0, 0.0, 0.0),
            omega_m: na::Vector3::new(0.0, 0.0, 0.0),
            // remember to invert
            c_bg_i: na::Rotation3::identity(),
            // remember to invert
            a_g_inv: na::Matrix3::from_row_slice(&[1.004999976946721,
                0.000096039954031,
                0.000192049170358,
               -0.000096015363867,
                0.994999987192623,
               -0.000128046099911,
               -0.000192061465440,
                0.000128027657288,
                1.002499973360656]).try_inverse().unwrap(),
            // a_g_inv: na::Matrix3::<f64>::identity(),
            bias: na::Vector3::new(0.000001, -0.000005, 0.0000075),
            bias_drift: 0.0,
        }
    }

    /// Generate a realistic gyroscope measurement
    /// 
    /// # Detailed Explanation
    /// 
    /// This function generates a realistic gyroscope measurement 
    /// following the below equation
    /// 
    /// omega_m = (a_g_inv * c_bg * omega_k) - bias
    /// 
    /// where omega_k is the actual angular velocity, c_bg is the 
    /// gyroscope misalignment, a_g is the internal misalignment 
    /// and scaling factor and bias is the gyroscope bias.
    /// 
    /// # Arguments
    /// 
    /// `self.omega_k: na::Vector3<f64>` - Actual angular velocity
    /// `self.c_bg: na::Rotation3<f64>` - Gyroscope misalignment
    /// `self.a_g: na::Matrix3<f64>` - Gyroscope scaling factor
    /// `self.bias: na::Vector3<f64>` - Gyroscope bias
    /// 
    /// # Results
    /// 
    /// `self.omega_m: na::Vector3<f64>` - Measured angular velocity
    pub fn generate_measurement(&mut self) {
        // generate a realistic gyroscope measurement
        self.omega_m = (self.a_g_inv * self.c_bg_i * self.omega_k) - self.bias;
                
    }   
}