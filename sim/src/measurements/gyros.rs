use nalgebra as na;
use rand::prelude::*;
use rand_distr::{Normal, Distribution};

#[derive(Clone)]
/// struct for the boresight gyroscope measurements generation
pub struct Gyro_bs {
    /// Actual boresight angular velocity
    pub omega_k: na::Vector3<f64>,
    /// Measured angular velocity (after offset and scale)
    pub omega_m: na::Vector3<f64>,
    /// sliding window filter omega
    pub omega_m_wnd: na::Vector3<f64>,
    /// Gyroscope offset (misalignment in mounting)
    pub c_bg_i: na::Rotation3<f64>,
    /// Gyroscope scale (gain and internal alignment offsets)
    pub a_g_inv: na::Matrix3<f64>,
    /// Gyroscope bias
    pub bias: na::Vector3<f64>,
    /// maximum bias drift magnitude
    pub bias_drift: f64,
    pub bias_dist: Normal<f64>,
    /// axial reading for flexible measurements (the axis being read)
    pub om_axial: f64,
    /// Standard deviation of the noise
    pub noise_std_dev: f64,
    pub noise_dist: Normal<f64>,
}

impl Gyro_bs {
    pub fn new() -> Gyro_bs {

        let dt: f64 = 2.0 / 1000.0;

        // // KVH Characteristics
        // let arw: f64 = 0.013;//deg/sqrt(hr)
        // let bi: f64 = 0.05;

        // Emcore 
        let arw: f64 = 0.002;//deg/sqrt(hr)
        let bi: f64 = 0.02; //deg/hour

        //IxBlue
        // let arw: f64 = 1.0*0.002;
        // let bi: f64 = 1.0*0.0065;



        // Convert ARW from deg/sqrt(hr) to rad/sqrt(Hz)
        let arw_noise = arw.to_radians() / 60.0;

        // Convert BI from deg/hr to deg/s
        let bi_noise = bi.to_radians() / 3600.0;

        let noise_dist = Normal::new(0.0, arw_noise / dt.sqrt()).unwrap();
        let bias_dist = Normal::new(0.0, bi_noise * dt.sqrt()).unwrap();

        let new_gyro = Gyro_bs {
            omega_k: na::Vector3::new(0.0, 0.0, 0.0),
            omega_m: na::Vector3::new(0.0, 0.0, 0.0),
            omega_m_wnd: na::Vector3::new(0.0,0.0,0.0),
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
            bias: na::Vector3::new(0.0001, -0.0005, 0.00075),
            bias_drift: 0.0,
            bias_dist,
            om_axial: 0.0,
            noise_std_dev: 0.0000,  // Example standard deviation for noise
            noise_dist
            
        };

        return new_gyro
    }
    
    //a default gyro with no bias and stuff
    pub fn new_def() -> Gyro_bs {

        let dt: f64 = 2.0 / 1000.0;

        // // KVH Characteristics
        // let arw: f64 = 0.013;//deg/sqrt(hr)
        // let bi: f64 = 0.05;

        // Emcore 
        let arw: f64 = 0.0;//deg/sqrt(hr)
        let bi: f64 = 0.0; //deg/hour

        //IxBlue
        // let arw: f64 = 1.0*0.002;
        // let bi: f64 = 1.0*0.0065;



        // Convert ARW from deg/sqrt(hr) to rad/sqrt(Hz)
        let arw_noise = arw.to_radians() / 60.0;

        // Convert BI from deg/hr to deg/s
        let bi_noise = bi.to_radians() / 3600.0;

        let noise_dist = Normal::new(0.0, arw_noise / dt.sqrt()).unwrap();
        let bias_dist = Normal::new(0.0, bi_noise * dt.sqrt()).unwrap();

        Gyro_bs {
            omega_k: na::Vector3::new(0.0, 0.0, 0.0),
            omega_m: na::Vector3::new(0.0, 0.0, 0.0),
            omega_m_wnd: na::Vector3::new(0.0,0.0,0.0),
            // remember to invert
            c_bg_i: na::Rotation3::identity(),
            // remember to invert
            a_g_inv: na::Matrix3::identity(),
            // a_g_inv: na::Matrix3::<f64>::identity(),
            bias: na::Vector3::new(0.0, 0.0, 0.0),
            bias_drift: 0.0,
            bias_dist, 
            om_axial: 0.0,
            noise_std_dev: 0.00000,  // No noise
            noise_dist
        }
    }

    /// Generate a realistic gyroscope measurement
    /// 
    /// # Detailed Explanation
    /// 
    /// This function generates a realistic gyroscope measurement 
    /// following the below equation
    /// 
    /// omega_m = (a_g_inv * c_bg * omega_k) - bias + noise
    /// 
    /// where omega_k is the actual angular velocity, c_bg is the 
    /// gyroscope misalignment, a_g is the internal misalignment 
    /// and scaling factor, bias is the gyroscope bias, and noise is 
    /// Gaussian noise.
    /// 
    /// # Arguments
    /// 
    /// `self.omega_k: na::Vector3<f64>` - Actual angular velocity
    /// `self.c_bg: na::Rotation3<f64>` - Gyroscope misalignment
    /// `self.a_g: na::Matrix3<f64>` - Gyroscope scaling factor
    /// `self.bias: na::Vector3<f64>` - Gyroscope bias
    /// `self.noise_std_dev: f64` - Standard deviation of the noise
    /// 
    /// # Results
    /// 
    /// `self.omega_m: na::Vector3<f64>` - Measured angular velocity
    pub fn generate_measurement(&mut self) {
        
        let mut rng = thread_rng();

        let bias_drift = na::Vector3::new(self.bias_dist.sample(&mut rng),self.bias_dist.sample(&mut rng),self.bias_dist.sample(&mut rng));

        self.bias = self.bias + 1.0*bias_drift;


        // generate a realistic gyroscope measurement
        self.omega_m = (self.a_g_inv * self.c_bg_i * self.omega_k) - self.bias;
        // self.omega_m = self.omega_k;

        
        // println!("Meas: Gyro before {}", self.omega_k);
        // Add Gaussian noise to each component of omega_m
        self.omega_m.x += (self.noise_dist.sample(&mut rng));
        self.omega_m.y += (self.noise_dist.sample(&mut rng));
        self.omega_m.z += (self.noise_dist.sample(&mut rng));
        // println!("Meas: Gyro After {}", self.omega_m);
    }   
}