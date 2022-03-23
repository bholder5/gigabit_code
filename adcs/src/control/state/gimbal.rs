// import crate for matrix math
extern crate nalgebra as na;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};

#[derive(Debug, Clone)]
pub struct Gimbal {
    pub roll: f64,
    pub pitch: f64,
    pub yaw: f64,
    pub rot: na::Rotation3<f64>,
    pub gmm: na::Matrix3<f64>,
}

impl Gimbal {
    pub fn calculate_rotation_matrix(&mut self){
        let rotx = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::x_axis(), self.roll).inverse();
        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), self.pitch).inverse();
        let rotz = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::z_axis(), self.yaw).inverse();

        // calculating the CbH using the gimbal euler 312 rotation.
        self.rot = roty * rotx * rotz;
        // println!("gimbal rot: {}", self.rot);
    }

    pub fn update_gimbal_coordinates(&mut self, vec: &[f64; 3]){
        trace!("update_gimbal_coordinates start");
        self.roll  = vec[0];
        self.pitch = vec[1];
        self.yaw  = vec[2];

        self.calculate_rotation_matrix();
        trace!("update_gimbal_coordinates end");
    }

    pub fn extract_gimbal_rpy(&mut self){
        trace!("extract_gimbal_rpy start");
        self.roll = self.rot[(1,2)].asin();
        self.pitch = (-self.rot[(0,2)]).atan2(self.rot[(2,2)]);
        self.yaw = (-self.rot[(1,0)]).atan2(self.rot[(1,1)]);
        // println!("roll: {}, pitch {}, yaw {}", self.roll, self.pitch, self.yaw);
        trace!("extract_gimbal_rpy end");
    }

    pub fn calculate_gimbal_mapping_matrix(&mut self){
        trace!("calculate_gimbal_mapping_matrix start");
        let rotx = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::x_axis(), self.roll).inverse();
        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), self.pitch).inverse();

        let col1 = roty * na::Vector3::<f64>::new(1.0, 0.0, 0.0);
        let col2 = na::Vector3::<f64>::new(0.0, 1.0, 0.0);
        let col3 = roty * rotx * na::Vector3::<f64>::new(0.0, 0.0, 1.0);
        let gmm = na::Matrix3::<f64>::from_columns(&[col1, col2, col3]);
        self.gmm = gmm;
        trace!("calculate_gimbal_mapping_matrix end");
    }

    pub fn new() -> Gimbal {
        let roll: f64 = -0.0259154632;
        let pitch: f64 = 0.750145461;
        let yaw: f64 = -0.0278;
        let rot = na::Rotation3::<f64>::identity();
        let gmm = na::Matrix3::<f64>::identity();
        let mut gimbal = Gimbal{
            roll,
            pitch,
            yaw,
            rot,
            gmm,
        };
        gimbal.calculate_rotation_matrix();
        gimbal.extract_gimbal_rpy();
        gimbal
    }
}