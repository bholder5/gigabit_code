extern crate nalgebra as na;
use crate::control::error as er;

pub mod servo;

pub type rwVector = na::SMatrix<f64, 1, 6>;
pub type Vector6 = na::SMatrix<f64, 6,1>;
pub type Matrix4x2 = na::SMatrix<f64,4,2>;

#[derive(Debug)]
pub struct ServoAxis {
    p: na::Matrix4<f64>,
    r: Matrix4x2,
    k: rwVector,
    dt: f64,
    rate: f64,
    tau_req: f64,
    tau_lim: f64,
    pub z: na::Vector4<f64>,
    max: f64,
}

impl ServoAxis {
    pub fn new(p: &na::Matrix4<f64>, r: &Matrix4x2, k: &rwVector, dt: &f64, rate: &f64, max: &f64) -> ServoAxis {
        let tau_req = 0.0;
        let tau_lim = 0.0;
        let z = na::Vector4::<f64>::zeros();
        let servoAxis = ServoAxis {
            p: p.clone(),
            r: r.clone(),
            k: k.clone(),
            dt: dt.clone(),
            rate: rate.clone(),
            tau_req,
            tau_lim,
            z,
            max: max.clone(),
        };
        return servoAxis
    }

    pub fn propogate(&mut self, err_pos: &f64, err_rate: &f64){
        let err = na::Vector2::<f64>::from_row_slice(&[err_pos.clone(), err_rate.clone()]);
        let zdot = (self.p*self.z + self.r*err)*self.dt;
        self.z = self.z + zdot;

    }
    pub fn get_control_input(&mut self, err_pos: &f64, err_rate: &f64) -> f64{
        let state_vec = [err_pos.clone(), err_rate.clone(), self.z[0], self.z[1], self.z[2], self.z[3]];
        let state = Vector6::from_row_slice(&state_vec);
        let mut tau_req = self.k * state/1.0;

        // println!("{}, {}", self.z[1], self.z[2]);
        
        if tau_req.norm() > self.max{
            let tau_bound = tau_req[0].signum() * self.max;
            tau_req = na::Matrix1::new(tau_bound);
        }
        
        self.tau_req = tau_req[0];

        return self.tau_req
    }

}

pub struct ServoControl {
    pub yaw: ServoAxis,
    pub roll: ServoAxis,
    pub pitch: ServoAxis,
    // dyanmic
    pub tau_y_lim: f64,
    pub tau_r_lim: f64,
    pub tau_p_lim: f64,
    pub enable: bool,
}

impl ServoControl {
    pub fn init() -> ServoControl{
        let p = servo::init_AMat1_pos();
        let r = servo::init_AMat2_pos();
        let k_yaw = servo::init_AMat3_pos();
        let k_roll = servo::init_AMat4_pos();
        let k_pitch = servo::init_AMat5_pos();
        let rate = 0.0;
        let dt = 0.001;
        let max = 10.0;

        let yaw = ServoAxis::new(&p, &r, &k_yaw, &dt, &rate, &max);
        let roll = ServoAxis::new(&p, &r, &k_roll, &dt, &rate, &max);
        let pitch = ServoAxis::new(&p, &r, &k_pitch, &dt, &rate, &max);
        let tau_y_lim = 0.0;
        let tau_r_lim = 0.0;
        let tau_p_lim = 0.0;
        let enable = true;
        

        let servoCont = ServoControl{
            yaw,
            roll,
            pitch,
            tau_y_lim,
            tau_r_lim,
            tau_p_lim,
            enable,
        };

        return servoCont
    }

    pub fn propogate(&mut self, err: &na::Vector3::<f64>, err_rate: &na::Vector3::<f64>){
        self.yaw.propogate(&err[2], &err_rate[2]);
        self.roll.propogate(&err[0], &err_rate[0]);
        self.pitch.propogate(&err[1], &err_rate[1]);
    }

    pub fn get_control_input(&mut self, err: &na::Vector3::<f64>, err_rate: &na::Vector3::<f64>){
        let yaw_req = self.yaw.get_control_input(&err[2], &err_rate[2]);
        let roll_req = self.roll.get_control_input(&err[0], &err_rate[0]);
        let pitch_req = self.pitch.get_control_input(&err[1], &err_rate[1]);

        // rate limit here

        //

        self.tau_y_lim = yaw_req/1.0;
        self.tau_r_lim = roll_req/1.0;
        self.tau_p_lim = pitch_req/1.0;
    }
}