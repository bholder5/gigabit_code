//! This module assembles a struct of structs,
//!  containing one for each type of measurement 
//! generation necessary, ex gyros, star cameras etc
pub mod gyros;
use crate::initialization::Params;
use adcs::control::state as st;
use crate::flex_sim;
use nalgebra as na;
#[derive(Clone)]
pub struct Meas {
    pub gyros_bs: gyros::Gyro_bs,
    pub gyro_bow: gyros::Gyro_bs,
    pub gyro_stern: gyros::Gyro_bs,
    pub gyro_port: gyros::Gyro_bs,
    pub gyro_sb: gyros::Gyro_bs,
    pub gyro_of: gyros::Gyro_bs,
    pub gps: st::gps::Gps,
    pub cbh: na::Rotation3<f64>,
    pub c8h: na::Rotation3<f64>,
    pub c7h: na::Rotation3<f64>,
    pub roll: f64,
    pub pitch: f64,
    /// yaw_p is azimuth with opposite sign, not the position of OF WRT pivot
    pub yaw_p: f64,
}

impl Meas {
    pub fn new() -> Meas {
        Meas {
            gyros_bs: gyros::Gyro_bs::new(),
            gyro_bow: gyros::Gyro_bs::new_def(),
            gyro_stern: gyros::Gyro_bs::new_def(),
            gyro_port: gyros::Gyro_bs::new_def(),
            gyro_sb: gyros::Gyro_bs::new_def(),
            gyro_of: gyros::Gyro_bs::new_def(),
            gps: st::gps::Gps::new(),
            cbh: na::Rotation3::<f64>::identity(),
            c8h: na::Rotation3::<f64>::identity(),
            c7h: na::Rotation3::<f64>::identity(),
            roll: 0.0,
            pitch: 0.0,
            yaw_p: 0.0,
        }
    }

    pub fn read_measurements(&mut self, bp: &Params, sim_state: &st::State, flex: &flex_sim::Flex_model) -> () {
        // READ GYROS/GENERATE GYRO BORESIGHT MEASUREMENT
        // bp.omega_m was the "measurement" before this function, 
        // hence it is actually omega_k
        self.gps.get_greenwhich_apparent_sidereal_time();

        let last = self.gps.gast + self.gps.lon;
        let arg_z = (last + 180.0).to_radians();

        let arg_y = (self.gps.lat - 90.0).to_radians();

        let roty = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::y_axis(), arg_y);
        let rotz = na::Rotation3::<f64>::from_axis_angle(&na::Vector3::z_axis(), arg_z);
        let ceh = rotz * roty;

        let w_di = na::Vector3::<f64>::new(0.0, 0.0, 15.04108/206265.0);
        let w_di_h = ceh.transpose() * w_di;
        // need to add in flexible affects.
        
        self.gyros_bs.omega_k = bp.omega_m + (self.cbh * w_di_h);

        // ENCODERS

        if flex.flex_enable{
            self.roll = bp.x[16] + flex.c_pos_out[1];
            self.pitch = bp.x[17] + flex.c_pos_out[3];
        } else {
            self.roll = bp.x[16];
            self.pitch = bp.x[17];
        }

        // frame based gyros
        // self.gyro_bow.omega_k = bp.omega_m_roll + (self.c8h * ceh.transpose() * w_di);
        // self.gyro_stern.omega_k = bp.omega_m_roll + (self.c8h * ceh.transpose() * w_di);
        // self.gyro_of.omega_k = bp.omega_m_yaw + (self.c7h * ceh.transpose() * w_di);
        // self.gyro_port.omega_k = self.gyros_bs.omega_k.clone();
        // self.gyro_sb.omega_k = self.gyros_bs.omega_k.clone();
        let omega_m_roll = bp.omega_m_roll.clone() + self.c8h * w_di_h;
        let omega_m_yaw = bp.omega_m_yaw.clone() + self.c7h * w_di_h;
        let omega_m_pitch = bp.omega_m.clone() + self.cbh * w_di_h;

        self.gyro_bow.omega_k = omega_m_roll;
        self.gyro_stern.omega_k = omega_m_roll;
        self.gyro_of.omega_k = omega_m_yaw;
        self.gyro_port.omega_k = omega_m_pitch;
        self.gyro_sb.omega_k = omega_m_pitch;
        

        // add flexible affects
        if flex.flex_enable{
            self.gyros_bs.omega_k = flex.g1_out + self.gyros_bs.omega_k;
        }

        // //generate measurements
        // self.gyros_bs.generate_measurement();
        // self.gyro_bow.generate_measurement();
        // self.gyro_stern.generate_measurement();
        // self.gyro_port.generate_measurement();
        // self.gyro_sb.generate_measurement();
        // self.gyro_of.generate_measurement();

        // axial measurements
        self.gyro_bow.omega_k[0] = self.gyro_bow.omega_k[0] + flex.c_out[1];
        self.gyro_stern.omega_k[0] = self.gyro_stern.omega_k[0] + flex.c_out[2];
        self.gyro_port.omega_k[1]  = self.gyro_port.omega_k[1] + flex.c_out[3];
        self.gyro_sb.omega_k[1] = self.gyro_sb.omega_k[1] + flex.c_out[4];
        self.gyro_of.omega_k[2] = self.gyro_of.omega_k[2] + flex.c_out[0];

        // println!("{}", self.gyro_bow.omega_m);

        self.gyros_bs.generate_measurement();
        self.gyro_bow.generate_measurement();
        self.gyro_stern.generate_measurement();
        self.gyro_port.generate_measurement();
        self.gyro_sb.generate_measurement();
        self.gyro_of.generate_measurement();

        // println!("{}", self.gyro_bow.omega_m);

        // axial measurements
        self.gyro_bow.om_axial = self.gyro_bow.omega_m[0];
        self.gyro_stern.om_axial = self.gyro_stern.omega_m[0];
        self.gyro_port.om_axial = self.gyro_port.omega_m[1];
        self.gyro_sb.om_axial = self.gyro_sb.omega_m[1];
        self.gyro_of.om_axial = self.gyro_of.omega_m[2];

        // println!("{} {} {} {} {}", self.gyro_bow.om_axial,self.gyro_stern.om_axial,self.gyro_port.om_axial,self.gyro_sb.om_axial,self.gyro_of.om_axial);

        
        
        
        //Cant tell if I should add the flex out because this is calculated from boresight which may already contain it?
        // BUT add flex into boresight gyros and thus boresight estimate. Not correction though?
        self.yaw_p = -sim_state.hor.az + flex.c_pos_out[0];

        // adjust measurements for 

    }
}