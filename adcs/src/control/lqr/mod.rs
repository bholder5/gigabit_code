extern crate nalgebra as na;
use crate::control::state as st;
use crate::control::error as er;

pub mod amat;
pub mod bmat;
pub mod pos;
pub mod neg;

pub type Vector8 = na::SMatrix<f64, 8, 1>;
pub type Matrix48 = na::SMatrix<f64, 4, 8>;
pub type Matrix8 = na::SMatrix<f64, 8, 8>;
pub type Matrix1 = na::SMatrix<f64, 1, 1>;

#[derive(Debug)]
pub struct LQRController {
    // LQR parameters
    lqr_gain_p: Matrix48, // LQR gain matrix
    lqr_gain_n: Matrix48, // LQR gain matrix
    pub control_input: na::Vector4<f64>, // Control input
    q: Matrix8, // State weighting matrix
    r: na::Matrix4<f64>, // Control weighting matrix
    ap: amat::AMat,
    bp: bmat::BMat,
    an: amat::AMat,
    bn: bmat::BMat,
    err: Vector8,
    roll: f64,
    roll_torque: f64,
    pitch: f64,
    pitch_k: f64,
}

impl LQRController {
    // Initialize the LQR controller with default values
    pub fn init() -> LQRController {
        let ap = amat::init_amat_pos();
        // println!("a {}", &a);
        let bp = bmat::init_bmat_pos();
        // println!("b {}", &b);
        let an = amat::init_amat_neg();
        // println!("a {}", &a);
        let bn = bmat::init_bmat_neg();

        let q = 1.0*Matrix8::from_diagonal(&Vector8::from_row_slice(
            &[
            100.0, //piv angl
            1950.0, // yaw angle
            1950.0, // roll angle
            150.0, // pitch angle
            75000.0, //yaw speed
            100000.0, // roll speed
            10000.0, // pitch speed
            0.00000001,
            ]) //momentum RW
        );

        let r =  1.0*na::Matrix4::from_diagonal(&na::Vector4::from_row_slice(
            &[
            0.02, //RW
            0.1, // Roll
            0.1, // Pitch
            100000000000.0 // piv
        ]));

        // let q = 80.2*Matrix10::from_diagonal(&Vector10::from_row_slice(
        //     &[
        //     0.0000001, //piv angle
        //     206265.0, // yaw angle
        //     406265.0, // roll angle
        //     206265.0, // pitch angle
        //     1000.0, //yaw speed
        //     1000.0, // roll speed
        //     1000.0, // pitch speed
        //     0.0000001,
        //     1000000.0,
        //     100000.0]) //momentum RW
        // );

        // let r =  100000.0*na::Matrix4::from_diagonal(&na::Vector4::from_row_slice(
        //     &[
        //     0.5, //RW
        //     2.0, // Roll
        //     2.0, // Pitch
        //     1000000.0 // piv
        // ]));
        let control_input = na::Vector4::<f64>::zeros();
        let lqr_gain_p = Matrix48::zeros();
        let lqr_gain_n = Matrix48::zeros();
        let err = Vector8::zeros();
        let roll = 0.0;
        let roll_torque = 6.260767;
        let pitch = 0.0;
        let pitch_k = 62.60767;
        let mut lqr = LQRController {
            lqr_gain_p, // Initialized to zeros, should be computed using a Riccati solver
            lqr_gain_n,
            control_input,
            q,
            r,
            ap,
            bp,
            an,
            bn,
            err,
            roll,
            roll_torque,
            pitch,
            pitch_k
        };

        lqr.care();
        // println!("gain {:.4}", &lqr.lqr_gain);
        return lqr
    }

    // Calculate the control input using LQR
    pub fn calculate_control(&mut self) {
        // Calculate the control input using the formula: u = -K * (x - x_ref)
        let vec_pos = na::Vector4::<f64>::from_row_slice(&[0.0,self.roll_torque,0.0, 0.0]);
        let vec_pitch = na::Vector4::<f64>::from_row_slice(&[0.0,0.0,self.pitch_k * (self.pitch+0.698132), 0.0]);
        let cont_pos = self.lqr_gain_p * (self.err) + vec_pos;
        let cont_neg = self.lqr_gain_n * self.err - vec_pos;
        
        let weight = (self.roll + 0.1)/0.2;
        let control_in = (weight * cont_pos) + ((1.0-weight)*cont_neg);
        self.control_input = control_in + vec_pitch;
    }

    // Get the control input
    pub fn get_control_input(&self) -> na::Vector4<f64> {
        self.control_input.clone()
    }

    /// assemble the measurement vector for the control calculation
    pub fn assem_err_vec(&mut self, state: &st::State, error: &er::Error){
        // pivot speed pos
        let mut err = Vector8::zeros();
        err.slice_mut((0,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[-state.piv_angle]));

        // yaw angle
        let yaw_err = error.err_comb_th.z + 0.0*state.piv_angle;
        err.slice_mut((1,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[yaw_err]));

        // roll angle
        let roll_err = error.err_comb_th.x;
        err.slice_mut((2,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[roll_err]));

        // pitch angle
        let pitch_err = error.err_comb_th.y;
        err.slice_mut((3,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[pitch_err]));

        // yaw anglular rate
        let yaw_rate_err = error.err_rate_lqr.z + 0.0*state.piv_speed;
        err.slice_mut((4,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[yaw_rate_err]));

        // roll anglular rate
        let roll_rate_err = error.err_rate_lqr.x;
        err.slice_mut((5,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[roll_rate_err]));

        // pitch anglular rate
        let pitch_rate_err = error.err_rate_lqr.y;
        err.slice_mut((6,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[pitch_rate_err]));

        // rw hs error
        let hs_err = state.rw_hs_nom - state.rw_hs;
        err.slice_mut((7,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[hs_err]));

        // flexure addition
        // err.slice_mut((8,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[error.err_fine_sum[0]]));
        // // flexure addition
        // err.slice_mut((9,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[error.err_fine_sum[2]]));
        // err.slice_mut((10,0), (1,1)).copy_from(&Matrix1::from_row_slice(&[error.err_fine_sum[1]]));

        self.err = err.clone();

        // update roll position
        self.roll = state.gmb_k.roll;
        self.pitch = state.gmb_k.pitch;
    }

    pub fn care(&mut self) {
        //
        // 
        //  POSITIVE
        // 
        let a = self.ap.clone();
        let b = self.bp.clone();
        let r = self.r.clone();
        let q = self.q.clone();

        let num_states = 8;
    
        let z12 = -&b * r.try_inverse().unwrap() * &b.transpose();
        let mut z = na::DMatrix::<f64>::zeros({2*num_states}, {2*(num_states)});

        z.slice_mut((0,0), (num_states,num_states)).copy_from(&a);
        z.slice_mut((0,num_states), (num_states,num_states)).copy_from(&z12);
        z.slice_mut((num_states,0), (num_states,num_states)).copy_from(&-q);
        z.slice_mut((num_states,num_states), (num_states,num_states)).copy_from(&-a.transpose());
    
        let decom = z.clone().schur();
    
        let eigen_vals = decom.complex_eigenvalues().clone();
    
        let mut z_c = na::DMatrix::<na::Complex<f64>>::zeros((2*num_states), (2*num_states));
    
        for (el_n, el) in z.clone().iter().enumerate() {
            z_c[el_n] = na::Complex::<f64>::new(el.clone(), 0.0);
        }
    
    
        // println!("{:.3} {:.3} {}", z, z_c, eigen_vals);
        let mut cnt = 0;
        let mut u11 = na::DMatrix::<na::Complex<f64>>::zeros({ num_states }, { num_states });
        let mut u12 = na::DMatrix::<na::Complex<f64>>::zeros({ num_states }, { num_states });
        
        for value in eigen_vals.iter(){
            // println!("{}", value);
            if value.re < 0.0{
    
                let mut temp = na::DMatrix::<na::Complex<f64>>::identity({ 2*num_states } , { 2*num_states });
    
    
                //Create the lamba * identity matrix
    
                for (loc, el) in temp.clone().iter().enumerate(){
                    temp[loc] = temp[loc] * value;
                }
    
                let temp2 = &z_c - temp;
    
                let svd_d = temp2.svd(true, true);
    
                let v = svd_d.v_t.unwrap().transpose();
                let mut eig_vn = na::DMatrix::<na::Complex<f64>>::zeros({ 2*num_states } , 1);
    
                for (col_cnt, col) in v.column_iter().enumerate(){
                    if col_cnt == ((2*num_states)-1){
                        eig_vn = na::DMatrix::<na::Complex<f64>>::from_column_slice({2*num_states}, 1, &col.as_slice());
                    }
                }
                u11.slice_mut((0,cnt), ({num_states},1)).copy_from(&eig_vn.slice((0, 0), (num_states,1)));
                u12.slice_mut((0,cnt), ({num_states},1)).copy_from(&eig_vn.slice((num_states, 0), (num_states,1)));
                
                cnt += 1;
    
            }
    
        }
    
        let Pc = u12 * u11.try_inverse().unwrap();
    
        let mut p = na::DMatrix::<f64>::zeros(num_states, num_states);
    
        let mut cn = 0;
        for el in Pc.iter(){
            p[cn] = el.re;
            cn +=1;
        }
    
        // println!("THIS IS P {:+.5e}", p);
        let k = r.try_inverse().unwrap() * b.transpose() * p;
        
        for (loc, el) in k.iter().enumerate(){
            self.lqr_gain_p[loc] = el.clone();
        }   
        //
        // 
        // 
        // 
        // NEGATIVE
        // 
        // 
        // 
        //
        let a = self.an.clone();
        let b = self.bn.clone();
        let r = self.r.clone();
        let q = self.q.clone();
    
        let z12 = -&b * r.try_inverse().unwrap() * &b.transpose();
        let mut z = na::DMatrix::<f64>::zeros({2*num_states}, {2*(num_states)});

        z.slice_mut((0,0), (num_states,num_states)).copy_from(&a);
        z.slice_mut((0,num_states), (num_states,num_states)).copy_from(&z12);
        z.slice_mut((num_states,0), (num_states,num_states)).copy_from(&-q);
        z.slice_mut((num_states,num_states), (num_states,num_states)).copy_from(&-a.transpose());
    
        let decom = z.clone().schur();
    
        let eigen_vals = decom.complex_eigenvalues().clone();
    
        let mut z_c = na::DMatrix::<na::Complex<f64>>::zeros((2*num_states), (2*num_states));
    
        for (el_n, el) in z.clone().iter().enumerate() {
            z_c[el_n] = na::Complex::<f64>::new(el.clone(), 0.0);
        }
    
    
        // println!("{:.3} {:.3} {}", z, z_c, eigen_vals);
        let mut cnt = 0;
        let mut u11 = na::DMatrix::<na::Complex<f64>>::zeros({ num_states }, { num_states });
        let mut u12 = na::DMatrix::<na::Complex<f64>>::zeros({ num_states }, { num_states });
        
        for value in eigen_vals.iter(){
            // println!("{}", value);
            if value.re < 0.0{
    
                let mut temp = na::DMatrix::<na::Complex<f64>>::identity({ 2*num_states } , { 2*num_states });
    
    
                //Create the lamba * identity matrix
    
                for (loc, el) in temp.clone().iter().enumerate(){
                    temp[loc] = temp[loc] * value;
                }
    
                let temp2 = &z_c - temp;
    
                let svd_d = temp2.svd(true, true);
    
                let v = svd_d.v_t.unwrap().transpose();
                let mut eig_vn = na::DMatrix::<na::Complex<f64>>::zeros({ 2*num_states } , 1);
    
                for (col_cnt, col) in v.column_iter().enumerate(){
                    if col_cnt == ((2*num_states)-1){
                        eig_vn = na::DMatrix::<na::Complex<f64>>::from_column_slice({2*num_states}, 1, &col.as_slice());
                    }
                }
                u11.slice_mut((0,cnt), ({num_states},1)).copy_from(&eig_vn.slice((0, 0), (num_states,1)));
                u12.slice_mut((0,cnt), ({num_states},1)).copy_from(&eig_vn.slice((num_states, 0), (num_states,1)));
                
                cnt += 1;
    
            }
    
        }
    
        let Pc = u12 * u11.try_inverse().unwrap();
    
        let mut p = na::DMatrix::<f64>::zeros(num_states, num_states);
    
        let mut cn = 0;
        for el in Pc.iter(){
            p[cn] = el.re;
            cn +=1;
        }
    
        // println!("THIS IS P {:+.5e}", p);
        let k = r.try_inverse().unwrap() * b.transpose() * p;
        
        for (loc, el) in k.iter().enumerate(){
            self.lqr_gain_n[loc] = el.clone();
        }        
    }
}
