extern crate nalgebra as na;

mod amat;
mod bmat;

pub use amat::*;
pub use bmat::*;

type BRigid = na::SMatrix<f64, 3, 5>;

pub struct PassiveControl {
    num_modes: usize,
    num_states: usize,
    a_flex: na::SMatrix<f64, 104, 104>,
    b_flex: na::SMatrix<f64, 104, 5>,
    a_rigid: na::Matrix3<f64>,
    b_rigid: BRigid,
    qr: na::DMatrix<f64>,
    r: na::SMatrix<f64, 5, 5>,
    ql: na::DMatrix<f64>,
    a_cl: na::DMatrix<f64>,
    p_lyap: na::DMatrix<f64>,
    b_c: na::DMatrix<f64>,
    c_c: na::DMatrix<f64>,
    state: na::DVector<f64>,
    pub u: na::DVector<f64>, //the requested torque
    pub u_rigid: na::DVector<f64>,
    pub enable: bool,
}

impl PassiveControl{
    pub fn init_pc() -> PassiveControl{
        let a_flex = init_amat();
        let b_flex = 1.0 *init_bmat();
        let a_rigid = na::Matrix3::<f64>::from_row_slice(&[978.1915, -15.2645, 0.0, -15.2645, 376.8085, 0.0, 0.0, 0.0, 134.0]);
        let b_rigid = a_rigid.try_inverse().unwrap() * BRigid::from_row_slice(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]);
        let num_modes: usize = 8;
        let num_states: usize = 2*num_modes;
        let r = 0.0001 * na::Matrix5::<f64>::identity();
        let mut ql = 1.0*na::DMatrix::<f64>::identity(num_states+3, num_states+3);
        let mut qr = 1.0 * na::DMatrix::<f64>::identity(num_states+3, num_states+3);

        let qr_rigid = 0.0015* a_rigid.clone();
        qr.slice_mut((0,0), (3, 3)).copy_from(&qr_rigid);
        qr.slice_mut((3,3), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 10.0)));
        qr.slice_mut((5,5), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 100.0)));
        qr.slice_mut((7,7), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 1.0)));
        qr.slice_mut((9,9), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 2.0)));
        qr.slice_mut((11,11), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 1.0)));
        qr.slice_mut((13,13), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 1.0)));
        qr.slice_mut((15,15), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 1.0)));
        qr.slice_mut((17,17), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.01, 1.0)));

        // println!("b matrix in passive control {}", &b_flex);

        let ql = 100.0 * &qr.clone().try_inverse().unwrap();


        let a_cl = na::DMatrix::<f64>::zeros(num_states+3, num_states+3);
        let p_lyap = na::DMatrix::<f64>::zeros(num_states+3, num_states+3);
        let b_c = na::DMatrix::<f64>::zeros(num_states, 5);
        let c_c = na::DMatrix::<f64>::zeros(5, num_states+3);
        let state = na::DVector::<f64>::zeros(num_states+3);
        let u = na::DVector::<f64>::zeros(5);
        let u_rigid = na::DVector::<f64>::zeros(3);
        let enable = true;  


        let passive_control = PassiveControl {
            a_flex,
            b_flex,
            a_rigid,
            b_rigid,
            num_modes,
            num_states,
            r,
            ql,
            qr,
            a_cl,
            p_lyap,
            b_c,
            c_c,
            state,
            u,
            u_rigid,
            enable
        };
        return passive_control

    }
    
    
    pub fn care(&mut self) -> (){

        // const num_modes: usize = 2;
        // const num_states: usize = 2*num_modes;
        
        // let a_full = init_amat();
        // let b_full = init_bmat();
        

        // let mut a = DMatrix::<f64>::zeros(self.num_states, self.num_states);
        // let mut b = DMatrix::<f64>::zeros(self.num_states,5);50

        // println!("num states {}", &self.num_states);

        let a_use = self.a_flex.slice((0, 0), (self.num_states,self.num_states)).clone();
        let b_use = self.b_flex.slice((0, 0), (self.num_states,5)).clone();

        let mut a = na::DMatrix::zeros(self.num_states + 3, self.num_states + 3);
        let mut b = na::DMatrix::zeros(self.num_states + 3, 5);

        a.slice_mut((3,3), (self.num_states, self.num_states)).copy_from(&a_use);
        // a.slice_mut((0,0), (3, 3)).copy_from(&self.a_rigid);

        b.slice_mut((0,0), (3, 5)).copy_from(&self.b_rigid);
        b.slice_mut((3,0), (self.num_states, 5)).copy_from(&b_use);

        // println!("The b matrix for comparison {}", &b);
        // println!("The a matrix for comparison {}", &a);


        let num_states = &self.num_states + 3;
        // let b = self.b_full.fixed_slice::<self.num_states,5>(0, 0).clone();
    
        let z12 = -&b * self.r.try_inverse().unwrap() * &b.transpose();
        let mut z = na::DMatrix::<f64>::zeros({2*num_states}, {2*(num_states)});
       
        let qr = self.qr.clone();
        // println!("qr is {}", &qr);

        z.slice_mut((0,0), (num_states,num_states)).copy_from(&a);
        z.slice_mut((0,num_states), (num_states,num_states)).copy_from(&z12);
        z.slice_mut((num_states,0), (num_states,num_states)).copy_from(&-qr);
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
        let k = self.r.try_inverse().unwrap() * b.transpose() * p;
        
        for (loc, el) in k.iter().enumerate(){
            self.c_c[loc] = el.clone();
        }
        // println!("k is : {:.6}, c is: {:.6}",k, self.c_c);
    
        //////////////////////////////////
        //
        //
        //
        //
        //
        //
        //
        //
        //
        
    
        self.a_cl = a - (b*k);
        // println!("a_cl = {:.5}", &self.a_cl);
        
    }
    
    pub fn lyap(&mut self) -> (){

        let num_states = self.num_states + 3;
        
        let nrows = num_states * num_states;
        let nrows_sq = num_states;
    
        let qv = na::DMatrix::<f64>::from_iterator(nrows, 1, self.ql.transpose().iter().cloned());
    
        let eye = na::DMatrix::<f64>::identity(nrows_sq as usize, nrows_sq as usize);
    
        let temp1 = eye.kronecker(&self.a_cl.transpose());
        let temp2 = self.a_cl.transpose().conjugate().kronecker(&eye);
        // println!("{} {}", temp1, temp2);
    
        let temp3 = temp1 + temp2;
    
        // println!("{}", temp3);
    
        let pv = temp3.try_inverse().unwrap() * qv;
        let mut p = na::DMatrix::<f64>::identity(nrows_sq as usize, nrows_sq as usize);
    
        for (pos,el) in pv.iter().enumerate(){
            p[pos] = el.clone();
        }
        self.p_lyap = -p.clone();

        self.b_c = -p.try_inverse().unwrap() * self.c_c.clone().transpose();
        
        // println!("p lyap {:.5}", &self.p_lyap);
        // println!("bc {:.5}", &self.b_c);
    
    }

    pub fn propogate_control_state(&mut self, gyro:  &[f64], dt: f64, num_steps: u16, gyro_des: &[f64]){
        let mut gyro_v = na::DVector::from_row_slice(gyro);
        let gyro_d = na::DVector::from_row_slice(gyro_des);

        // let mut gyro_v = &gyro;

        // if &gyro_v.norm() > &0.0001 {
        //     gyro_v = gyro_v * 0.0;
        // }
        // println!("{}", self.eta);
        // let now = Instant::now();
        let num_states = self.num_states+3;
        let mut k1 = na::DVector::<f64>::zeros(num_states);
        let mut k2 = na::DVector::<f64>::zeros(num_states);
        let mut k3 = na::DVector::<f64>::zeros(num_states);
        let mut k4 = na::DVector::<f64>::zeros(num_states);

        let tau_contrib = self.b_c.clone() * (&gyro_d - &gyro_v);

        for step in 0..num_steps{
            let state0 = self.state.clone();
            let a_cl = self.a_cl.clone();
            k1 = ((&a_cl*&state0) + &tau_contrib) * dt;
            k2 = ((&a_cl*(&state0+(&k1/2.0))) + &tau_contrib) * dt;
            k3 = ((&a_cl*(&state0+(&k2/2.0))) + &tau_contrib) * dt;
            k4 = ((&a_cl*(&state0+&k3)) + &tau_contrib) * dt;
            

            let state_dot = (k1+(2.0*k2)+(2.0*k3)+k4)/6.0;
            self.state = state0 + state_dot;
        }
        let tau_des = -&self.c_c * &self.state;
        // println!("state {} tau_des {} gyros {} gyro_des {}", &self.state, &tau_des, &gyro_v, &gyro_d);
        self.u = tau_des.clone();
        // let u_r = na::DVector::from_row_slice(&[tau_des[0], ((tau_des[1] + tau_des[2])/2.0), ((tau_des[3] + tau_des[4])/2.0)]);
        // self.u_rigid = u_r;
        // println!("tau in fc: {} \n tau rigid: {}", tau_des.clone(), self.u_rigid.clone());
    }
}