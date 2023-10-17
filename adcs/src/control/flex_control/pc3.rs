extern crate nalgebra as na;
use crate::control::gains::{Gains};

use crate::control::flex_control as fc;
use crate::control::lqr;

type BRigid = na::SMatrix<f64, 9, 4>;
type BRigid2 = na::SMatrix<f64, 3, 4>;
type Matrix1 = na::SMatrix<f64,1,1>;
type Vector9 = na::SMatrix<f64, 9, 1>;
type Matrix49 = na::SMatrix<f64, 4, 9>;
type Matrix9 = na::SMatrix<f64, 9, 9>;


#[derive(Clone)]
pub struct PassiveControl {
    num_modes: usize,
    num_states: usize,
    a_flex: na::SMatrix<f64, 104, 104>,
    b_flex: na::SMatrix<f64, 104, 5>,
    a: na::SMatrix<f64,9,9>,
    b: na::SMatrix<f64, 9,4>,
    qr: na::DMatrix<f64>,
    r: na::SMatrix<f64, 4, 4>,
    ql: na::DMatrix<f64>,
    a_cl: na::DMatrix<f64>,
    p_lyap: na::DMatrix<f64>,
    b_c: na::DMatrix<f64>,
    c_c: na::DMatrix<f64>,
    state: na::DVector<f64>,
    tau_max: f64,
    pub u: na::DVector<f64>, //the requested torque
    pub u_rigid: na::DVector<f64>,
    pub enable: bool,
    pub update: bool,
}

impl PassiveControl{
    pub fn init_pc_p() -> PassiveControl{
        let a_flex = fc::amat::init_amat();
        let b_flex = 1000.0 *fc::bmat::init_bmat();
        let a = lqr::amat::init_amat_3();
        // println!("a {}", &a);
        let b = lqr::bmat::init_bmat_3();
        // println!("b {}", &b);

        // let qr_rigid = Matrix9::from_diagonal(&Vector9::from_row_slice(
        //     &[
        //         99993.6809181096,	99999.9638945158,	0.00100000000000000,
        //         96371.8875195649,	22735.6443000617,	12605.1625984226,
        //         0.00100000000000000,	19121.1224821913,	0.00100000000000000
        //     ]));

        let qr_rigid = Matrix9::from_diagonal(&Vector9::from_row_slice(
            &[
                0.00001,	
                454.0,
                454.0,
                100.0,	
                0.00001,
                1000000.0,
                1000000.0,
                10000.0,
                0.0000010
            ]));


        let r =  na::Matrix4::from_diagonal(&na::Vector4::from_row_slice(
            &[
            0.01, //RW
            0.05, // Roll
            0.1, // Pitch
            100000.0 // piv
        ]));
   
        // println!("b_rigid: {}", &b_rigid);
        let num_modes: usize = 0;
        let num_rigid: usize = 9;
        let num_states: usize = 2*num_modes+num_rigid;

        let mut ql = 1.0*na::DMatrix::<f64>::identity(num_states, num_states);
        let mut qr = 1.0 * na::DMatrix::<f64>::identity(num_states, num_states);

        qr.slice_mut((0,0), (9, 9)).copy_from(&(qr_rigid));
        // println!("check 1 ");
        // qr.slice_mut((3,3), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 100.0)));
        // qr.slice_mut((5,5), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1000.0)));
        // qr.slice_mut((6,6), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((8,8), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((10,10), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((12,12), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((14,14), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((16,16), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((18,18), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((20,20), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((22,22), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((24,24), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));

        // println!("b matrix in passive control {}", &b_flex);
        // let ql_rigid = 0.1*Matrix9::from_diagonal(&Vector9::from_row_slice(
            //     &[5384.60616938582,	35080.9199512616,	3644.14765342167,
            //     100000.0,	7326.42098644794,	46169.9196801448,
            //     89699.3233076463,	68705.1995568645,	99847.0993272210
            //     ]));
    
    
            let ql_rigid = 1.0*Matrix9::from_diagonal(&Vector9::from_row_slice(
                    &[0.00100000000000000,
                    707.372313457858,
                    0.00100000000000000,
                    0.00100000000000000,
                    0.00100000000000000,
                    0.00100000000000000,
                    6366.01251438075,
                    0.00100000000000000,
                    1721.27686908652
                    ]));

            ql.slice_mut((0,0), (9, 9)).copy_from(&(qr_rigid));
        // println!("{}", &ql);


        let a_cl = na::DMatrix::<f64>::zeros(num_states, num_states);
        let p_lyap = na::DMatrix::<f64>::zeros(num_states, num_states);
        let b_c = na::DMatrix::<f64>::zeros(num_states, 4);
        let c_c = na::DMatrix::<f64>::zeros(4, num_states);
        let state = na::DVector::<f64>::zeros(num_states);
        let u = na::DVector::<f64>::zeros(4);
        let u_rigid = na::DVector::<f64>::zeros(3);
        let tau_max = 40.0;
        let enable = false;
        let update = false;  
        let ql2 = ql.clone();
        let qr2 = qr.clone();

        let r_sc_slew = 1.0;
        let qr_sc_slew = 0.01;
        let ql_sc_slew = 100.0;

        let mut passive_control = PassiveControl {
            a_flex,
            b_flex,
            a,
            b,
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
            tau_max,
            u,
            u_rigid,
            enable,
            update
        };
        passive_control.recalc(&r, &ql2, &qr2);

        return passive_control

    }

    pub fn init_pc_n() -> PassiveControl{
        let a_flex = fc::amat::init_amat();
        let b_flex = 1000.0 *fc::bmat::init_bmat();
        let a = lqr::amat::init_amat_3n();
        // println!("a {}", &a);
        let b = lqr::bmat::init_bmat_3n();
        // println!("b {}", &b);


        // let qr_rigid = Matrix9::from_diagonal(&Vector9::from_row_slice(
        //     &[
        //         99993.6809181096,	99999.9638945158,	0.00100000000000000,
        //         96371.8875195649,	22735.6443000617,	12605.1625984226,
        //         0.00100000000000000,	19121.1224821913,	0.00100000000000000
        //     ]));

        let qr_rigid = Matrix9::from_diagonal(&Vector9::from_row_slice(
            &[
                0.00001,	
                454.0,
                454.0,
                100.0,	
                0.00001,
                1000000.0,
                1000000.0,
                10000.0,
                0.0000010
            ]));

        let r =  na::Matrix4::from_diagonal(&na::Vector4::from_row_slice(
            &[
            0.01, //RW
            0.05, // Roll
            0.1, // Pitch
            100000.0 // piv
        ]));
   
        // println!("b_rigid: {}", &b_rigid);
        let num_modes: usize = 0;
        let num_rigid: usize = 9;
        let num_states: usize = 2*num_modes+num_rigid;

        let mut ql = 1.0*na::DMatrix::<f64>::identity(num_states, num_states);
        let mut qr = 1.0 * na::DMatrix::<f64>::identity(num_states, num_states);

        qr.slice_mut((0,0), (9, 9)).copy_from(&(qr_rigid));
        // println!("check 1 ");
        // qr.slice_mut((3,3), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 100.0)));
        // qr.slice_mut((5,5), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1000.0)));
        // qr.slice_mut((6,6), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((8,8), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((10,10), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((12,12), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((14,14), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((16,16), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((18,18), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((20,20), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((22,22), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));
        // qr.slice_mut((24,24), (2, 2)).copy_from(&na::Matrix2::<f64>::from_diagonal(&na::Vector2::new(0.000001, 1.0)));

        // println!("b matrix in passive control {}", &b_flex);
        // let ql_rigid = 0.1*Matrix9::from_diagonal(&Vector9::from_row_slice(
        //     &[5384.60616938582,	35080.9199512616,	3644.14765342167,
        //     100000.0,	7326.42098644794,	46169.9196801448,
        //     89699.3233076463,	68705.1995568645,	99847.0993272210
        //     ]));


        let ql_rigid = Matrix9::from_diagonal(&Vector9::from_row_slice(
                &[0.00100000000000000,
                707.372313457858,
                0.00100000000000000,
                0.00100000000000000,
                0.00100000000000000,
                0.00100000000000000,
                6366.01251438075,
                0.00100000000000000,
                1721.27686908652
                ]));

            ql.slice_mut((0,0), (9, 9)).copy_from(&(qr_rigid));
        // println!("{}", &ql);
        

        let a_cl = na::DMatrix::<f64>::zeros(num_states, num_states);
        let p_lyap = na::DMatrix::<f64>::zeros(num_states, num_states);
        let b_c = na::DMatrix::<f64>::zeros(num_states, 4);
        let c_c = na::DMatrix::<f64>::zeros(4, num_states);
        let state = na::DVector::<f64>::zeros(num_states);
        let u = na::DVector::<f64>::zeros(4);
        let u_rigid = na::DVector::<f64>::zeros(3);
        let tau_max = 10.0;
        let enable = false;
        let update = false;  
        let ql2 = ql.clone();
        let qr2 = qr.clone();

        let r_sc_slew = 1.0;
        let qr_sc_slew = 0.01;
        let ql_sc_slew = 100.0;

        let mut passive_control = PassiveControl {
            a_flex,
            b_flex,
            a,
            b,
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
            tau_max,
            u,
            u_rigid,
            enable,
            update
        };
        passive_control.recalc(&r, &ql2, &qr2);

        return passive_control

    }


    pub fn reset(&mut self){
        self.state = na::DVector::<f64>::zeros(self.num_states);
    }
    
    pub fn recalc(&mut self, r: &na::SMatrix<f64, 4, 4>, ql: &na::DMatrix<f64>, qr: &na::DMatrix<f64>) {
        self.care(r, qr);
        self.lyap(ql);
    }
    
    pub fn care(&mut self, r: &na::SMatrix<f64, 4, 4>, qr: &na::DMatrix<f64>) {

        // const num_modes: usize = 2;
        // const num_states: usize = 2*num_modes;
        
        // let a_full = init_amat();
        // let b_full = init_bmat();
        

        // let mut a = DMatrix::<f64>::zeros(self.num_states, self.num_states);
        // let mut b = DMatrix::<f64>::zeros(self.num_states,5);50

        // println!("num states {}", &self.num_states);

        // let a_use = self.a_flex.slice((0, 0), (self.num_modes*2,self.num_modes * 2)).clone();
        // let b_use = self.b_flex.slice((0, 0), (self.num_states,5)).clone();

        let mut a = na::DMatrix::zeros(self.num_states, self.num_states);
        let mut b = na::DMatrix::zeros(self.num_states, 4);

        // a.slice_mut((8,8), (self.num_states, self.num_states)).copy_from(&a_use);
        a.slice_mut((0,0), (9, 9)).copy_from(&self.a);

        for loc in 0..self.num_modes{
            let loc2 = (loc*2)+1;
            let val = a[(loc2,loc2)].clone()/10.0;
            a.slice_mut((loc2,loc2), (1,1)).copy_from(&(Matrix1::from_row_slice(&[val])));
        }
        // println!("{}", &a);

        b.slice_mut((0,0), (9, 4)).copy_from(&self.b);
        // b.slice_mut((6,0), (self.num_states, 5)).copy_from(&b_use);

        // println!("The b matrix for comparison {}", &b);
        // println!("The a matrix for comparison {}", &a);


        let num_states = self.num_states.clone();
        // let b = self.b_full.fixed_slice::<self.num_states,5>(0, 0).clone();
    
        let z12 = -&b * r.try_inverse().unwrap() * &b.transpose();
        let mut z = na::DMatrix::<f64>::zeros({2*num_states}, {2*(num_states)});
       
        let qr = qr.clone();
        println!("qr is {}", &qr);

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
        let k = r.try_inverse().unwrap() * b.transpose() * p;
        
        for (loc, el) in k.iter().enumerate(){
            self.c_c[loc] = el.clone();
        }
        println!("k is : {:.3}, c is: {:.3}",k, self.c_c);
    
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
        println!("a_cl = {:.5}", &self.a_cl);
        
    }
    
    pub fn lyap(&mut self, ql: &na::DMatrix<f64>) {

        let num_states = self.num_states.clone();
        
        let nrows = num_states * num_states;
        let nrows_sq = num_states;
    
        let qv = na::DMatrix::<f64>::from_iterator(nrows, 1, ql.transpose().iter().cloned());
    
        let eye = na::DMatrix::<f64>::identity(nrows_sq as usize, nrows_sq as usize);
    
        let temp1 = eye.kronecker(&self.a_cl.transpose());
        let temp2 = self.a_cl.transpose().conjugate().kronecker(&eye);
        // println!("{} {}", temp1, temp2);
        // println!("Testing\n");
        let temp3 = temp1 + temp2;
    
        // println!("{}", temp3);
    
        let pv = temp3.try_inverse().unwrap() * qv;
        let mut p = na::DMatrix::<f64>::identity(nrows_sq as usize, nrows_sq as usize);
        println!("Testing\n");
        for (pos,el) in pv.iter().enumerate(){
            p[pos] = el.clone();
        }
        self.p_lyap = -p.clone();
        // println!("Testing\n");
        self.b_c = -p.try_inverse().unwrap() * self.c_c.clone().transpose();
        
        println!("p lyap {:.3}", &self.p_lyap);
        println!("bc {:.3}", &self.b_c);
    
    }

    pub fn propogate_control_state(&mut self, gyro:  &[f64], dt: f64, num_steps: u16, gyro_des: &[f64]){
        let mut gyro_v = na::DVector::from_row_slice(gyro);
        let gyro_d = na::DVector::from_row_slice(gyro_des);
        // println!("Testing\n");
        // let mut gyro_v = &gyro;

        // if &gyro_v.norm() > &0.0001 {
        //     gyro_v = gyro_v * 0.0;
        // }
        // println!("{}", self.eta);
        // let now = Instant::now();
        let num_states = self.num_states;
        let mut k1 = na::DVector::<f64>::zeros(num_states);
        let mut k2 = na::DVector::<f64>::zeros(num_states);
        let mut k3 = na::DVector::<f64>::zeros(num_states);
        let mut k4 = na::DVector::<f64>::zeros(num_states);
        // println!("Testing\n");
        let tau_contrib = -self.b_c.clone() * (&gyro_v - &gyro_d);
        // let tau_contrib = self.b_c.clone() * (&gyro_d - (&gyro_v));

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
        // let tau_des = -&self.c_c * &self.state;
        let tau_des = &self.c_c * &self.state;
        // println!("Testing\n");
        // println!("state {} tau_des {} gyros {} gyro_des {}", &self.state, &tau_des, &gyro_v, &gyro_d);
        ////// BOUND TORQUE
        let mut tau_bound = tau_des.clone();
        // println!("trque before: {}", &tau_des);
        for (loc) in (0..4){
            let torque = tau_des[loc].clone();
            let mag = torque.abs();
            if mag > self.tau_max {
                tau_bound[loc] = torque.signum() * self.tau_max;
            }
            if loc == 0 {
                tau_bound[loc] = -tau_bound[loc];
            }

        }
        // println!("Testing\n");
        // println!("trque after: {}", &tau_bound);
        ///////
        self.u = tau_bound.clone();
        // let u_r = na::DVector::from_row_slice(&[tau_des[0], ((tau_des[1] + tau_des[2])/2.0), ((tau_des[3] + tau_des[4])/2.0)]);
        // self.u_rigid = u_r;
        // println!("tau in fc: {} \n tau rigid: {}", tau_des.clone(), self.u_rigid.clone());
    }
}