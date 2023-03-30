extern crate nalgebra as na;

mod amat;
mod bmat;

pub use amat::*;
pub use bmat::*;

pub struct PassiveControl {
    num_modes: usize,
    num_states: usize,
    a_full: na::SMatrix<f64, 104, 104>,
    b_full: na::SMatrix<f64, 104, 5>,
    qr: na::DMatrix<f64>,
    r: na::SMatrix<f64, 5, 5>,
    ql: na::DMatrix<f64>,
    a_cl: na::DMatrix<f64>,
    p_lyap: na::DMatrix<f64>,
    b_c: na::DMatrix<f64>,
    c_c: na::DMatrix<f64>,
    state: na::DVector<f64>,
    u: na::DVector<f64>,
}

impl PassiveControl{
    pub fn init_pc() -> PassiveControl{
        let a_full = init_amat();
        let b_full = init_bmat();
        let num_modes: usize = 2;
        let num_states: usize = 2*num_modes;
        let r = na::Matrix5::<f64>::identity();
        let ql = na::DMatrix::<f64>::identity(num_states, num_states);
        let qr = na::DMatrix::<f64>::identity(num_states, num_states);
        let a_cl = na::DMatrix::<f64>::zeros(num_states, num_states);
        let p_lyap = na::DMatrix::<f64>::zeros(num_states, num_states);
        let b_c = na::DMatrix::<f64>::zeros(num_states, 5);
        let c_c = na::DMatrix::<f64>::zeros(5, num_states);
        let state = na::DVector::<f64>::zeros(num_states);
        let u = na::DVector::<f64>::zeros(5);
        


        let passive_control = PassiveControl {
            a_full,
            b_full,
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
            u
        };
        return passive_control

    }
    
    
    pub fn care(&mut self) -> (){

        // const num_modes: usize = 2;
        // const num_states: usize = 2*num_modes;
        
        // let a_full = init_amat();
        // let b_full = init_bmat();
        

        // let mut a = DMatrix::<f64>::zeros(self.num_states, self.num_states);
        // let mut b = DMatrix::<f64>::zeros(self.num_states,5);

            
        let a = self.a_full.slice((0, 0), (self.num_states,self.num_states)).clone();
        let b = self.b_full.slice((0, 0), (self.num_states,5)).clone();

        // let b = self.b_full.fixed_slice::<self.num_states,5>(0, 0).clone();
    
        let z12 = -b * self.r.try_inverse().unwrap() * b.transpose();
        let mut z = na::DMatrix::<f64>::zeros({2*self.num_states}, {2*self.num_states});
       
        let qr = self.qr.clone();

        z.slice_mut((0,0), (self.num_states,self.num_states)).copy_from(&a);
        z.slice_mut((0,self.num_states), (self.num_states,self.num_states)).copy_from(&z12);
        z.slice_mut((self.num_states,0), (self.num_states,self.num_states)).copy_from(&-qr);
        z.slice_mut((self.num_states,self.num_states), (self.num_states,self.num_states)).copy_from(&-a.transpose());
    
        let decom = z.clone().schur();
    
        let eigen_vals = decom.complex_eigenvalues().clone();
    
        let mut z_c = na::DMatrix::<na::Complex<f64>>::zeros((2*self.num_states), (2*self.num_states));
    
        for (el_n, el) in z.clone().iter().enumerate() {
            z_c[el_n] = na::Complex::<f64>::new(el.clone(), 0.0);
        }
    
    
        // println!("{:.3} {:.3} {}", z, z_c, eigen_vals);
        let mut cnt = 0;
        let mut u11 = na::DMatrix::<na::Complex<f64>>::zeros({ self.num_states }, { self.num_states });
        let mut u12 = na::DMatrix::<na::Complex<f64>>::zeros({ self.num_states }, { self.num_states });
        
        for value in eigen_vals.iter(){
            // println!("{}", value);
            if value.re < 0.0{
    
                let mut temp = na::DMatrix::<na::Complex<f64>>::identity({ 2*self.num_states } , { 2*self.num_states });
    
    
                //Create the lamba * identity matrix
    
                for (loc, el) in temp.clone().iter().enumerate(){
                    temp[loc] = temp[loc] * value;
                }
    
                let temp2 = &z_c - temp;
    
                let svd_d = temp2.svd(true, true);
    
                let v = svd_d.v_t.unwrap().transpose();
                let mut eig_vn = na::DMatrix::<na::Complex<f64>>::zeros({ 2*self.num_states } , 1);
    
                for (col_cnt, col) in v.column_iter().enumerate(){
                    if col_cnt == ((2*self.num_states)-1){
                        eig_vn = na::DMatrix::<na::Complex<f64>>::from_column_slice({2*self.num_states}, 1, &col.as_slice());
                    }
                }
                u11.slice_mut((0,cnt), ({self.num_states},1)).copy_from(&eig_vn.slice((0, 0), (self.num_states,1)));
                u12.slice_mut((0,cnt), ({self.num_states},1)).copy_from(&eig_vn.slice((self.num_states, 0), (self.num_states,1)));
                
                cnt += 1;
    
            }
    
        }
    
        let Pc = u12 * u11.try_inverse().unwrap();
    
        let mut p = na::DMatrix::<f64>::zeros(self.num_states, self.num_states);
    
        let mut cn = 0;
        for el in Pc.iter(){
            p[cn] = el.re;
            cn +=1;
        }
    
        // println!("THIS IS P {:+.5e}", P);
        let k = self.r.try_inverse().unwrap() * b.transpose() * p;
        
        for (loc, el) in k.iter().enumerate(){
            self.c_c[loc] = -el;
        }
        // println!("k is : {}, c is: {}",k, self.c_c);
    
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
        
    
        self.a_cl = a + (b*k);
        
    }
    
    pub fn lyap(&mut self) -> (){
        
        let nrows = self.num_states * self.num_states;
        let nrows_sq = self.num_states;
    
        let qv = na::DMatrix::<f64>::from_iterator(nrows, 1, self.ql.transpose().iter().cloned());
    
        let eye = na::DMatrix::<f64>::identity(nrows_sq as usize, nrows_sq as usize);
    
        let temp1 = eye.kronecker(&self.a_cl);
        let temp2 = self.a_cl.conjugate().kronecker(&eye);
        // println!("{} {}", temp1, temp2);
    
        let temp3 = temp1 + temp2;
    
        // println!("{}", temp3);
    
        let pv = temp3.try_inverse().unwrap() * qv;
        let mut p = na::DMatrix::<f64>::identity(nrows_sq as usize, nrows_sq as usize);
    
        for (pos,el) in pv.iter().enumerate(){
            p[pos] = el.clone();
        }
        self.p_lyap = p.clone();

        self.b_c = p.try_inverse().unwrap() * self.c_c.clone().transpose();
        
        // println!("{}", p);
    
    }

    pub fn propogate_flex(&mut self, gyro:  &[f64; 5], dt: f64, num_steps: u16){
        let gyro_v = na::DVector::from_row_slice(gyro);
        // println!("{}", self.eta);
        // let now = Instant::now();
        let mut k1 = na::DVector::<f64>::zeros(self.num_states);
        let mut k2 = na::DVector::<f64>::zeros(self.num_states);
        let mut k3 = na::DVector::<f64>::zeros(self.num_states);
        let mut k4 = na::DVector::<f64>::zeros(self.num_states);

        let tau_contrib = self.b_c.clone() * gyro_v;

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
    }

}