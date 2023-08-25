extern crate nalgebra as na;

mod amat;
mod bmat;

pub type Vector8 = na::SMatrix<f64, 8, 1>;
pub type Matrix48 = na::SMatrix<f64, 4, 8>;
pub type Matrix8 = na::SMatrix<f64, 8, 8>;

#[derive(Debug)]
pub struct LQRController {
    // LQR parameters
    lqr_gain: Matrix48, // LQR gain matrix
    control_input: na::Vector4<f64>, // Control input
    reference: Vector8, // Reference state
    q: Matrix8, // State weighting matrix
    r: na::Matrix4<f64>, // Control weighting matrix
    a: amat::AMat,
    b: bmat::BMat,
}

impl LQRController {
    // Initialize the LQR controller with default values
    pub fn init() -> LQRController {
        let a = amat::init_amat();
        println!("a {}", &a);
        let b = bmat::init_bmat();
        println!("b {}", &b);
        let q = Matrix8::identity();
        let r = na::Matrix4::<f64>::identity();
        let control_input = na::Vector4::<f64>::zeros();
        let lqr_gain = Matrix48::zeros();
        let reference = Vector8::zeros();

        let mut lqr = LQRController {
            lqr_gain, // Initialized to zeros, should be computed using a Riccati solver
            control_input,
            reference,
            q,
            r,
            a,
            b,
        };
  

        lqr.care();
        println!("gain {:.4}", &lqr.lqr_gain);
        return lqr
    }

    // Set the reference state
    pub fn set_reference(&mut self, reference: &[f64]) {
        self.reference = Vector8::from_row_slice(reference);
    }

    // Calculate the control input using LQR
    pub fn calculate_control(&mut self, est_state: &[f64]) {
        // Calculate the control input using the formula: u = -K * (x - x_ref)
        let estimated_state = Vector8::from_row_slice(est_state);
        self.control_input = -self.lqr_gain * (estimated_state - &self.reference);
    }

    // Get the control input
    pub fn get_control_input(&self) -> na::Vector4<f64> {
        self.control_input.clone()
    }

    pub fn care(&mut self) {
        let a = self.a.clone();
        let b = self.b.clone();
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
            self.lqr_gain[loc] = el.clone();
        }       
    }
}
