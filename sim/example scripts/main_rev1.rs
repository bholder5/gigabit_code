#![allow(unused_imports)]
// import crate for matrix math
use ndarray::{arr2, Array2};
// import for writing csv
extern crate csv;
extern crate serde;
#[macro_use]
extern crate serde_derive;
use std::error::Error;
use csv::WriterBuilder;
// use std::fs::File;
use std::process;
// use std::io;
use std::env;
use std::ffi::OsString;
use std::fs::OpenOptions;
use std::path::Path;

/// linear systems crates etc for automatica
extern crate automatica;

use num_complex::Complex;

use automatica::{signals::continuous, Seconds, Ss, TfMatrix};


struct SSmodel {
    _mass_mat:  Array2<f32>,
    _stiff_mat: Array2<f32>,
    // _A_mat: Array2<f32>,
    // _B_mat: Array2<f32>,
    // _C_mat: Array2<f32>,
}

#[derive(Serialize)]
#[serde(rename_all = "PascalCase")]
struct AutomaticaWrite {
    time: String,
    x1: f64,
    x2: f64,
    y1: f64,
    y2: f64,
}
//Assemble and return the mass matrix based on the rigid body and 
// flexible moi constants
fn build_mass(moi_ir: f32, moi_i1: f32, moi_i2: f32)-> Array2<f32>{
    let _mass_mat: Array2<f32>= arr2(&[[(moi_ir+moi_i1+moi_i2), moi_i1, moi_i2], 
                                       [                moi_i1, moi_i1,    0.0],
                                       [                moi_i2,    0.0, moi_i2]]);
        return _mass_mat;
}
//Assemble and return the stiffness matrix based on the  
// flexible mode k values.
fn build_stiffness(k1: f32, k2: f32)-> Array2<f32>{
    let _stiff_mat: Array2<f32>= arr2(&[[0.0, 0.0, 0.0], 
                                       [0.0,  k1, 0.0],
                                       [0.0, 0.0,  k2]]);
        return _stiff_mat;
}

#[allow(unused_assignments)]
#[allow(unused_variables)]
#[allow(dead_code)]
// for linear systems automatica toolbox examples
#[allow(clippy::many_single_char_names)]
#[allow(clippy::non_ascii_literal)]
fn main() {
    // Initializing the model constants for simulating a rigid body mode and 2 elastic mode yaw system
    const MOI_IR:f32 = 1.0;
    const MOI_I1:f32 = 2.0;
    const MOI_I2:f32 = 3.0;

    const K1:f32 = 4.0;
    const K2:f32 = 5.0;

    // initializing the model of the system
    let _model = SSmodel {
        _mass_mat: build_mass(MOI_IR, MOI_I1, MOI_I2),
        _stiff_mat: build_stiffness(K1, K2)};

    println!("{},\n {}", _model._mass_mat, _model._stiff_mat);

    // Check matrix multiplication results using identity matrix
    let _identity: Array2<f32>= arr2(&[[ 1.0, 0.0, 0.0], 
                                           [ 0.0, 1.0, 0.0],
                                           [ 0.0, 0.0,1.0]]);

    let _mat_mult_result: Array2<f32>= _model._mass_mat.dot(&_identity);

    // println!("{} x \n {} = \n {}", _model._mass_mat, _identity, _mat_mult_result);

    // This is my test print function, commented out and printing results of rk sim from automatica
    // if env::args_os().len() > 1{
    //     let _test_vec_x: Vec<i32> = vec![1,2,3,4,5,6,7,8,9,10];
    //     let _test_vec_y: Vec<i32> = vec![2,4,6,8,10,12,14,16,18,20];
    //     for loc in 0.._test_vec_x.len(){
    //         // serves to debug the csv write
    //         // println!("{} line", loc);
    //         if let Err(err) = write_function(_test_vec_x[loc], _test_vec_y[loc]) {
    //             println!("{}", err);
    //             process::exit(1);
    //         }
    //     }
    // } else {
    //     println!("No arugment provided, csv not written");
    // }

// a linear systems tool box test/ example code
// 
// 
// 
// 
//

    let a = [-1., 1., -1., 0.25];
    let b = [1., 0.25];
    let c = [0., 1., -1., 1.];
    let d = [0., 1.];

    let sys = Ss::new_from_slice(2, 1, 2, &a, &b, &c, &d);
    let poles = sys.poles();

    println!("\nUnitary step response:");
    let step = continuous::step(1., 1);

    let rk2 = sys.rk2(&step, &[0., 0.], Seconds(0.1), 150);

    
    // Change to 'true' to print the result
    if env::args_os().len() > 1{ 
        for i in sys.rk2(&step, &[0., 0.], Seconds(0.1), 150) {
            if let Err(err) = write_func_automatica(&i) {
                println!("{}", err);
                process::exit(1);
            }
        }
    } else {
        println!("No arugment provided, csv not written");
    }

}

// a function to grab my csv file pointed to as an input on the command line
fn get_first_arg() -> Result<OsString, Box<dyn Error>> {
    // the first argument here is a path to the target hence actually taking the second argument
    match env::args_os().nth(1) {
        None => Err(From::from("Expected 1 argument, but got none")),
        Some(file_path) => Ok(file_path),
    }
}

fn write_func_automatica(s: &automatica::linear_system::solver::Step<f64>) -> Result<(), Box<dyn Error>> {
    let file_path = get_first_arg()?;

    if Path::new(&file_path.clone()).exists() {
        let file = OpenOptions::new().append(true).open(&file_path.clone()).unwrap();
        let mut wtr = WriterBuilder::new().has_headers(false).from_writer(file);
        wtr.serialize(AutomaticaWrite {
            time: s.time().to_string(),
            x1: s.state()[0],
            x2: s.state()[1],
            y1: s.output()[0],
            y2: s.output()[1],
        })?;
        wtr.flush()?;
    }
    else {
        let file = OpenOptions::new().append(true).create(true).open(&file_path.clone()).unwrap();
        let mut wtr = WriterBuilder::new().has_headers(true).from_writer(file);
        wtr.serialize(AutomaticaWrite {
            time: s.time().to_string(),
            x1: s.state()[0],
            x2: s.state()[1],
            y1: s.output()[0],
            y2: s.output()[1],
        })?;
        wtr.flush()?;
    }    
    Ok(())
}