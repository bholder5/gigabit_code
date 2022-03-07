
// import crate for matrix math
// pub use ndarray::{arr2, Array2};

pub use std::ffi::OsString;
pub use std::fs::OpenOptions;
pub use std::path::Path;
pub use std::env;
use std::error::Error;
use csv::WriterBuilder;
use crate::types::AutomaticaWrite;
pub use nalgebra::{DMatrix, Matrix3, Matrix};

//Assemble and return the mass matrix based on the rigid body and 
// flexible moi constants
pub fn build_mass(moi_ir: f64, moi_i1: f64, moi_i2: f64)-> DMatrix<f64>{
    let _mass_mat: DMatrix<f64>= DMatrix::<f64>::from_vec(3,3, vec![moi_ir+moi_i1+moi_i2, moi_i1, moi_i2,
                                                                          moi_i1,               moi_i1,    0.0,
                                                                          moi_i2,                  0.0, moi_i2]);
        return _mass_mat;
}

//Assemble and return the stiffness matrix based on the  
// flexible mode k values.
pub fn build_stiffness(k1: f64, k2: f64)-> DMatrix<f64>{
    let _stiff_mat: DMatrix<f64>= DMatrix::<f64>::from_vec(3,3,vec![0.0, 0.0, 0.0, 
                                                                          0.0, k1,  0.0,
                                                                          0.0, 0.0, k2]);
        return _stiff_mat;
}

// a function to insert a small matrix (sub matrix) in block format
// into a larger matrix
pub fn mat_block_change(mut _mat_s: DMatrix<f64>, mut _mat_l: DMatrix<f64>, _blk: usize)-> DMatrix<f64>{
    // find idices
    let _dim_mat = (_mat_l.nrows(), _mat_l.ncols());
    let _dim_sub = (_mat_s.nrows(), _mat_s.ncols());
    let mut rw: Vec<usize> = (0.._dim_sub.0).collect();
    let mut cl: Vec<usize> = (0.._dim_sub.1).collect();

    //debug purposes, should make this safer with proper error message when dims arent compatible
    // println!("{:?}, {:?}, {:?}, {:?}", _dim_mat, _dim_sub, rw, cl);
    
    if _blk == 2{
        cl = (_dim_mat.1 - _dim_sub.1.._dim_mat.1).collect();
    } else if _blk == 3{
        rw = (_dim_mat.0 - _dim_sub.0.._dim_mat.0).collect();
    } else if _blk == 4{
        cl = (_dim_mat.1 - _dim_sub.1.._dim_mat.1).collect();
        rw = (_dim_mat.0 - _dim_sub.0.._dim_mat.0).collect();    
    }
    
    for x in rw.clone(){
        for y in cl.clone(){
            _mat_l.row_mut(x)[y] = _mat_s.row_mut(x-rw[0])[y - cl[0]].clone();
        }
    }
    // println!("{}\n {}", _mat_l, _mat_s);
    return _mat_l;
}

// For writing csv
//
//
// a function to grab my csv file pointed to as an input on the command line
pub fn get_first_arg() -> Result<OsString, Box<dyn Error>> {
    // the first argument here is a path to the target hence actually taking the second argument
    match env::args_os().nth(1) {
        None => Err(From::from("Expected 1 argument, but got none")),
        Some(file_path) => Ok(file_path),
    }
}

// specifically writes csv for this simulation
pub fn write_func_automatica(s: &automatica::linear_system::solver::Step<f64>, u: f64, t0: f64) -> Result<(), Box<dyn Error>> {
    let file_path = get_first_arg()?;
    let mut t= s.time().to_string().parse().unwrap();
    t = t + t0;

    if Path::new(&file_path.clone()).exists() {
        let file = OpenOptions::new().append(true).open(&file_path.clone()).unwrap();
        let mut wtr = WriterBuilder::new().has_headers(false).from_writer(file);
        wtr.serialize(AutomaticaWrite {
            time: t,
            x1: s.state()[0],
            x2: s.state()[1],
            x3: s.state()[2],
            x4: s.state()[3],
            x5: s.state()[4],
            x6: s.state()[5],
            y1: s.output()[0],
            u: u,
        })?;
        wtr.flush()?;
    }
    else {
        let file = OpenOptions::new().append(true).create(true).open(&file_path.clone()).unwrap();
        let mut wtr = WriterBuilder::new().has_headers(true).from_writer(file);
        wtr.serialize(AutomaticaWrite {
            time: t,
            x1: s.state()[0],
            x2: s.state()[1],
            x3: s.state()[2],
            x4: s.state()[3],
            x5: s.state()[4],
            x6: s.state()[5],
            y1: s.output()[0],
            u: u,
        })?;
        wtr.flush()?;
    }    
    Ok(())
}