// import crate for matrix math
// pub use ndarray::{arr2, Array2};
pub use nalgebra::{DMatrix, Matrix3};
// Output string to write the simulation data to csv file
#[derive(Serialize)]
#[serde(rename_all = "PascalCase")]
pub struct AutomaticaWrite {
    pub time: f64,
    pub x1: f64,
    pub x2: f64,
    pub x3: f64,
    pub x4: f64,
    pub x5: f64,
    pub x6: f64,
    pub y1: f64,
    pub u: f64,
}

// Defining a struct to hold the definitions of the system
// matrices 
pub struct SSmodel {
    pub _mass_mat:  DMatrix<f64>,
    pub _stiff_mat: DMatrix<f64>,
    pub _a_mat: DMatrix<f64>,
    pub _b_mat: DMatrix<f64>,
    pub _c_mat: DMatrix<f64>,
    pub _d_mat: DMatrix<f64>,
    
}