//! The Miscellaneous submodule details functions used throught the control and estimation modeles/algorithms.
//!
//! This crate implements a a series of function that are used through much of the control and estimation
//! process. THese functions include things such as `unxmat` and `grab_vec3`
extern crate nalgebra as na;

/// Function to grab and copy the contents of an `na::Vector3<f64>` type into another location.
///
/// # Arguments
///
/// `v2: &na::Vector3<f64>` - the vector slice to be copied
/// `v1: &mut na::Vector3<f64>` - and mutable vector slice to be copied into
#[allow(dead_code)]
pub fn grab_vec3(v1: &mut na::Vector3<f64>, v2: &na::Vector3<f64>) -> () {
    for _k in 0..3 as usize {
        v1[_k] = v2[_k];
    }
}
/// Function to deconstruct a cross matrix.
///
/// # Details
/// This function takes in a crossed matrix and reconstructs the original vector being crossed
///
/// # Arguments
///
/// `rot: &na::Matrix3<f64>` - the matrix to be deconstructed
/// `phi: &mut na::Vector3<f64>` - the vector to be reconstructed
pub fn unxmat(rot: &na::Matrix3<f64>, phi: &mut na::Vector3<f64>) -> () {
    phi[0] = (rot[(2, 1)] - rot[(1, 2)]) / 2.0;
    phi[1] = (rot[(0, 2)] - rot[(2, 0)]) / 2.0;
    phi[2] = (rot[(1, 0)] - rot[(0, 1)]) / 2.0;
}
/// Function to grenerate the cross matric from a vector
///
/// # Details
///
/// This function takes in a vector and generates the cross matrix from it
///
/// # Arguments
///
/// `phi: &na::Vector3<f64>` - the vector to be crossed
///
/// # Results
///
/// `rot: na::Matrix3<f64>` - the cross matrix
pub fn xmat(phi: &na::Vector3<f64>, rot: &mut na::Matrix3<f64>) {
    rot[(0, 1)] = -phi[2];
    rot[(1, 0)] = phi[2];
    rot[(0, 2)] = phi[1];
    rot[(2, 0)] = -phi[1];
    rot[(1, 2)] = -phi[0];
    rot[(2, 1)] = phi[0];
    rot[(0, 0)] = 0.0;
    rot[(1, 1)] = 0.0;  
    rot[(2, 2)] = 0.0;
}
