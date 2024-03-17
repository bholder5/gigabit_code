extern crate nalgebra as na;

pub type AMat1 = na::SMatrix<f64, 2, 2>;

pub fn init_AMat1_pos() -> AMat1 {
    let a = AMat1::from_row_slice(&[
0.0000000000000000,1.0000000000000000,
-0.0429919967711452,-0.0000000000000000,
    ]);
    return a;
}

pub type AMat2 = na::SMatrix<f64, 2, 2>;

pub fn init_AMat2_pos() -> AMat2 {
    let a = AMat2::from_row_slice(&[
0.0000000000000000,1.0000000000000000,
1.0000000000000000,0.0000000000000000,
    ]);
    return a;
}

pub type AMat3 = na::SMatrix<f64, 1, 4>;

pub fn init_AMat3_pos() -> AMat3 {
    let a = AMat3::from_row_slice(&[
0.0010283592132166,1.4165598180369721,-0.0000000715310370,-0.0000492536153273,
    ]);
    return a;
}

pub type AMat4 = na::SMatrix<f64, 1, 4>;

pub fn init_AMat4_pos() -> AMat4 {
    let a = AMat4::from_row_slice(&[
0.0000000197911747,0.0031560216655633,0.0000000010325766,0.0004925458262128,
    ]);
    return a;
}

pub type AMat5 = na::SMatrix<f64, 1, 4>;

pub fn init_AMat5_pos() -> AMat5 {
    let a = AMat5::from_row_slice(&[
0.0000000084817130,0.0015880896574600,0.0000000015247988,0.0004925460100207,
    ]);
    return a;
}

