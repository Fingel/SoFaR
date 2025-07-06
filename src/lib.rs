pub mod constants;
pub mod vml;

// Type aliases to reduce finger strain
///"Position" or "p" vectors (or just plain 3-vectors) have dimension
///(3) in Fortran and [3] in C.
type Pvector = [f64; 3];

///"Rotation" or "r" matrices have dimensions (3,3) in Fortran and [3][3]
///in C.  When used for rotation, they are "orthogonal";  the inverse of
///such a matrix is equal to the transpose.  Most of the routines in
///this library do not assume that r-matrices are necessarily orthogonal
///and in fact work on any 3x3 matrix.
type Rmatrix = [[f64; 3]; 3];

///"Position/velocity" or "pv" vectors have dimensions (3,2) in Fortran
///and [2][3] in C.
type PVvector = [[f64; 3]; 2];

//TODO: Remove this in favor of round() as it compiles down to
// a single instruction on x86_64: ROUNDSD
/// Round to nearest whole number (double)
/// This is equivalent to Rust `x.round()`
fn dnint(x: f64) -> f64 {
    if x.abs() < 0.5 {
        0.0
    } else if x < 0.0 {
        (x - 0.5).ceil()
    } else {
        (x + 0.5).floor()
    }
}

//TODO: Remove this in favor of trunc() as it compiles down to
// a single instruction on x86_64: FISTTP
/// Truncate to the nearest whole number towards zero (double)
/// This is equivalent to Rust `x.trunc()`
fn dint(x: f64) -> f64 {
    if x < 0.0 { x.ceil() } else { x.floor() }
}

// Angles
// Angle conversions
pub use vml::angles::to_sexagesimal::days_to_hms as d2tf;
pub use vml::angles::to_sexagesimal::radians_to_hms as a2tf;

// Pvectors and Rmatrices
// Initialize
pub use vml::pvrm::initialize::identity_rmatrix as ir;
pub use vml::pvrm::initialize::zero_pvector as zp;
pub use vml::pvrm::initialize::zero_rmatrix as zr;

// Rotations
pub use vml::pvrm::rotations::rotate_rmatrix_about_x as rx;
pub use vml::pvrm::rotations::rotate_rmatrix_about_y as ry;
pub use vml::pvrm::rotations::rotate_rmatrix_about_z as rz;

// Spherical/Cartesian conversions
pub use vml::pvrm::sphere_cart_conv::p_vector_to_spherical as p2s;
pub use vml::pvrm::sphere_cart_conv::spherical_to_pvector as s2p;
pub use vml::pvrm::sphere_cart_conv::spherical_to_unit_vector as s2c;
pub use vml::pvrm::sphere_cart_conv::unit_vector_to_spherical as c2s;

// Operations on vectors
pub use vml::pvrm::vec_ops::pvector_cross_product as pxp;
pub use vml::pvrm::vec_ops::pvector_dot_product as pdp;
pub use vml::pvrm::vec_ops::pvector_minus_pvector as pmp;
pub use vml::pvrm::vec_ops::pvector_modulus as pm;
pub use vml::pvrm::vec_ops::pvector_multiply_scalar as sxp;
pub use vml::pvrm::vec_ops::pvector_normalize as pn;
pub use vml::pvrm::vec_ops::pvector_plus_pvector as ppp;
pub use vml::pvrm::vec_ops::pvector_plus_scaled_pvector as ppsp;

// Operations on matrices
pub use vml::pvrm::matrix_ops::rmatrix_multiply as rxr;
pub use vml::pvrm::matrix_ops::transpose_rmatrix as tr;

// Matrix-vector products
pub use vml::pvrm::matrix_vec_products::rmatrix_pvector_product as rxp;
pub use vml::pvrm::matrix_vec_products::transpose_rmatrix_pvector_product as trxp;

// Separation and position-angle
pub use vml::pvrm::sep_position_angle::angular_separation_pvector as sepp;
pub use vml::pvrm::sep_position_angle::angular_separation_spherical as seps;
pub use vml::pvrm::sep_position_angle::position_angle_from_pvector as pap;
pub use vml::pvrm::sep_position_angle::position_angle_from_spherical as pas;

// Rotation vectors
pub use vml::pvrm::rotation_vectors::rmatrix_to_rvector as rm2v;
pub use vml::pvrm::rotation_vectors::rvector_to_rmatrix as rv2m;

// PV Vectors
// Initialize
pub use vml::pvvec::initialize::zero_pvvector as zpv;

// Copy, extend, extract
pub use vml::pvvec::copy_extend_extract::append_zvelocity_pvvector as p2pv;
pub use vml::pvvec::copy_extend_extract::copy_pvvector as cpv;
pub use vml::pvvec::copy_extend_extract::discard_velocity_pvvector as pv2p;
