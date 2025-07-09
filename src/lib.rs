use std::fmt;

pub mod calendars;
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

#[derive(Debug)]
pub struct Warning {
    pub code: i32,
    pub message: String,
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} = {}", self.code, self.message)
    }
}

// Macro to quickly create warning
#[macro_export]
macro_rules! warning {
    ($code:expr, $($arg:tt)*) => {
        $crate::Warning{code: $code, message: format!($($arg)*)}
    };
}

/// Most SOFA functions accept return arguments in the form of mutable pointers.
/// The actual return value of these functions is either absent, or a status code.
/// 0 means success, and any other value indicates an error or warning. Many functions
/// return warning codes, but still complete their intended operations.
///
/// Because we aren't using return arguments, we need a way to indicate these warning
/// states while still returning a value. Instead of simply returning a tuple, we define
/// this warning type so that we can also include context with the error code. Most
/// importantly, it removes the ambiguity from the return type of tuple (value, status)
/// vs functions that return (value, value) and have no warnings.
#[derive(Debug)]
pub struct Warned<T> {
    pub value: T,
    pub warning: Option<Warning>,
}

impl<T> Warned<T> {
    /// Value with no warnings associated
    pub fn new(value: T) -> Self {
        Self {
            value,
            warning: None,
        }
    }

    /// Contains a warning, but also a computed value
    pub fn with_warning(value: T, warning: Warning) -> Self {
        Self {
            value,
            warning: Some(warning),
        }
    }

    /// Return 0 if no warning, otherwise return the warning code
    pub fn ok_code(&self) -> i32 {
        self.warning.as_ref().map_or(0, |w| w.code)
    }
    /// Returns an error if there is a warning
    pub fn value_safe(&self) -> Result<&T, &Warning> {
        self.warning.as_ref().map_or(Ok(&self.value), Err)
    }
}

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

// Calendar
pub use crate::calendars::{
    besselian_e_to_jd as epb2jd, cal_to_jd as cal2jd, jd_to_besselian_e as epb,
    jd_to_cal as jd2cal, jd_to_cal_fmt as jdcalf, jd_to_je as epj, je_to_jd as epj2jd,
};

// Vector Math

// Angles
// Wrap
pub use vml::angles::wrap::angle_normalize_pm as anpm;
pub use vml::angles::wrap::angle_normalize_positive as anp;

// Angle conversions
pub use vml::angles::to_sexagesimal::days_to_hms as d2tf;
pub use vml::angles::to_sexagesimal::radians_to_dms as a2af;
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

// Sphereical/Cartesian conversions
pub use vml::pvvec::sphere_cart_conv::pvvector_to_spherical as pv2s;
pub use vml::pvvec::sphere_cart_conv::spherical_to_pvvector as s2pv;

// Operations on pv-vectors
pub use vml::pvvec::pvvector_ops::pvvector_cross_pvvector as pvxpv;
pub use vml::pvvec::pvvector_ops::pvvector_dot_pvvector as pvdpv;
pub use vml::pvvec::pvvector_ops::pvvector_minus_pvvector as pvmpv;
pub use vml::pvvec::pvvector_ops::pvvector_modulus as pvm;
pub use vml::pvvec::pvvector_ops::pvvector_multiply_scalar as sxpv;
pub use vml::pvvec::pvvector_ops::pvvector_multiply_two_scalar as s2xpv;
pub use vml::pvvec::pvvector_ops::pvvector_plus_pvvector as pvppv;
pub use vml::pvvec::pvvector_ops::pvvector_update as pvu;
pub use vml::pvvec::pvvector_ops::pvvector_update_discard_velocity as pvup;

// Matrix-vector products
pub use vml::pvvec::matrix_vector_products::rmatrix_multiply_pvvector as rxpv;
pub use vml::pvvec::matrix_vector_products::rmatrix_multiply_pvvector_transpose as trxpv;
