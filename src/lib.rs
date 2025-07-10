use std::fmt;

pub mod astrometry;
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

// Astrometry
pub use crate::astrometry::{
    apply_pm_px as pmpx, apply_pm_safe as pmsafe, astro_icrs_to_cirs_quick as atciqz,
    catalog_to_astrometric as atcc13, catalog_to_astrometric_quick as atccq,
    catalog_to_cirs as atci13, cirs_to_icrs as atic13, cirs_to_icrs_quick as aticq,
    cirs_to_icrs_quick_def as aticqn, cirs_to_observed as atio13, cirs_to_observed_quick as atioq,
    earth_rotation as aper, earth_rotation_13 as aper13, icrs_to_cirs_quick as atciq,
    icrs_to_cirs_quick_def as atciqn, icrs_to_observed as atco13, light_deflection as ld,
    light_deflection_n as ldn, light_delfection_sun as ldsun, observatory_pv as pvtob,
    observed_to_astro_icrs as atoc13, observed_to_cirs as atoi13, observed_to_cirs_quick as atoiq,
    prep_c2o_terrestrial as apio, prep_c2o_terrestrial_13 as apio13, prep_i2c_space as apcs,
    prep_i2c_space_13 as apcs13, prep_i2c_terrestrial as apci, prep_i2c_terrestrial_13 as apci13,
    prep_i2g_geocentric as apcg, prep_i2g_geocentric_13 as apcg13, prep_i2obs_terrestrial as apco,
    prep_i2obs_terrestrial_13 as apco13, pv_to_star_catalog as pvstar,
    refraction_constants as refco, star_pm as starpm, star_pv as starpv, stellar_aberration as ab,
};

// Vector Math

// Angles
// Wrap
pub use vml::angles::wrap::{angle_normalize_pm as anpm, angle_normalize_positive as anp};

// Angle conversions
pub use vml::angles::to_sexagesimal::{
    days_to_hms as d2tf, radians_to_dms as a2af, radians_to_hms as a2tf,
};

// Pvectors and Rmatrices
// Initialize
pub use vml::pvrm::initialize::{identity_rmatrix as ir, zero_pvector as zp, zero_rmatrix as zr};

// Rotations
pub use vml::pvrm::rotations::{
    rotate_rmatrix_about_x as rx, rotate_rmatrix_about_y as ry, rotate_rmatrix_about_z as rz,
};

// Spherical/Cartesian conversions
pub use vml::pvrm::sphere_cart_conv::{
    p_vector_to_spherical as p2s, spherical_to_pvector as s2p, spherical_to_unit_vector as s2c,
    unit_vector_to_spherical as c2s,
};

// Operations on vectors
pub use vml::pvrm::vec_ops::{
    pvector_cross_product as pxp, pvector_dot_product as pdp, pvector_minus_pvector as pmp,
    pvector_modulus as pm, pvector_multiply_scalar as sxp, pvector_normalize as pn,
    pvector_plus_pvector as ppp, pvector_plus_scaled_pvector as ppsp,
};

// Operations on matrices
pub use vml::pvrm::matrix_ops::{rmatrix_multiply as rxr, transpose_rmatrix as tr};

// Matrix-vector products
pub use vml::pvrm::matrix_vec_products::{
    rmatrix_pvector_product as rxp, transpose_rmatrix_pvector_product as trxp,
};

// Separation and position-angle
pub use vml::pvrm::sep_position_angle::{
    angular_separation_pvector as sepp, angular_separation_spherical as seps,
    position_angle_from_pvector as pap, position_angle_from_spherical as pas,
};

// Rotation vectors
pub use vml::pvrm::rotation_vectors::{rmatrix_to_rvector as rm2v, rvector_to_rmatrix as rv2m};

// PV Vectors
// Initialize
pub use vml::pvvec::initialize::zero_pvvector as zpv;

// Copy, extend, extract
pub use vml::pvvec::copy_extend_extract::{
    append_zvelocity_pvvector as p2pv, copy_pvvector as cpv, discard_velocity_pvvector as pv2p,
};

// Sphereical/Cartesian conversions
pub use vml::pvvec::sphere_cart_conv::{
    pvvector_to_spherical as pv2s, spherical_to_pvvector as s2pv,
};

// Operations on pv-vectors
pub use vml::pvvec::pvvector_ops::{
    pvvector_cross_pvvector as pvxpv, pvvector_dot_pvvector as pvdpv,
    pvvector_minus_pvvector as pvmpv, pvvector_modulus as pvm, pvvector_multiply_scalar as sxpv,
    pvvector_multiply_two_scalar as s2xpv, pvvector_plus_pvvector as pvppv, pvvector_update as pvu,
    pvvector_update_discard_velocity as pvup,
};

// Matrix-vector products
pub use vml::pvvec::matrix_vector_products::{
    rmatrix_multiply_pvvector as rxpv, rmatrix_multiply_pvvector_transpose as trxpv,
};
