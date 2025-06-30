pub mod constants;
pub mod vectormatrix;

/// Round to nearest whole number (double)
/// This is equivalent to Rust `x.round()`
pub fn dnint(x: f64) -> f64 {
    if x.abs() < 0.5 {
        0.0
    } else if x < 0.0 {
        (x - 0.5).ceil()
    } else {
        (x + 0.5).floor()
    }
}

/// Truncate to the nearest whole number towards zero (double)
/// This is equivalent to Rust `x.trunc()`
pub fn dint(x: f64) -> f64 {
    if x < 0.0 { x.ceil() } else { x.floor() }
}
