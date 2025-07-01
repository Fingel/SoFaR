pub mod constants;
pub mod vml;

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
