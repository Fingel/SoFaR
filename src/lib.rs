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
