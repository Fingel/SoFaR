//! Operations involving p-vectors and r-matrices

//   Initialize
//   ZP        zero p-vector
//   ZR        initialize r-matrix to null
//   IR        initialize r-matrix to identity

///  - - - - - -
///   i a u Z p
///  - - - - - -
///
///  Zero a p-vector.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  vector/matrix support function.
///
///  Returned:
///     p        double[3]      zero p-vector
///
///  This revision:  2021 May 11
///
///  SOFA release 2023-10-11
///
///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
pub fn zp() -> [f64; 3] {
    [0.0, 0.0, 0.0]
}
pub use zp as zero_p_vector;
///  - - - - - -
///   i a u Z r
///  - - - - - -
///
///  Initialize an r-matrix to the null matrix.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  vector/matrix support function.
///
///  Returned:
///     r        double[3][3]    r-matrix
///
///  This revision:  2021 May 11
///
///  SOFA release 2023-10-11
///
///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
pub fn zr() -> [[f64; 3]; 3] {
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
}
pub use zr as zero_r_matrix;

///  - - - - - -
///   i a u I r
///  - - - - - -
///
///  Initialize an r-matrix to the identity matrix.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  vector/matrix support function.
///
///  Returned:
///     r       double[3][3]    r-matrix
///
///  This revision:  2021 May 11
///
///  SOFA release 2023-10-11
///
///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
pub fn ir() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
pub use ir as identity_r_matrix;

#[cfg(test)]
mod tests {
    use super::*;

    /// t_sofa.c t_zp
    #[test]
    #[allow(unused_assignments)]
    fn test_zp() {
        let mut p = [0.3, 1.2, -2.5];
        p = zp();
        assert_eq!(p, [0.0, 0.0, 0.0]);
    }

    /// t_sofa.c t_zr
    #[test]
    #[allow(unused_assignments)]
    fn test_zr() {
        let mut r = [[2.0, 3.0, 3.0], [3.0, 2.0, 4.0], [2.0, 3.0, 5.0]];
        r = zr();
        assert_eq!(r, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
    }

    /// t_sofa.c t_ir
    #[test]
    #[allow(unused_assignments)]
    fn test_ir() {
        let mut r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
        r = ir();
        assert_eq!(r, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    }
}
