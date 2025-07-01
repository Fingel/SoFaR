//! Operations involving p-vectors and r-matrices

///   Initialize
///   - ZP        zero p-vector
///   - ZR        initialize r-matrix to null
///   - IR        initialize r-matrix to identity
pub mod initialize {
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
}

///   Copy
///   - CP        copy p-vector
///   - CR        copy r-matrix
pub mod copy {
    ///  Copy a p-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     p        double[3]     p-vector to be copied
    ///
    ///  Returned:
    ///     c        double[3]     copy
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn cp(p: &[f64; 3]) -> [f64; 3] {
        //TODO: This is pointless as arrays implement Copy trait
        [p[0], p[1], p[2]]
    }

    ///  Copy an r-matrix.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     r        double[3][3]    r-matrix to be copied
    ///
    ///  Returned:
    ///     c        double[3][3]    copy
    ///
    ///  Called:
    ///     iauCp        copy p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn cr(r: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        //TODO: This is pointless as arrays implement Copy trait
        [
            [r[0][0], r[0][1], r[0][2]],
            [r[1][0], r[1][1], r[1][2]],
            [r[2][0], r[2][1], r[2][2]],
        ]
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// t_sofa.c t_cp
        #[test]
        fn test_cp() {
            let p = [0.3, 1.2, -2.5];
            let c = cp(&p);
            //TODO let c = p;
            assert_eq!(c, p);
        }

        /// t_sofa.c t_cr
        #[test]
        fn test_cr() {
            let r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let c = cr(&r);
            //TODO let c = r;
            assert_eq!(c, r);
        }
    }
}
