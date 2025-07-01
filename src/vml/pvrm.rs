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

///  Build rotations
/// - RX        rotate r-matrix about x
/// - RY        rotate r-matrix about y
/// - RZ        rotate r-matrix about z
pub mod rotations {

    ///  Rotate an r-matrix about the x-axis.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     phi    double          angle (radians)
    ///
    ///  Given and returned:
    ///     r      double[3][3]    r-matrix, rotated
    ///
    ///  Notes:
    ///
    ///  1) Calling this function with positive phi incorporates in the
    ///     supplied r-matrix r an additional rotation, about the x-axis,
    ///     anticlockwise as seen looking towards the origin from positive x.
    ///
    ///  2) The additional rotation can be represented by this matrix:
    ///
    ///(  1        0            0      )
    ///(                               )
    ///(  0   + cos(phi)   + sin(phi)  )
    ///(                               )
    ///(  0   - sin(phi)   + cos(phi)  )
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn rx(phi: f64, r: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let s = phi.sin();
        let c = phi.cos();
        let a10 = c * r[1][0] + s * r[2][0];
        let a11 = c * r[1][1] + s * r[2][1];
        let a12 = c * r[1][2] + s * r[2][2];
        let a20 = -s * r[1][0] + c * r[2][0];
        let a21 = -s * r[1][1] + c * r[2][1];
        let a22 = -s * r[1][2] + c * r[2][2];
        [
            [r[0][0], r[0][1], r[0][2]],
            [a10, a11, a12],
            [a20, a21, a22],
        ]
    }

    ///  Rotate an r-matrix about the y-axis.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     theta  double          angle (radians)
    ///
    ///  Given and returned:
    ///     r      double[3][3]    r-matrix, rotated
    ///
    ///  Notes:
    ///
    ///  1) Calling this function with positive theta incorporates in the
    ///     supplied r-matrix r an additional rotation, about the y-axis,
    ///     anticlockwise as seen looking towards the origin from positive y.
    ///
    ///  2) The additional rotation can be represented by this matrix:
    ///
    ///(  + cos(theta)     0      - sin(theta)  )
    ///(                                        )
    ///(       0           1           0        )
    ///(                                        )
    ///(  + sin(theta)     0      + cos(theta)  )
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn ry(theta: f64, r: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let s = theta.sin();
        let c = theta.cos();
        let a00 = c * r[0][0] - s * r[2][0];
        let a01 = c * r[0][1] - s * r[2][1];
        let a02 = c * r[0][2] - s * r[2][2];
        let a20 = s * r[0][0] + c * r[2][0];
        let a21 = s * r[0][1] + c * r[2][1];
        let a22 = s * r[0][2] + c * r[2][2];

        [
            [a00, a01, a02],
            [r[1][0], r[1][1], r[1][2]],
            [a20, a21, a22],
        ]
    }

    ///  Rotate an r-matrix about the z-axis.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     psi    double          angle (radians)
    ///
    ///  Given and returned:
    ///     r      double[3][3]    r-matrix, rotated
    ///
    ///  Notes:
    ///
    ///  1) Calling this function with positive psi incorporates in the
    ///     supplied r-matrix r an additional rotation, about the z-axis,
    ///     anticlockwise as seen looking towards the origin from positive z.
    ///
    ///  2) The additional rotation can be represented by this matrix:
    ///
    ///(  + cos(psi)   + sin(psi)     0  )
    ///(                                 )
    ///(  - sin(psi)   + cos(psi)     0  )
    ///(                                 )
    ///(       0            0         1  )
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn rz(psi: f64, r: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let s = psi.sin();
        let c = psi.cos();
        let a00 = c * r[0][0] + s * r[1][0];
        let a01 = c * r[0][1] + s * r[1][1];
        let a02 = c * r[0][2] + s * r[1][2];
        let a10 = -s * r[0][0] + c * r[1][0];
        let a11 = -s * r[0][1] + c * r[1][1];
        let a12 = -s * r[0][2] + c * r[1][2];

        [
            [a00, a01, a02],
            [a10, a11, a12],
            [r[2][0], r[2][1], r[2][2]],
        ]
    }

    #[allow(clippy::excessive_precision)]
    #[cfg(test)]
    mod tests {
        use assert_approx_eq::assert_approx_eq;

        use super::*;

        /// t_sofa.c t_rx
        #[test]
        fn test_rx() {
            //TODO: would be amazing to have a Vec3 type with approx equal trait...
            let phi = 0.3456789;
            let r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let result = rx(phi, r);
            assert_eq!(result[0][0], 2.0);
            assert_eq!(result[0][1], 3.0);
            assert_eq!(result[0][2], 2.0);

            assert_approx_eq!(result[1][0], 3.839043388235612460, 1e-12);
            assert_approx_eq!(result[1][1], 3.237033249594111899, 1e-12);
            assert_approx_eq!(result[1][2], 4.516714379005982719, 1e-12);

            assert_approx_eq!(result[2][0], 1.806030415924501684, 1e-12);
            assert_approx_eq!(result[2][1], 3.085711545336372503, 1e-12);
            assert_approx_eq!(result[2][2], 3.687721683977873065, 1e-12);
        }

        /// t_sofa.c t_ry
        #[test]
        fn test_ry() {
            let theta = 0.3456789;
            let r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let result = ry(theta, r);
            assert_approx_eq!(result[0][0], 0.8651847818978159930, 1e-12);
            assert_approx_eq!(result[0][1], 1.467194920539316554, 1e-12);
            assert_approx_eq!(result[0][2], 0.1875137911274457342, 1e-12);

            assert_approx_eq!(result[1][0], 3.0, 1e-12);
            assert_approx_eq!(result[1][1], 2.0, 1e-12);
            assert_approx_eq!(result[1][2], 3.0, 1e-12);

            assert_approx_eq!(result[2][0], 3.500207892850427330, 1e-12);
            assert_approx_eq!(result[2][1], 4.779889022262298150, 1e-12);
            assert_approx_eq!(result[2][2], 5.381899160903798712, 1e-12);
        }

        /// t_sofa.c t_rz
        #[test]
        fn test_rz() {
            let psi = 0.3456789;
            let r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let result = rz(psi, r);
            assert_approx_eq!(result[0][0], 2.898197754208926769, 1e-12);
            assert_approx_eq!(result[0][1], 3.500207892850427330, 1e-12);
            assert_approx_eq!(result[0][2], 2.898197754208926769, 1e-12);

            assert_approx_eq!(result[1][0], 2.144865911309686813, 1e-12);
            assert_approx_eq!(result[1][1], 0.865184781897815993, 1e-12);
            assert_approx_eq!(result[1][2], 2.144865911309686813, 1e-12);

            assert_approx_eq!(result[2][0], 3.0, 1e-12);
            assert_approx_eq!(result[2][1], 4.0, 1e-12);
            assert_approx_eq!(result[2][2], 5.0, 1e-12);
        }
    }
}
