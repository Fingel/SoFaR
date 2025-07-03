//! Operations involving p-vectors and r-matrices

///   Initialize
///   - ZP        zero p-vector
///   - ZR        initialize r-matrix to null
///   - IR        initialize r-matrix to identity
pub mod initialize {
    use crate::{Pvector, Rmatrix};

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
    pub fn zp() -> Pvector {
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
    pub fn zr() -> Rmatrix {
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
    pub fn ir() -> Rmatrix {
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
    use crate::{Pvector, Rmatrix};
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
    pub fn cp(p: &Pvector) -> Pvector {
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
    pub fn cr(r: &Rmatrix) -> Rmatrix {
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
    use crate::Rmatrix;

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
    pub fn rx(phi: f64, r: &Rmatrix) -> Rmatrix {
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
    pub fn ry(theta: f64, r: &Rmatrix) -> Rmatrix {
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
    pub fn rz(psi: f64, r: &Rmatrix) -> Rmatrix {
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
            let result = rx(phi, &r);
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
            let result = ry(theta, &r);
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
            let result = rz(psi, &r);
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

/// Spherical/Cartesian conversions
/// - S2C       spherical to unit vector
/// - C2S       unit vector to spherical
/// - S2P       spherical to p-vector
/// - P2S       p-vector to spherical
pub mod sphere_cart_conv {
    use crate::Pvector;
    use crate::vml::pvrm::vec_ops::{pm, sxp};

    ///  Convert spherical coordinates to Cartesian.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     theta    double       longitude angle (radians)
    ///     phi      double       latitude angle (radians)
    ///
    ///  Returned:
    ///     c        double[3]    direction cosines
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn s2c(theta: f64, phi: f64) -> Pvector {
        let cp = phi.cos();
        [theta.cos() * cp, theta.sin() * cp, phi.sin()]
    }
    pub use s2c as spherical_to_unit_vector;

    ///  P-vector to spherical coordinates.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     p      double[3]    p-vector
    ///
    ///  Returned:
    ///     theta  double       longitude angle (radians)
    ///     phi    double       latitude angle (radians)
    ///
    ///  Notes:
    ///
    ///  1) The vector p can have any magnitude; only its direction is used.
    ///
    ///  2) If p is null, zero theta and phi are returned.
    ///
    ///  3) At either pole, zero theta is returned.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn c2s(p: &Pvector) -> (f64, f64) {
        let x = p[0];
        let y = p[1];
        let z = p[2];
        let d2 = x * x + y * y;
        let theta = if d2 == 0.0 { 0.0 } else { y.atan2(x) };
        let phi = if z == 0.0 { 0.0 } else { z.atan2(d2.sqrt()) };
        (theta, phi)
    }
    pub use c2s as unit_vector_to_spherical;

    ///  Convert spherical polar coordinates to p-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     theta   double       longitude angle (radians)
    ///     phi     double       latitude angle (radians)
    ///     r       double       radial distance
    ///
    ///  Returned:
    ///     p       double[3]    Cartesian coordinates
    ///
    ///  Called:
    ///     iauS2c       spherical coordinates to unit vector
    ///     iauSxp       multiply p-vector by scalar
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn s2p(theta: f64, phi: f64, r: f64) -> Pvector {
        let u = s2c(theta, phi);
        sxp(r, &u)
    }
    pub use s2p as spherical_to_p_vector;

    ///  P-vector to spherical polar coordinates.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     p        double[3]    p-vector
    ///
    ///  Returned:
    ///     theta    double       longitude angle (radians)
    ///     phi      double       latitude angle (radians)
    ///     r        double       radial distance
    ///
    ///  Notes:
    ///
    ///  1) If P is null, zero theta, phi and r are returned.
    ///
    ///  2) At either pole, zero theta is returned.
    ///
    ///  Called:
    ///     iauC2s       p-vector to spherical
    ///     iauPm        modulus of p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn p2s(p: &Pvector) -> (f64, f64, f64) {
        let (theta, phi) = c2s(p);
        let r = pm(p);
        (theta, phi, r)
    }
    pub use p2s as p_vector_to_spherical;

    #[cfg(test)]
    mod tests {
        use assert_approx_eq::assert_approx_eq;

        use super::*;

        /// t_sofa.c t_s2c
        #[test]
        fn test_s2c() {
            let c = s2c(3.0123, -0.999);
            assert_approx_eq!(c[0], -0.5366267667260523906, 1e-12);
            assert_approx_eq!(c[1], 0.0697711109765145365, 1e-12);
            assert_approx_eq!(c[2], -0.8409302618566214041, 1e-12);
        }

        /// t_sofa.c t_c2s
        #[test]
        fn test_c2s() {
            let p = [100.0, -50.0, 25.0];
            let (theta, phi) = c2s(&p);
            assert_approx_eq!(theta, -0.4636476090008061162, 1e-14);
            assert_approx_eq!(phi, 0.2199879773954594463, 1e-14);
        }

        /// t_sofa.c t_s2p
        #[test]
        fn test_s2p() {
            let p = s2p(-3.21, 0.123, 0.456);
            assert_approx_eq!(p[0], -0.4514964673880165228, 1e-12);
            assert_approx_eq!(p[1], 0.0309339427734258688, 1e-12);
            assert_approx_eq!(p[2], 0.0559466810510877933, 1e-12);
        }

        /// t_sofa.c t_p2s
        #[test]
        fn test_p2s() {
            let p = [100.0, -50.0, 25.0];
            let (theta, phi, r) = p2s(&p);
            assert_approx_eq!(theta, -0.4636476090008061162, 1e-12);
            assert_approx_eq!(phi, 0.2199879773954594463, 1e-12);
            assert_approx_eq!(r, 114.5643923738960002, 1e-9);
        }
    }
}

///Operations on vectors
///- PPP       p-vector plus p-vector
///- PMP       p-vector minus p-vector
///- PPSP      p-vector plus scaled p-vector
///- PDP       inner (=scalar=dot) product of two p-vectors
///- PXP       outer (=vector=cross) product of two p-vectors
///- PM        modulus of p-vector
///- PN        normalize p-vector returning modulus
///- SXP       multiply p-vector by scalar
pub mod vec_ops {
    use crate::Pvector;
    use crate::vml::pvrm::initialize::zp;
    //TODO: SIMD would probably provide gains for this module

    ///  P-vector addition.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double[3]      first p-vector
    ///     b        double[3]      second p-vector
    ///
    ///  Returned:
    ///     apb      double[3]      a + b
    ///
    ///  Note:
    ///     It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn ppp(a: &Pvector, b: &Pvector) -> Pvector {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }
    pub use ppp as pvector_plus_pvector;

    ///  P-vector subtraction.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double[3]      first p-vector
    ///     b        double[3]      second p-vector
    ///
    ///  Returned:
    ///     amb      double[3]      a - b
    ///
    ///  Note:
    ///     It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn pmp(a: &Pvector, b: &Pvector) -> Pvector {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }
    pub use pmp as pvector_minus_pvector;

    ///  P-vector plus scaled p-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a      double[3]     first p-vector
    ///     s      double        scalar (multiplier for b)
    ///     b      double[3]     second p-vector
    ///
    ///  Returned:
    ///     apsb   double[3]     a + s*b
    ///
    ///  Note:
    ///     It is permissible for any of a, b and apsb to be the same array.
    ///
    ///  Called:
    ///     iauSxp       multiply p-vector by scalar
    ///     iauPpp       p-vector plus p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end
    pub fn ppsp(a: &Pvector, s: f64, b: &Pvector) -> Pvector {
        let sb = sxp(s, b);
        ppp(a, &sb)
    }
    pub use ppsp as pvector_plus_scaled_pvector;

    ///  p-vector inner (=scalar=dot) product.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a      double[3]     first p-vector
    ///     b      double[3]     second p-vector
    ///
    ///  Returned (function value):
    ///            double        a . b
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn pdp(a: &Pvector, b: &Pvector) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
    pub use pdp as pvector_dot_product;

    ///  p-vector outer (=vector=cross) product.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double[3]      first p-vector
    ///     b        double[3]      second p-vector
    ///
    ///  Returned:
    ///     axb      double[3]      a x b
    ///
    ///  Note:
    ///     It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn pxp(a: &Pvector, b: &Pvector) -> Pvector {
        let xa = a[0];
        let ya = a[1];
        let za = a[2];
        let xb = b[0];
        let yb = b[1];
        let zb = b[2];
        [ya * zb - za * yb, za * xb - xa * zb, xa * yb - ya * xb]
    }
    pub use pxp as pvector_cross_product;

    ///  Modulus of p-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     p      double[3]     p-vector
    ///
    ///  Returned (function value):
    ///            double        modulus
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn pm(p: &Pvector) -> f64 {
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
    }
    pub use pm as pvector_modulus;

    ///  Convert a p-vector into modulus and unit vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     p        double[3]      p-vector
    ///
    ///  Returned:
    ///     r        double         modulus
    ///     u        double[3]      unit vector
    ///
    ///  Notes:
    ///
    ///  1) If p is null, the result is null.  Otherwise the result is a unit
    ///     vector.
    ///
    ///  2) It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  Called:
    ///     iauPm        modulus of p-vector
    ///     iauZp        zero p-vector
    ///     iauSxp       multiply p-vector by scalar
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn pn(p: &Pvector) -> (f64, Pvector) {
        // Obtain the modulus and test for zero.
        let w = pm(p);
        if w == 0.0 {
            // Null vector.
            let u = zp();
            (w, u)
        } else {
            // Unit vector.
            let u = sxp(1.0 / w, p);
            (w, u)
        }
    }
    pub use pn as pvector_normalize;

    ///  Multiply a p-vector by a scalar.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     s      double        scalar
    ///     p      double[3]     p-vector
    ///
    ///  Returned:
    ///     sp     double[3]     s * p
    ///
    ///  Note:
    ///     It is permissible for p and sp to be the same array.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn sxp(s: f64, p: &Pvector) -> Pvector {
        [s * p[0], s * p[1], s * p[2]]
    }
    pub use sxp as pvector_multiply_scalar;

    #[cfg(test)]
    mod tests {
        use assert_approx_eq::assert_approx_eq;

        use super::*;

        /// t_sofa.c t_ppp
        #[test]
        fn test_ppp() {
            let a = [2.0, 2.0, 3.0];
            let b = [1.0, 3.0, 4.0];
            let apb = ppp(&a, &b);
            assert_eq!(apb, [3.0, 5.0, 7.0]);
        }

        /// t_sofa.c t_pmp
        #[test]
        fn test_pmp() {
            let a = [2.0, 2.0, 3.0];
            let b = [1.0, 3.0, 4.0];
            let amb = pmp(&a, &b);
            assert_eq!(amb, [1.0, -1.0, -1.0]);
        }

        /// t_sofa.c t_ppsp
        #[test]
        fn test_ppsp() {
            let a = [2.0, 2.0, 3.0];
            let s = 5.0;
            let b = [1.0, 3.0, 4.0];
            let apsb = ppsp(&a, s, &b);
            assert_eq!(apsb, [7.0, 17.0, 23.0]);
        }

        /// t_sofa.c t_pdp
        #[test]
        fn test_pdp() {
            let a = [2.0, 2.0, 3.0];
            let b = [1.0, 3.0, 4.0];
            let apb = pdp(&a, &b);
            assert_eq!(apb, 20.0);
        }

        /// t_sofa.c t_pxp
        #[test]
        fn test_pxp() {
            let a = [2.0, 2.0, 3.0];
            let b = [1.0, 3.0, 4.0];
            let axb = pxp(&a, &b);
            assert_eq!(axb, [-1.0, -5.0, 4.0]);
        }

        /// t_sofa.c t_pm
        #[test]
        fn test_pm() {
            let p = [0.3, 1.2, -2.5];
            let r = pm(&p);
            assert_approx_eq!(r, 2.789265136196270604, 1e-12);
        }

        /// t_sofa.c t_pn
        #[test]
        fn test_pn() {
            let p = [0.3, 1.2, -2.5];
            let (r, u) = pn(&p);
            assert_approx_eq!(r, 2.789265136196270604, 1e-12);

            assert_approx_eq!(u[0], 0.1075552109073112058, 1e-12);
            assert_approx_eq!(u[1], 0.4302208436292448232, 1e-12);
            assert_approx_eq!(u[2], -0.8962934242275933816, 1e-12);
        }

        #[test]
        fn test_pn_zero() {
            let p = [0.0, 0.0, 0.0];
            let (r, u) = pn(&p);
            assert_eq!(r, 0.0);
            assert_eq!(u, [0.0, 0.0, 0.0]);
        }

        /// t_sofa.c t_sxp
        #[test]
        fn test_sxp() {
            let s = 2.0;
            let p = [0.3, 1.2, -2.5];
            let sp = sxp(s, &p);
            assert_eq!(sp, [0.6, 2.4, -5.0]);
        }
    }
}

/// Operations on matrices
///
/// RXR       r-matrix multiply
/// TR        transpose r-matrix
pub mod matrix_ops {
    use crate::Rmatrix;
    use crate::vml::pvrm::initialize::zr;

    // TODO: SIMD
    ///  Multiply two r-matrices.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double[3][3]    first r-matrix
    ///     b        double[3][3]    second r-matrix
    ///
    ///  Returned:
    ///     atb      double[3][3]    a * b
    ///
    ///  Note:
    ///     It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  Called:
    ///     iauCr        copy r-matrix
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    #[allow(clippy::needless_range_loop)]
    pub fn rxr(a: &Rmatrix, b: &Rmatrix) -> Rmatrix {
        //TODO: naive mmultiply implementation
        let mut atb = zr();
        let mut w;
        for i in 0..3 {
            for j in 0..3 {
                w = 0.0;
                for k in 0..3 {
                    w += a[i][k] * b[k][j];
                }
                atb[i][j] = w;
            }
        }
        atb
    }
    pub use rxr as rmatrix_mutiply;

    ///  Transpose an r-matrix.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     r        double[3][3]    r-matrix
    ///
    ///  Returned:
    ///     rt       double[3][3]    transpose
    ///
    ///  Note:
    ///     It is permissible for r and rt to be the same array.
    ///
    ///  Called:
    ///     iauCr        copy r-matrix
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    #[allow(clippy::needless_range_loop)]
    pub fn tr(r: &Rmatrix) -> Rmatrix {
        let mut tr = zr();
        for i in 0..3 {
            for j in 0..3 {
                tr[i][j] = r[j][i];
            }
        }
        tr
    }
    pub use tr as transpose_rmatrix;

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::vml::pvrm::initialize::zr;

        /// t_sofa.c t_rxr
        #[test]
        fn test_rxr() {
            let a = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let b = [[1.0, 2.0, 2.0], [4.0, 1.0, 1.0], [3.0, 0.0, 1.0]];
            let atb = [[20.0, 7.0, 9.0], [20.0, 8.0, 11.0], [34.0, 10.0, 15.0]];
            let result = rxr(&a, &b);
            assert_eq!(result, atb);
        }

        #[test]
        fn test_rxr_parity() {
            use rsofa::iauRxr;
            let mut a = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let mut b = [[1.0, 2.0, 2.0], [4.0, 1.0, 1.0], [3.0, 0.0, 1.0]];
            let atb = rxr(&a, &b);
            let mut atb_iau = zr();
            unsafe {
                iauRxr(a.as_mut_ptr(), b.as_mut_ptr(), atb_iau.as_mut_ptr());
            }
            assert_eq!(atb, atb_iau);
        }

        /// t_sofa.c t_tr
        #[test]
        fn test_tr() {
            let r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let rt = [[2.0, 3.0, 3.0], [3.0, 2.0, 4.0], [2.0, 3.0, 5.0]];
            let result = tr(&r);
            assert_eq!(result, rt);
        }

        #[test]
        fn test_tr_parity() {
            use rsofa::iauTr;
            let mut r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            // let rt = [[2.0, 3.0, 3.0], [3.0, 2.0, 4.0], [2.0, 3.0, 5.0]];
            let rt = tr(&r);
            let mut rt_iau = zr();
            unsafe {
                iauTr(r.as_mut_ptr(), rt_iau.as_mut_ptr());
            }
            assert_eq!(rt, rt_iau);
        }
    }
}

/// Matrix-vector products
/// - RXP       product of r-matrix and p-vector
/// - TRXP      product of transpose of r-matrix and p-vector
pub mod matrix_vec_products {
    use crate::vml::pvrm::initialize::zp;
    use crate::vml::pvrm::matrix_ops::tr;
    use crate::{Pvector, Rmatrix};

    ///  Multiply a p-vector by an r-matrix.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     r        double[3][3]    r-matrix
    ///     p        double[3]       p-vector
    ///
    ///  Returned:
    ///     rp       double[3]       r * p
    ///
    ///  Note:
    ///     It is permissible for p and rp to be the same array.
    ///
    ///  Called:
    ///     iauCp        copy p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    #[allow(clippy::needless_range_loop)]
    pub fn rxp(r: &Rmatrix, p: &Pvector) -> Pvector {
        //TODO No need for mut here
        let mut w;
        let mut rp = zp();
        for j in 0..3 {
            w = 0.0;
            for i in 0..3 {
                w += r[j][i] * p[i];
            }
            rp[j] = w;
        }
        rp
    }
    pub use rxp as rmatrix_pvector_product;

    ///  Multiply a p-vector by the transpose of an r-matrix.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     r        double[3][3]   r-matrix
    ///     p        double[3]      p-vector
    ///
    ///  Returned:
    ///     trp      double[3]      r^T * p
    ///
    ///  Note:
    ///     It is permissible for p and trp to be the same array.
    ///
    ///  Called:
    ///     iauTr        transpose r-matrix
    ///     iauRxp       product of r-matrix and p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn trxp(r: &Rmatrix, p: &Pvector) -> Pvector {
        // Transpose of matrix r.
        let transpose = tr(r);

        // Matrix tr * vector p -> vector trp
        rxp(&transpose, p)
    }
    pub use trxp as transpose_rmatrix_pvector_product;

    #[cfg(test)]
    mod tests {
        use assert_approx_eq::assert_approx_eq;

        use super::*;

        /// t_sofa.c t_rxp
        #[test]
        fn test_rxp() {
            let r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let p = [0.2, 1.5, 0.1];
            let rp = rxp(&r, &p);
            assert_approx_eq!(rp[0], 5.1, 1e-12);
            assert_approx_eq!(rp[1], 3.9, 1e-12);
            assert_approx_eq!(rp[2], 7.1, 1e-12);
        }

        #[test]
        fn test_rxp_parity() {
            use rsofa::iauRxp;
            let mut r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let mut p = [0.2, 1.5, 0.1];
            let rp = rxp(&r, &p);

            let mut rp_iau = zp();
            unsafe {
                iauRxp(r.as_mut_ptr(), p.as_mut_ptr(), rp_iau.as_mut_ptr());
            }
            assert_eq!(rp, rp_iau);
        }

        /// t_sofa.c t_trxp
        #[test]
        fn test_trxp() {
            let r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let p = [0.2, 1.5, 0.1];
            let trp = trxp(&r, &p);
            assert_approx_eq!(trp[0], 5.2, 1e-12);
            assert_approx_eq!(trp[1], 4.0, 1e-12);
            assert_approx_eq!(trp[2], 5.4, 1e-12);
        }

        #[test]
        fn test_trxp_parity() {
            use rsofa::iauTrxp;
            let mut r = [[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
            let mut p = [0.2, 1.5, 0.1];
            let trp = trxp(&r, &p);
            let mut trp_iau = zp();
            unsafe {
                iauTrxp(r.as_mut_ptr(), p.as_mut_ptr(), trp_iau.as_mut_ptr());
            }
            assert_eq!(trp, trp_iau);
        }
    }
}

/// Separation and position-angle
/// - SEPP      angular separation from p-vectors
/// - SEPS      angular separation from spherical coordinates
/// - PAP       position-angle from p-vectors
/// - PAS       position-angle from spherical coordinates
pub mod sep_position_angle {
    use crate::{
        Pvector,
        vml::pvrm::{
            initialize::zp,
            sphere_cart_conv::s2c,
            vec_ops::{pdp, pm, pmp, pn, pxp},
        },
    };

    ///  Angular separation between two p-vectors.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a      double[3]    first p-vector (not necessarily unit length)
    ///     b      double[3]    second p-vector (not necessarily unit length)
    ///
    ///  Returned (function value):
    ///            double       angular separation (radians, always positive)
    ///
    ///  Notes:
    ///
    ///  1) If either vector is null, a zero result is returned.
    ///
    ///  2) The angular separation is most simply formulated in terms of
    ///     scalar product.  However, this gives poor accuracy for angles
    ///     near zero and pi.  The present algorithm uses both cross product
    ///     and dot product, to deliver full accuracy whatever the size of
    ///     the angle.
    ///
    ///  Called:
    ///     iauPxp       vector product of two p-vectors
    ///     iauPm        modulus of p-vector
    ///     iauPdp       scalar product of two p-vectors
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn sepp(a: &Pvector, b: &Pvector) -> f64 {
        // Sine of angle between the vectors, multiplied by the two moduli.
        let axb = pxp(a, b);
        let ss = pm(&axb);

        // Cosine of the angle, multiplied by the two moduli.
        let cs = pdp(a, b);

        // The angle.
        if ss != 0.0 || cs != 0.0 {
            ss.atan2(cs)
        } else {
            0.0
        }
    }
    pub use sepp as angular_separation_pvector;

    ///  Angular separation between two sets of spherical coordinates.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     al     double       first longitude (radians)
    ///     ap     double       first latitude (radians)
    ///     bl     double       second longitude (radians)
    ///     bp     double       second latitude (radians)
    ///
    ///  Returned (function value):
    ///            double       angular separation (radians)
    ///
    ///  Called:
    ///     iauS2c       spherical coordinates to unit vector
    ///     iauSepp      angular separation between two p-vectors
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn seps(al: f64, ap: f64, bl: f64, bp: f64) -> f64 {
        // Spherical to Cartesian.
        let ac = s2c(al, ap);
        let bc = s2c(bl, bp);

        // Angle between the vectors.
        sepp(&ac, &bc)
    }
    pub use seps as angular_separation_spherical;

    ///  Position-angle from two p-vectors.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a      double[3]  direction of reference point
    ///     b      double[3]  direction of point whose PA is required
    ///
    ///  Returned (function value):
    ///            double     position angle of b with respect to a (radians)
    ///
    ///  Notes:
    ///
    ///  1) The result is the position angle, in radians, of direction b with
    ///     respect to direction a.  It is in the range -pi to +pi.  The
    ///     sense is such that if b is a small distance "north" of a the
    ///     position angle is approximately zero, and if b is a small
    ///     distance "east" of a the position angle is approximately +pi/2.
    ///
    ///  2) The vectors a and b need not be of unit length.
    ///
    ///  3) Zero is returned if the two directions are the same or if either
    ///     vector is null.
    ///
    ///  4) If vector a is at a pole, the result is ill-defined.
    ///
    ///  Called:
    ///     iauPn        decompose p-vector into modulus and direction
    ///     iauPm        modulus of p-vector
    ///     iauPxp       vector product of two p-vectors
    ///     iauPmp       p-vector minus p-vector
    ///     iauPdp       scalar product of two p-vectors
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn pap(a: &Pvector, b: &Pvector) -> f64 {
        let st: f64;
        let mut ct: f64;
        let mut eta = zp();

        // Modulus and direction of the a vector.
        let (am, au) = pn(a);

        // Modulus of the b vector.
        let bm = pm(b);

        // Deal with the case of a null vector
        if am == 0.0 || bm == 0.0 {
            st = 0.0;
            ct = 1.0;
        } else {
            // The "north" axis tangential from a (arbitrary length).
            let xa = a[0];
            let ya = a[1];
            let za = a[2];
            eta[0] = -xa * za;
            eta[1] = -ya * za;
            eta[2] = xa * xa + ya * ya;

            // The "east" axis tangential from a (same length).
            let xi = pxp(&eta, &au);

            // The vector from a to b.
            let a2b = pmp(b, a);

            // Resolve into components along the north and east axes.
            st = pdp(&a2b, &xi);
            ct = pdp(&a2b, &eta);

            // Deal with degenerate cases
            if st == 0.0 && ct == 0.0 {
                ct = 1.0;
            }
        }

        // Position angle
        st.atan2(ct)
    }
    pub use pap as postition_angle_from_pvector;

    ///  Position-angle from spherical coordinates.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     al     double     longitude of point A (e.g. RA) in radians
    ///     ap     double     latitude of point A (e.g. Dec) in radians
    ///     bl     double     longitude of point B
    ///     bp     double     latitude of point B
    ///
    ///  Returned (function value):
    ///            double     position angle of B with respect to A
    ///
    ///  Notes:
    ///
    ///  1) The result is the bearing (position angle), in radians, of point
    ///     B with respect to point A.  It is in the range -pi to +pi.  The
    ///     sense is such that if B is a small distance "east" of point A,
    ///     the bearing is approximately +pi/2.
    ///
    ///  2) Zero is returned if the two points are coincident.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end
    pub fn pas(al: f64, ap: f64, bl: f64, bp: f64) -> f64 {
        let dl = bl - al;
        let y = dl.sin() * bp.cos();
        let x = bp.sin() * ap.cos() - bp.cos() * ap.sin() * dl.cos();
        if x != 0.0 || y != 0.0 {
            y.atan2(x)
        } else {
            0.0
        }
    }
    pub use pas as postition_angle_from_spherical;

    #[cfg(test)]
    mod tests {
        use assert_approx_eq::assert_approx_eq;

        use super::*;

        /// t_sofa.c t_sepp
        #[test]
        fn test_sepp() {
            let a = [1.0, 0.1, 0.2];
            let b = [-3.0, 1e-3, 0.2];
            let s = sepp(&a, &b);
            assert_approx_eq!(s, 2.860391919024660768, 1e-12);
        }

        #[test]
        fn test_sepp_parity() {
            use rsofa::iauSepp;
            let mut a = [1.0, 0.1, 0.2];
            let mut b = [-3.0, 1e-3, 0.2];
            let s = sepp(&a, &b);
            unsafe {
                let s_iau = iauSepp(a.as_mut_ptr(), b.as_mut_ptr());
                assert_eq!(s, s_iau);
            }
        }

        /// t_sofa.c t_seps
        #[test]
        fn test_seps() {
            let al = 1.0;
            let ap = 0.1;

            let bl = 0.2;
            let bp = -3.0;
            let s = seps(al, ap, bl, bp);
            assert_approx_eq!(s, 2.346722016996998842, 1e-14);
        }

        #[test]
        fn test_seps_parity() {
            use rsofa::iauSeps;
            let al = 1.0;
            let ap = 0.1;

            let bl = 0.2;
            let bp = -3.0;
            let s = seps(al, ap, bl, bp);
            unsafe {
                let s_iau = iauSeps(al, ap, bl, bp);
                assert_eq!(s, s_iau);
            }
        }

        /// t_sofa.c t_pap
        #[test]
        fn test_pap() {
            let a = [1.0, 0.1, 0.2];
            let b = [-3.0, 1e-3, 0.2];
            let theta = pap(&a, &b);
            assert_approx_eq!(theta, 0.3671514267841113674, 1e-12);
        }

        #[test]
        fn test_pap_parity() {
            use rsofa::iauPap;
            let mut a = [1.0, 0.1, 0.2];
            let mut b = [-3.0, 1e-3, 0.2];
            let theta = pap(&a, &b);
            unsafe {
                let theta_iau = iauPap(a.as_mut_ptr(), b.as_mut_ptr());
                assert_eq!(theta, theta_iau);
            }
        }

        /// t_sofa.c t_pas
        #[test]
        fn test_pas() {
            let al = 1.0;
            let ap = 0.1;
            let bl = 0.2;
            let bp = -1.0;

            let theta = pas(al, ap, bl, bp);
            assert_approx_eq!(theta, -2.724544922932270424, 1e-12);
        }

        #[test]
        fn test_pas_parity() {
            use rsofa::iauPas;
            let al = 1.0;
            let ap = 0.1;
            let bl = 0.2;
            let bp = -1.0;

            let theta = pas(al, ap, bl, bp);
            unsafe {
                let theta_iau = iauPas(al, ap, bl, bp);
                assert_eq!(theta, theta_iau);
            }
        }
    }
}

/// Rotation vectors
/// - RV2M      r-vector to r-matrix
/// - RM2V      r-matrix to r-vector
pub mod rotation_vectors {
    use crate::{Pvector, Rmatrix, vml::pvrm::initialize::zr};
    ///  Form the r-matrix corresponding to a given r-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     w        double[3]      rotation vector (Note 1)
    ///
    ///  Returned:
    ///     r        double[3][3]    rotation matrix
    ///
    ///  Notes:
    ///
    ///  1) A rotation matrix describes a rotation through some angle about
    ///     some arbitrary axis called the Euler axis.  The "rotation vector"
    ///     supplied to This function has the same direction as the Euler
    ///     axis, and its magnitude is the angle in radians.
    ///
    ///  2) If w is null, the identity matrix is returned.
    ///
    ///  3) The reference frame rotates clockwise as seen looking along the
    ///     rotation vector from the origin.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn rv2m(w: &Pvector) -> Rmatrix {
        // Euler angle (magnitude of rotation vector) and functions.
        let mut x = w[0];
        let mut y = w[1];
        let mut z = w[2];
        let phi = (x * x + y * y + z * z).sqrt();
        let s = phi.sin();
        let c = phi.cos();
        let f = 1.0 - c;

        // Euler axis (direction of rotation vector), perhaps null.
        if phi > 0.0 {
            x /= phi;
            y /= phi;
            z /= phi;
        }

        // Form the rotation matrix.
        let mut r = zr();
        r[0][0] = x * x * f + c;
        r[0][1] = x * y * f + z * s;
        r[0][2] = x * z * f - y * s;
        r[1][0] = y * x * f - z * s;
        r[1][1] = y * y * f + c;
        r[1][2] = y * z * f + x * s;
        r[2][0] = z * x * f + y * s;
        r[2][1] = z * y * f - x * s;
        r[2][2] = z * z * f + c;

        r
    }
    pub use rv2m as rvector_to_rmatrix;

    ///  Express an r-matrix as an r-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     r        double[3][3]    rotation matrix
    ///
    ///  Returned:
    ///     w        double[3]       rotation vector (Note 1)
    ///
    ///  Notes:
    ///
    ///  1) A rotation matrix describes a rotation through some angle about
    ///     some arbitrary axis called the Euler axis.  The "rotation vector"
    ///     returned by this function has the same direction as the Euler axis,
    ///     and its magnitude is the angle in radians.  (The magnitude and
    ///     direction can be separated by means of the function iauPn.)
    ///
    ///  2) If r is null, so is the result.  If r is not a rotation matrix
    ///     the result is undefined;  r must be proper (i.e. have a positive
    ///     determinant) and real orthogonal (inverse = transpose).
    ///
    ///  3) The reference frame rotates clockwise as seen looking along
    ///     the rotation vector from the origin.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end
    pub fn rm2v(r: &Rmatrix) -> Pvector {
        // TODO: Type for Rvector
        let x = r[1][2] - r[2][1];
        let y = r[2][0] - r[0][2];
        let z = r[0][1] - r[1][0];
        let s2 = (x * x + y * y + z * z).sqrt();
        if s2 > 0.0 {
            let c2 = r[0][0] + r[1][1] + r[2][2] - 1.0;
            let phi = s2.atan2(c2);
            let f = phi / s2;
            [x * f, y * f, z * f]
        } else {
            [0.0, 0.0, 0.0]
        }
    }
    pub use rm2v as rmatrix_to_rvector;

    #[cfg(test)]
    mod tests {
        use assert_approx_eq::assert_approx_eq;

        use super::*;

        /// t_sofa.c t_rv2m
        #[test]
        fn test_rv2m() {
            let w = [0.0, 1.41371669, -1.88495559];
            let r = rv2m(&w);

            assert_approx_eq!(r[0][0], -0.7071067782221119905, 1e-14);
            assert_approx_eq!(r[0][1], -0.5656854276809129651, 1e-14);
            assert_approx_eq!(r[0][2], -0.4242640700104211225, 1e-14);

            assert_approx_eq!(r[1][0], 0.5656854276809129651, 1e-14);
            assert_approx_eq!(r[1][1], -0.0925483394532274246, 1e-14);
            assert_approx_eq!(r[1][2], -0.8194112531408833269, 1e-14);

            assert_approx_eq!(r[2][0], 0.4242640700104211225, 1e-14);
            assert_approx_eq!(r[2][1], -0.8194112531408833269, 1e-14);
            assert_approx_eq!(r[2][2], 0.3854415612311154341, 1e-14);
        }

        #[test]
        fn test_rv2m_parity() {
            use rsofa::iauRv2m;
            let mut w = [0.0, 1.41371669, -1.88495559];
            let r = rv2m(&w);

            let mut iau_r = zr();
            unsafe {
                iauRv2m(w.as_mut_ptr(), iau_r.as_mut_ptr());
            }
            assert_eq!(r, iau_r);
        }

        /// t_sofa.c t_rm2v
        #[test]
        fn test_rm2v() {
            let r = [[0.0, -0.8, -0.6], [0.8, -0.36, 0.48], [0.6, 0.48, -0.64]];
            let w = rm2v(&r);

            assert_approx_eq!(w[0], 0.0, 1e-12);
            assert_approx_eq!(w[1], 1.413716694115406957, 1e-12);
            assert_approx_eq!(w[2], -1.884955592153875943, 1e-12);
        }

        #[test]
        fn test_rm2v_parity() {
            use rsofa::iauRm2v;

            let mut r = [[0.0, -0.8, -0.6], [0.8, -0.36, 0.48], [0.6, 0.48, -0.64]];
            let w = rm2v(&r);

            let mut w_iau = [0.0, 0.0, 0.0];
            unsafe {
                iauRm2v(r.as_mut_ptr(), w_iau.as_mut_ptr());
            }
            assert_eq!(w, w_iau);
        }
    }
}
