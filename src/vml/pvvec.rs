//! Operations involving PV-Vectors

/// Initialize
/// - zero_pvvector (ZPV)       zero pv-vector
pub mod initialize {
    use crate::PVvector;

    ///  Zero a pv-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Returned:
    ///     pv       double[2][3]      zero pv-vector
    ///
    ///  Called:
    ///     iauZp        zero p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn zero_pvvector() -> PVvector {
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// t_sofa.c t_zpv
        #[allow(unused_variables)]
        #[test]
        fn test_zpv() {
            let pv = [[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]];
            let pv = zero_pvvector();
            assert_eq!(pv, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        }

        #[allow(unused_variables)]
        #[test]
        fn test_zpv_parity() {
            use rsofa::iauZpv;
            let pv = [[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]];
            let pv = zero_pvvector();

            let mut pv_iau = [[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]];
            unsafe { iauZpv(pv_iau.as_mut_ptr()) }
            assert_eq!(pv, pv_iau);
        }
    }
}

///  Copy/extend/extract
///
///- copy_pvvector (CPV)       copy pv-vector
///- append_zvelocity_pvvector (P2PV)      append zero velocity to p-vector
///- discard_velocity_pvvector (PV2P)      discard velocity component of pv-vector
pub mod copy_extend_extract {
    use crate::PVvector;
    use crate::Pvector;
    use crate::vml::pvrm::copy::cp;
    use crate::vml::pvrm::initialize::zero_pvector;

    ///  Copy a position/velocity vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     pv     double[2][3]    position/velocity vector to be copied
    ///
    ///  Returned:
    ///     c      double[2][3]    copy
    ///
    ///  Called:
    ///     iauCp        copy p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn copy_pvvector(pv: &PVvector) -> PVvector {
        // TODO: pointless as vectors implement the copy trait
        [cp(&pv[0]), cp(&pv[1])]
    }

    ///  Extend a p-vector to a pv-vector by appending a zero velocity.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     p        double[3]       p-vector
    ///
    ///  Returned:
    ///     pv       double[2][3]    pv-vector
    ///
    ///  Called:
    ///     iauCp        copy p-vector
    ///     iauZp        zero p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn append_zvelocity_pvvector(p: &Pvector) -> PVvector {
        [cp(p), zero_pvector()]
    }

    ///  Discard velocity component of a pv-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     pv      double[2][3]     pv-vector
    ///
    ///  Returned:
    ///     p       double[3]        p-vector
    ///
    ///  Called:
    ///     iauCp        copy p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn discard_velocity_pvvector(pv: &PVvector) -> Pvector {
        cp(&pv[0])
    }

    #[cfg(test)]
    pub mod tests {
        use super::*;
        use crate::vml::pvvec::initialize::zero_pvvector;

        /// t_sofa.c t_cpv
        #[test]
        fn test_cpv() {
            let pv = [[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]];
            let copied_pv = copy_pvvector(&pv);
            assert_eq!(pv, copied_pv);
        }

        #[test]
        fn test_cpv_parity() {
            use rsofa::iauCpv;
            let mut pv = [[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]];
            let copied_pv = copy_pvvector(&pv);

            let mut copied_pv_iau = zero_pvvector();
            unsafe {
                iauCpv(pv.as_mut_ptr(), copied_pv_iau.as_mut_ptr());
            }
            assert_eq!(copied_pv, copied_pv_iau);
        }

        /// t_sofac. t_p2pv
        #[test]
        fn test_p2pv() {
            let p = [0.25, 1.2, 3.0];
            let pv = append_zvelocity_pvvector(&p);
            assert_eq!(pv[0], p);
            assert_eq!(pv[1], [0.0; 3]);
        }

        #[test]
        fn test_p2pv_parity() {
            use rsofa::iauP2pv;
            let mut p = [0.25, 1.2, 3.0];
            let pv = append_zvelocity_pvvector(&p);

            let mut pv_iau = zero_pvvector();
            unsafe {
                iauP2pv(p.as_mut_ptr(), pv_iau.as_mut_ptr());
            }
            assert_eq!(pv, pv_iau);
        }

        /// t_sofa.c t_pv2p
        #[test]
        fn test_pv2p() {
            let pv = [[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]];
            let p = discard_velocity_pvvector(&pv);
            assert_eq!(p, pv[0]);
        }

        #[test]
        fn test_pv2p_parity() {
            use rsofa::iauPv2p;
            let mut pv = [[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]];
            let p = discard_velocity_pvvector(&pv);

            let mut p_iau = zero_pvector();
            unsafe {
                iauPv2p(pv.as_mut_ptr(), p_iau.as_mut_ptr());
            }
            assert_eq!(p, p_iau);
        }
    }
}
///Spherical/Cartesian conversions
///
/// - spherical_to_pvvector (S2PV)      spherical to pv-vector
/// - pvvector_to_spherical (PV2S)      pv-vector to spherical
pub mod sphere_cart_conv {
    use crate::PVvector;

    ///  Convert position/velocity from spherical to Cartesian coordinates.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     theta    double          longitude angle (radians)
    ///     phi      double          latitude angle (radians)
    ///     r        double          radial distance
    ///     td       double          rate of change of theta
    ///     pd       double          rate of change of phi
    ///     rd       double          rate of change of r
    ///
    ///  Returned:
    ///     pv       double[2][3]    pv-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn spherical_to_pvvector(
        theta: f64,
        phi: f64,
        r: f64,
        td: f64,
        pd: f64,
        rd: f64,
    ) -> PVvector {
        let st = theta.sin();
        let ct = theta.cos();
        let sp = phi.sin();
        let cp = phi.cos();
        let rcp = r * cp;
        let x = rcp * ct;
        let y = rcp * st;
        let rpd = r * pd;
        let w = rpd * sp - cp * rd;

        [
            [x, y, r * sp],
            [-y * td - w * ct, x * td - w * st, rpd * cp + sp * rd],
        ]
    }

    ///  Convert position/velocity from Cartesian to spherical coordinates.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     pv       double[2][3]  pv-vector
    ///
    ///  Returned:
    ///     theta    double        longitude angle (radians)
    ///     phi      double        latitude angle (radians)
    ///     r        double        radial distance
    ///     td       double        rate of change of theta
    ///     pd       double        rate of change of phi
    ///     rd       double        rate of change of r
    ///
    ///  Notes:
    ///
    ///  1) If the position part of pv is null, theta, phi, td and pd
    ///     are indeterminate.  This is handled by extrapolating the
    ///     position through unit time by using the velocity part of
    ///     pv.  This moves the origin without changing the direction
    ///     of the velocity component.  If the position and velocity
    ///     components of pv are both null, zeroes are returned for all
    ///     six results.
    ///
    ///  2) If the position is a pole, theta, td and pd are indeterminate.
    ///     In such cases zeroes are returned for all three.
    ///
    ///  This revision:  2021 May 11
    pub fn pvvector_to_spherical(pv: &PVvector) -> (f64, f64, f64, f64, f64, f64) {
        // Return values
        let theta: f64;
        let phi: f64;
        #[allow(clippy::needless_late_init)]
        let r: f64;
        let td: f64;
        let pd: f64;
        #[allow(clippy::needless_late_init)]
        let rd: f64;

        // Components of position/velocity vector.
        let mut x = pv[0][0];
        let mut y = pv[0][1];
        let mut z = pv[0][2];
        let xd = pv[1][0];
        let yd = pv[1][1];
        let zd = pv[1][2];

        // Component of r in XY plane squared.
        let mut rxy2 = x * x + y * y;

        // Modulus squared.
        let mut r2 = rxy2 + z * z;

        // Modulus.
        let rtrue = r2.sqrt();

        // If null vector, move the origin along the direction of movement.
        let mut rw = rtrue;
        if rtrue == 0.0 {
            x = xd;
            y = yd;
            z = zd;
            rxy2 = x * x + y * y;
            r2 = rxy2 + z * z;
            rw = r2.sqrt();
        }

        // Position and velocity in spherical coordinates.
        let rxy = rxy2.sqrt();
        let xyp = x * xd + y * yd;
        if rxy2 != 0.0 {
            theta = y.atan2(x);
            phi = z.atan2(rxy);
            td = (x * yd - y * xd) / rxy2;
            pd = (zd * rxy2 - z * xyp) / (r2 * rxy);
        } else {
            theta = 0.0;
            phi = if z != 0.0 { z.atan2(rxy) } else { 0.0 };
            td = 0.0;
            pd = 0.0;
        }
        r = rtrue;
        rd = if rw != 0.0 { (xyp + z * zd) / rw } else { 0.0 };

        (theta, phi, r, td, pd, rd)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::vml::pvvec::initialize::zero_pvvector;
        use assert_approx_eq::assert_approx_eq;

        /// t_sofa.c t_s2pv
        #[test]
        fn test_s2pv() {
            let pv = spherical_to_pvvector(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5);
            assert_approx_eq!(pv[0][0], -0.4514964673880165228, 1e-12);
            assert_approx_eq!(pv[0][1], 0.0309339427734258688, 1e-12);
            assert_approx_eq!(pv[0][2], 0.0559466810510877933, 1e-12);

            assert_approx_eq!(pv[1][0], 0.1292270850663260170e-4, 1e-16);
            assert_approx_eq!(pv[1][1], 0.2652814182060691422e-5, 1e-16);
            assert_approx_eq!(pv[1][2], 0.2568431853930292259e-5, 1e-16);
        }

        #[test]
        fn test_s2pv_parity() {
            use rsofa::iauS2pv;
            let pv = spherical_to_pvvector(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5);
            let mut pv_iau = zero_pvvector();
            unsafe {
                iauS2pv(
                    -3.21,
                    0.123,
                    0.456,
                    -7.8e-6,
                    9.01e-6,
                    -1.23e-5,
                    pv_iau.as_mut_ptr(),
                );
            }
            assert_eq!(pv, pv_iau);
        }

        /// t_sofa.c t_pv2s
        #[test]
        fn test_pv2s() {
            let pv = [
                [
                    -0.4514964673880165,
                    0.03093394277342585,
                    0.05594668105108779,
                ],
                [
                    1.292270850663260e-5,
                    2.652814182060692e-6,
                    2.568431853930293e-6,
                ],
            ];
            let (theta, phi, r, td, pd, rd) = pvvector_to_spherical(&pv);
            assert_approx_eq!(theta, 3.073185307179586515, 1e-12);
            assert_approx_eq!(phi, 0.1229999999999999992, 1e-12);
            assert_approx_eq!(r, 0.4559999999999999757, 1e-12);
            assert_approx_eq!(td, -0.7800000000000000364e-5, 1e-16);
            assert_approx_eq!(pd, 0.9010000000000001639e-5, 1e-16);
            assert_approx_eq!(rd, -0.1229999999999999832e-4, 1e-16);
        }

        #[test]
        fn test_pv2s_parity() {
            use rsofa::iauPv2s;
            let mut pv = [
                [
                    -0.4514964673880165,
                    0.03093394277342585,
                    0.05594668105108779,
                ],
                [
                    1.292270850663260e-5,
                    2.652814182060692e-6,
                    2.568431853930293e-6,
                ],
            ];
            let (theta, phi, r, td, pd, rd) = pvvector_to_spherical(&pv);

            let mut theta_iau = 0.0;
            let mut phi_iau = 0.0;
            let mut r_iau = 0.0;
            let mut td_iau = 0.0;
            let mut pd_iau = 0.0;
            let mut rd_iau = 0.0;
            unsafe {
                iauPv2s(
                    pv.as_mut_ptr(),
                    &mut theta_iau,
                    &mut phi_iau,
                    &mut r_iau,
                    &mut td_iau,
                    &mut pd_iau,
                    &mut rd_iau,
                );
            }
            assert_eq!(theta, theta_iau);
            assert_eq!(phi, phi_iau);
            assert_eq!(r, r_iau);
            assert_eq!(td, td_iau);
            assert_eq!(pd, pd_iau);
            assert_eq!(rd, rd_iau);
        }
    }
}
///Operations on pv-vectors
///
/// - pvvector_plus_pvvector (PVPPV)     pv-vector plus pv-vector
/// - pvvector_minus_pvvector (PVMPV)     pv-vector minus pv-vector
/// - pvvector_dot_pvvector (PVDPV)     inner (=scalar=dot) product of two pv-vectors
/// - pvvector_cross_pvvector (PVXPV)     outer (=vector=cross) product of two pv-vectors
/// - pvvector_modulus (PVM)       modulus of pv-vector
/// - pvvector_multiply_scalar (SXPV)      multiply pv-vector by scalar
/// - pvvector_multiply_two_scalar (S2XPV)     multiply pv-vector by two scalars
/// - pvvector_update (PVU)       update pv-vector
/// - pvvector_update_discard_velocity (PVUP)      update pv-vector discarding velocity
pub mod pvvector_ops {
    use crate::vml::pvrm::vec_ops::{
        pvector_cross_product, pvector_dot_product, pvector_minus_pvector, pvector_modulus,
        pvector_multiply_scalar, pvector_plus_pvector, pvector_plus_scaled_pvector,
    };
    use crate::{PVvector, Pvector, cpv};

    ///  Add one pv-vector to another.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double[2][3]      first pv-vector
    ///     b        double[2][3]      second pv-vector
    ///
    ///  Returned:
    ///     apb      double[2][3]      a + b
    ///
    ///  Note:
    ///     It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  Called:
    ///     iauPpp       p-vector plus p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_plus_pvvector(a: &PVvector, b: &PVvector) -> PVvector {
        [
            pvector_plus_pvector(&a[0], &b[0]),
            pvector_plus_pvector(&a[1], &b[1]),
        ]
    }

    ///  Subtract one pv-vector from another.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a       double[2][3]      first pv-vector
    ///     b       double[2][3]      second pv-vector
    ///
    ///  Returned:
    ///     amb     double[2][3]      a - b
    ///
    ///  Note:
    ///     It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  Called:
    ///     iauPmp       p-vector minus p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_minus_pvvector(a: &PVvector, b: &PVvector) -> PVvector {
        [
            pvector_minus_pvector(&a[0], &b[0]),
            pvector_minus_pvector(&a[1], &b[1]),
        ]
    }

    ///  Inner (=scalar=dot) product of two pv-vectors.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double[2][3]      first pv-vector
    ///     b        double[2][3]      second pv-vector
    ///
    ///  Returned:
    ///     adb      double[2]         a . b (see note)
    ///
    ///  Note:
    ///
    ///     If the position and velocity components of the two pv-vectors are
    ///     ( ap, av ) and ( bp, bv ), the result, a . b, is the pair of
    ///     numbers ( ap . bp , ap . bv + av . bp ).  The two numbers are the
    ///     dot-product of the two p-vectors and its derivative.
    ///
    ///  Called:
    ///     iauPdp       scalar product of two p-vectors
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_dot_pvvector(a: &PVvector, b: &PVvector) -> [f64; 2] {
        // a . b = constant part of result.
        let adb_const = pvector_dot_product(&a[0], &b[0]);

        // a . bdot.
        let adbd = pvector_dot_product(&a[0], &b[1]);

        // adot . b
        let addb = pvector_dot_product(&a[1], &b[0]);

        // Velocity part of result.
        let adb_vel = adbd + addb;

        [adb_const, adb_vel]
    }

    ///  Outer (=vector=cross) product of two pv-vectors.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double[2][3]      first pv-vector
    ///     b        double[2][3]      second pv-vector
    ///
    ///  Returned:
    ///     axb      double[2][3]      a x b
    ///
    ///  Notes:
    ///
    ///  1) If the position and velocity components of the two pv-vectors are
    ///     ( ap, av ) and ( bp, bv ), the result, a x b, is the pair of
    ///     vectors ( ap x bp, ap x bv + av x bp ).  The two vectors are the
    ///     cross-product of the two p-vectors and its derivative.
    ///
    ///  2) It is permissible to re-use the same array for any of the
    ///     arguments.
    ///
    ///  Called:
    ///     iauCpv       copy pv-vector
    ///     iauPxp       vector product of two p-vectors
    ///     iauPpp       p-vector plus p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_cross_pvvector(a: &PVvector, b: &PVvector) -> PVvector {
        // Make copies of the inputs.
        let wa = cpv(a);
        let wb = cpv(b);

        // a x b = position part of result.
        let axb_pos = pvector_cross_product(&wa[0], &wb[0]);

        // a x bdot + adot x b = velocity part of result.
        let axbd = pvector_cross_product(&wa[0], &wb[1]);
        let adxb = pvector_cross_product(&wa[1], &wb[0]);
        let axb_vel = pvector_plus_pvector(&axbd, &adxb);

        [axb_pos, axb_vel]
    }

    ///  Modulus of pv-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     pv     double[2][3]   pv-vector
    ///
    ///  Returned:
    ///     r      double         modulus of position component
    ///     s      double         modulus of velocity component
    ///
    ///  Called:
    ///     iauPm        modulus of p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_modulus(pv: &PVvector) -> (f64, f64) {
        // Distance
        let r = pvector_modulus(&pv[0]);

        // Speed
        let s = pvector_modulus(&pv[1]);

        (r, s)
    }
    ///  Multiply a pv-vector by a scalar.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     s       double          scalar
    ///     pv      double[2][3]    pv-vector
    ///
    ///  Returned:
    ///     spv     double[2][3]    s * pv
    ///
    ///  Note:
    ///     It is permissible for pv and spv to be the same array.
    ///
    ///  Called:
    ///     iauS2xpv     multiply pv-vector by two scalars
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_multiply_scalar(s: f64, pv: &PVvector) -> PVvector {
        pvvector_multiply_two_scalar(s, s, pv)
    }

    ///  Multiply a pv-vector by two scalars.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     s1     double         scalar to multiply position component by
    ///     s2     double         scalar to multiply velocity component by
    ///     pv     double[2][3]   pv-vector
    ///
    ///  Returned:
    ///     spv    double[2][3]   pv-vector: p scaled by s1, v scaled by s2
    ///
    ///  Note:
    ///     It is permissible for pv and spv to be the same array.
    ///
    ///  Called:
    ///     iauSxp       multiply p-vector by scalar
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_multiply_two_scalar(s1: f64, s2: f64, pv: &PVvector) -> PVvector {
        [
            pvector_multiply_scalar(s1, &pv[0]),
            pvector_multiply_scalar(s2, &pv[1]),
        ]
    }
    ///  Update a pv-vector.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     dt       double           time interval
    ///     pv       double[2][3]     pv-vector
    ///
    ///  Returned:
    ///     upv      double[2][3]     p updated, v unchanged
    ///
    ///  Notes:
    ///
    ///  1) "Update" means "refer the position component of the vector
    ///     to a new date dt time units from the existing date".
    ///
    ///  2) The time units of dt must match those of the velocity.
    ///
    ///  3) It is permissible for pv and upv to be the same array.
    ///
    ///  Called:
    ///     iauPpsp      p-vector plus scaled p-vector
    ///     iauCp        copy p-vector
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_update(dt: f64, pv: PVvector) -> PVvector {
        [pvector_plus_scaled_pvector(&pv[0], dt, &pv[1]), pv[1]]
    }

    ///  Update a pv-vector, discarding the velocity component.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     dt       double            time interval
    ///     pv       double[2][3]      pv-vector
    ///
    ///  Returned:
    ///     p        double[3]         p-vector
    ///
    ///  Notes:
    ///
    ///  1) "Update" means "refer the position component of the vector to a
    ///     new date dt time units from the existing date".
    ///
    ///  2) The time units of dt must match those of the velocity.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn pvvector_update_discard_velocity(dt: f64, pv: &PVvector) -> Pvector {
        [
            pv[0][0] + dt * pv[1][0],
            pv[0][1] + dt * pv[1][1],
            pv[0][2] + dt * pv[1][2],
        ]
    }

    #[cfg(test)]
    mod tests {
        use crate::vml::pvrm::initialize::zero_pvector;
        use crate::vml::pvvec::initialize::zero_pvvector;
        use assert_approx_eq::assert_approx_eq;

        use super::*;

        // t_sofa.c t_pvxpv
        #[test]
        fn test_pvxpv() {
            let a = [[2.0, 2.0, 3.0], [6.0, 0.0, 4.0]];
            let b = [[1.0, 3.0, 4.0], [0.0, 2.0, 8.0]];
            let axb = pvvector_cross_pvvector(&a, &b);
            assert_eq!(axb[0], [-1.0, -5.0, 4.0]);
            assert_eq!(axb[1], [-2.0, -36.0, 22.0]);
        }

        #[test]
        fn test_pvxpv_parity() {
            use rsofa::iauPvxpv;
            let mut a = [[2.0, 2.0, 3.0], [6.0, 0.0, 4.0]];
            let mut b = [[1.0, 3.0, 4.0], [0.0, 2.0, 8.0]];
            let axb = pvvector_cross_pvvector(&a, &b);

            let mut axb_iau = zero_pvvector();
            unsafe {
                iauPvxpv(a.as_mut_ptr(), b.as_mut_ptr(), axb_iau.as_mut_ptr());
            }
            assert_eq!(axb, axb_iau);
        }

        // t_sofa.c t_pvdpv
        #[test]
        fn test_pvdpv() {
            let a = [[2.0, 2.0, 3.0], [6.0, 0.0, 4.0]];
            let b = [[1.0, 3.0, 4.0], [0.0, 2.0, 8.0]];
            let adb = pvvector_dot_pvvector(&a, &b);
            assert_eq!(adb[0], 20.0);
            assert_eq!(adb[1], 50.0);
        }

        #[test]
        fn test_pvdpv_parity() {
            use rsofa::iauPvdpv;
            let mut a = [[2.0, 2.0, 3.0], [6.0, 0.0, 4.0]];
            let mut b = [[1.0, 3.0, 4.0], [0.0, 2.0, 8.0]];
            let adb = pvvector_dot_pvvector(&a, &b);

            let mut adb_iau = [0.0; 2];
            unsafe {
                iauPvdpv(a.as_mut_ptr(), b.as_mut_ptr(), adb_iau.as_mut_ptr());
            }
            assert_eq!(adb, adb_iau);
        }

        // t_sofa.c t_pvmpv
        #[test]
        fn test_pvmpv() {
            let a = [[2.0, 2.0, 3.0], [5.0, 6.0, 3.0]];
            let b = [[1.0, 3.0, 4.0], [3.0, 2.0, 1.0]];
            let amb = pvvector_minus_pvvector(&a, &b);
            assert_eq!(amb[0][0], 1.0);
            assert_eq!(amb[0][1], -1.0);
            assert_eq!(amb[0][2], -1.0);

            assert_eq!(amb[1][0], 2.0);
            assert_eq!(amb[1][1], 4.0);
            assert_eq!(amb[1][2], 2.0);
        }

        #[test]
        fn test_pvmpv_parity() {
            use rsofa::iauPvmpv;
            let mut a = [[2.0, 2.0, 3.0], [5.0, 6.0, 3.0]];
            let mut b = [[1.0, 3.0, 4.0], [3.0, 2.0, 1.0]];
            let amb = pvvector_minus_pvvector(&a, &b);

            let mut amb_iau = zero_pvvector();
            unsafe {
                iauPvmpv(a.as_mut_ptr(), b.as_mut_ptr(), amb_iau.as_mut_ptr());
            }
            assert_eq!(amb, amb_iau);
        }

        // t_sofa.c t_pvm
        #[test]
        fn test_pvm() {
            let pv = [[0.3, 1.2, -2.5], [0.45, -0.25, 1.1]];
            let (r, s) = pvvector_modulus(&pv);
            assert_approx_eq!(r, 2.789265136196270604, 1e-12);
            assert_approx_eq!(s, 1.214495780149111922, 1e-12);
        }

        #[test]
        fn test_pvm_parity() {
            use rsofa::iauPvm;
            let mut pv = [[0.3, 1.2, -2.5], [0.45, -0.25, 1.1]];
            let (r, s) = pvvector_modulus(&pv);
            let mut r_iau = 0.0;
            let mut s_iau = 0.0;
            unsafe {
                iauPvm(pv.as_mut_ptr(), &mut r_iau, &mut s_iau);
            }
            assert_eq!(r, r_iau);
            assert_eq!(s, s_iau);
        }

        // t_sofa.c t_sxpv
        #[test]
        fn test_sxpv() {
            let s = 2.0;
            let pv = [[0.3, 1.2, -2.5], [0.5, 3.2, -0.7]];

            let spv = pvvector_multiply_scalar(s, &pv);
            assert_eq!(spv, [[0.6, 2.4, -5.0], [1.0, 6.4, -1.4]]);
        }

        #[test]
        fn test_sxpv_parity() {
            use rsofa::iauSxpv;
            let s = 2.0;
            let mut pv = [[0.3, 1.2, -2.5], [0.5, 3.2, -0.7]];

            let spv = pvvector_multiply_scalar(s, &pv);
            let mut spv_iau = zero_pvvector();
            unsafe {
                iauSxpv(s, pv.as_mut_ptr(), spv_iau.as_mut_ptr());
            }
            assert_eq!(spv, spv_iau);
        }
        // t_sofa.c t_s2xpv
        #[test]
        fn test_s2xpv() {
            let s1 = 2.0;
            let s2 = 3.0;
            let pv = [[0.3, 1.2, -2.5], [0.5, 2.3, -0.4]];
            let spv = pvvector_multiply_two_scalar(s1, s2, &pv);
            assert_approx_eq!(spv[0][0], 0.6, 1e-12);
            assert_approx_eq!(spv[0][1], 2.4, 1e-12);
            assert_approx_eq!(spv[0][2], -5.0, 1e-12);

            assert_approx_eq!(spv[1][0], 1.5, 1e-12);
            assert_approx_eq!(spv[1][1], 6.9, 1e-12);
            assert_approx_eq!(spv[1][2], -1.2, 1e-12);
        }

        #[test]
        fn test_s2xpv_parity() {
            use rsofa::iauS2xpv;
            let s1 = 2.0;
            let s2 = 3.0;
            let mut pv = [[0.3, 1.2, -2.5], [0.5, 2.3, -0.4]];
            let spv = pvvector_multiply_two_scalar(s1, s2, &pv);

            let mut spv_iau = zero_pvvector();
            unsafe {
                iauS2xpv(s1, s2, pv.as_mut_ptr(), spv_iau.as_mut_ptr());
            }
            assert_eq!(spv, spv_iau);
        }

        // t_sofa.c t_pvppv
        #[test]
        fn test_pvppv() {
            let a = [[2.0, 2.0, 3.0], [5.0, 6.0, 3.0]];
            let b = [[1.0, 3.0, 4.0], [3.0, 2.0, 1.0]];
            let apb = pvvector_plus_pvvector(&a, &b);
            assert_eq!(apb[0][0], 3.0);
            assert_eq!(apb[0][1], 5.0);
            assert_eq!(apb[0][2], 7.0);

            assert_eq!(apb[1][0], 8.0);
            assert_eq!(apb[1][1], 8.0);
            assert_eq!(apb[1][2], 4.0);
        }
        #[test]
        fn test_pvppv_parity() {
            use rsofa::iauPvppv;
            let mut a = [[2.0, 2.0, 3.0], [5.0, 6.0, 3.0]];
            let mut b = [[1.0, 3.0, 4.0], [3.0, 2.0, 1.0]];
            let apb = pvvector_plus_pvvector(&a, &b);

            let mut apb_iau = zero_pvvector();
            unsafe {
                iauPvppv(a.as_mut_ptr(), b.as_mut_ptr(), apb_iau.as_mut_ptr());
            }
            assert_eq!(apb, apb_iau);
        }

        // t_sofa.c t_pvu
        #[test]
        fn test_pvu() {
            let pv = [
                [
                    126668.5912743160734,
                    2136.792716839935565,
                    -245251.2339876830229,
                ],
                [
                    -0.4051854035740713039e-2,
                    -0.6253919754866175788e-2,
                    0.1189353719774107615e-1,
                ],
            ];
            let upv = pvvector_update(2920.0, pv);
            assert_approx_eq!(upv[0][0], 126656.7598605317105, 1e-6);
            assert_approx_eq!(upv[0][1], 2118.531271155726332, 1e-8);
            assert_approx_eq!(upv[0][2], -245216.5048590656190, 1e-6);

            assert_approx_eq!(upv[1][0], -0.4051854035740713039e-2, 1e-12);
            assert_approx_eq!(upv[1][1], -0.6253919754866175788e-2, 1e-12);
            assert_approx_eq!(upv[1][2], 0.1189353719774107615e-1, 1e-12);
        }

        #[test]
        fn test_pvu_parity() {
            use rsofa::iauPvu;
            let mut pv = [
                [
                    126668.5912743160734,
                    2136.792716839935565,
                    -245251.2339876830229,
                ],
                [
                    -0.4051854035740713039e-2,
                    -0.6253919754866175788e-2,
                    0.1189353719774107615e-1,
                ],
            ];
            let upv = pvvector_update(2920.0, pv);

            let mut upv_iau = zero_pvvector();
            unsafe {
                iauPvu(2920.0, pv.as_mut_ptr(), upv_iau.as_mut_ptr());
            }
            assert_eq!(upv, upv_iau);
        }
        // t_sofa.c t_pvup
        #[test]
        fn test_pvup() {
            let pv = [
                [
                    126668.5912743160734,
                    2136.792716839935565,
                    -245251.2339876830229,
                ],
                [
                    -0.4051854035740713039e-2,
                    -0.6253919754866175788e-2,
                    0.1189353719774107615e-1,
                ],
            ];
            let p = pvvector_update_discard_velocity(2920.0, &pv);
            assert_approx_eq!(p[0], 126656.7598605317105, 1e-6);
            assert_approx_eq!(p[1], 2118.531271155726332, 1e-8);
            assert_approx_eq!(p[2], -245216.5048590656190, 1e-6);
        }

        #[test]
        fn test_pvup_parity() {
            use rsofa::iauPvup;
            let mut pv = [
                [
                    126668.5912743160734,
                    2136.792716839935565,
                    -245251.2339876830229,
                ],
                [
                    -0.4051854035740713039e-2,
                    -0.6253919754866175788e-2,
                    0.1189353719774107615e-1,
                ],
            ];
            let p = pvvector_update_discard_velocity(2920.0, &pv);
            let mut p_iau = zero_pvector();
            unsafe {
                iauPvup(2920.0, pv.as_mut_ptr(), p_iau.as_mut_ptr());
            }
            assert_eq!(p, p_iau);
        }
    }
}
