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
