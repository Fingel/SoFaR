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
