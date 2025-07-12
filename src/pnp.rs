//!Precession, nutation, polar motion
//! - bias_2000                 BI00      frame bias components, IAU 2000
//! - bias_precession_2000      BP00      frame bias and precession matrices, IAU 2000
//! - bias_precession_2006      BP06      frame bias and precession matrices, IAU 2006
//! - npb_to_xy                 BPN2XY    extract CIP X,Y coordinates from NPB matrix
//! - c2i_2000a                 C2I00A    celestial-to-intermediate matrix, IAU 2000A
//! - c2i_2000b                 C2I00B    celestial-to-intermediate matrix, IAU 2000B
//! - c2i_2006                  C2I06A    celestial-to-intermediate matrix, IAU 2006/2000A
//! - c2i_npb                   C2IBPN    celestial-to-intermediate matrix, given NPB matrix, IAU 2000
//! - c2i_xy                    C2IXY     celestial-to-intermediate matrix, given X,Y, IAU 2000
//! - c2i_xys                   C2IXYS    celestial-to-intermediate matrix, given X,Y and s
//! - c2t_2000a                 C2T00A    celestial-to-terrestrial matrix, IAU 2000A
//! - c2t_2000b                 C2T00B    celestial-to-terrestrial matrix, IAU 2000B
//! - c2t_2006                  C2T06A    celestial-to-terrestrial matrix, IAU 2006/2000A
//! - c2t_cio                   C2TCIO    form CIO-based celestial-to-terrestrial matrix
//! - c2t_eqx                   C2TEQX    form equinox-based celestial-to-terrestrial matrix
//! - c2t_nutation              C2TPE     celestial-to-terrestrial matrix given nutation, IAU 2000
//! - xy_to_c2t                 C2TXY     celestial-to-terrestrial matrix given CIP, IAU 2000
//! - eo_2006                   EO06A     equation of the origins, IAU 2006/2000A
//! - eo_npb_s                  EORS      equation of the origins, given NPB matrix and s
//! - fw_to_rmatrix             FW2M      Fukushima-Williams angles to r-matrix
//! - fw_to_xy                  FW2XY     Fukushima-Williams angles to X,Y
//! - long_term_precession      LTP       long-term precession matrix
//! - long_term_precession_b    LTPB      long-term precession matrix, including ICRS frame bias
//! - long_term_precession_ecl  LTPECL    long-term precession of the ecliptic
//! - long_term_precession_eq   LTPEQU    long-term precession of the equator
//! - nutation_matrix_2000a     NUM00A    nutation matrix, IAU 2000A
//! - nutation_matrix_2000b     NUM00B    nutation matrix, IAU 2000B
//! - nutation_matrix_2006      NUM06A    nutation matrix, IAU 2006/2000A
//! - nutation_matrix           NUMAT     form nutation matrix
//! - nutation_2000a            NUT00A    nutation, IAU 2000A
//! - nutation_2000b            NUT00B    nutation, IAU 2000B
//! - nutation_2006             NUT06A    nutation, IAU 2006/2000A
//! - nutation_1980             NUT80     nutation, IAU 1980
//! - nutation_matrix_1980      NUTM80    nutation matrix, IAU 1980
//! - obliquity_2006            OBL06     mean obliquity, IAU 2006
//! - obliquity_1980            OBL80     mean obliquity, IAU 1980
//! - precession_bias_2006      PB06      zeta,z,theta precession angles, IAU 2006, including bias
//! - precession_fw_2006        PFW06     bias-precession Fukushima-Williams angles, IAU 2006
//! - precession_matrix_2000    PMAT00    precession matrix (including frame bias), IAU 2000
//! - precession_matrix_2006    PMAT06    PB matrix, IAU 2006
//! - precession_matrix_1976    PMAT76    precession matrix, IAU 1976
//! - pn_2000                   PN00      bias/precession/nutation results, IAU 2000
//! - pn_2000a                  PN00A     bias/precession/nutation, IAU 2000A
//! - pn_2000b                  PN00B     bias/precession/nutation, IAU 2000B
//! - pn_2006                   PN06      bias/precession/nutation results, IAU 2006
//! - pn_2006a                  PN06A     bias/precession/nutation results, IAU 2006/2000A
//! - pnm_2000a                 PNM00A    classical NPB matrix, IAU 2000A
//! - pnm_2000b                 PNM00B    classical NPB matrix, IAU 2000B
//! - pnm_2006a                 PNM06A    classical NPB matrix, IAU 2006/2000A
//! - pnm_1980                  PNM80     precession/nutation matrix, IAU 1976/1980
//! - p_2006_equinox            P06E      precession angles, IAU 2006, equinox based
//! - polar_motion_matrix       POM00     polar motion matrix
//! - pr_2000                   PR00      IAU 2000 precession adjustments
//! - prec_1976                 PREC76    accumulated precession angles, IAU 1976
//! - s_2000                    S00       the CIO locator s, given X,Y, IAU 2000A
//! - s_2000a                   S00A      the CIO locator s, IAU 2000A
//! - s_2000b                   S00B      the CIO locator s, IAU 2000B
//! - s_2006                    S06       the CIO locator s, given X,Y, IAU 2006
//! - s_2006a                   S06A      the CIO locator s, IAU 2006/2000A
//! - sp_2000                   SP00      the TIO locator s', IERS 2003
//! - xy_2006a                  XY06      CIP, IAU 2006/2000A, from series
//! - xys_2000a                 XYS00A    CIP and s, IAU 2000A
//! - xys_2000b                 XYS00B    CIP and s, IAU 2000B
//! - xys_2006a                 XYS06A    CIP and s, IAU 2006/2000A
pub fn bias_2000() {
    todo!();
}
pub fn bias_precession_2000() {
    todo!();
}
pub fn bias_precession_2006() {
    todo!();
}
pub fn npb_to_xy() {
    todo!();
}
pub fn c2i_2000a() {
    todo!();
}
pub fn c2i_2000b() {
    todo!();
}
pub fn c2i_2006() {
    todo!();
}
pub fn c2i_npb() {
    todo!();
}
pub fn c2i_xy() {
    todo!();
}
pub fn c2i_xys() {
    todo!();
}
pub fn c2t_2000a() {
    todo!();
}
pub fn c2t_2000b() {
    todo!();
}
pub fn c2t_2006() {
    todo!();
}
pub fn c2t_cio() {
    todo!();
}
pub fn c2t_eqx() {
    todo!();
}
pub fn c2t_nutation() {
    todo!();
}
pub fn xy_to_c2t() {
    todo!();
}
pub fn eo_2006() {
    todo!();
}
pub fn eo_npb_s() {
    todo!();
}
pub fn fw_to_rmatrix() {
    todo!();
}
pub fn fw_to_xy() {
    todo!();
}
pub fn long_term_precession() {
    todo!();
}
pub fn long_term_precession_b() {
    todo!();
}
pub fn long_term_precession_ecl() {
    todo!();
}
pub fn long_term_precession_eq() {
    todo!();
}
pub fn nutation_matrix_2000a() {
    todo!();
}
pub fn nutation_matrix_2000b() {
    todo!();
}
pub fn nutation_matrix_2006() {
    todo!();
}
pub fn nutation_matrix() {
    todo!();
}
pub fn nutation_2000b() {
    todo!();
}
pub fn nutation_2000a() {
    todo!();
}
pub fn nutation_2006() {
    todo!();
}
pub fn nutation_1980() {
    todo!();
}
pub fn nutation_matrix_1980() {
    todo!();
}
pub fn obliquity_2006() {
    todo!();
}
pub fn obliquity_1980() {
    todo!();
}
pub fn precession_bias_2006() {
    todo!();
}
pub fn precession_fw_2006() {
    todo!();
}
pub fn precession_matrix_2000() {
    todo!();
}
pub fn precession_matrix_2006() {
    todo!();
}
pub fn precession_matrix_1976() {
    todo!();
}
pub fn pn_2000() {
    todo!();
}
pub fn pn_2000a() {
    todo!();
}
pub fn pn_2000b() {
    todo!();
}
pub fn pn_2006() {
    todo!();
}
pub fn pn_2006a() {
    todo!();
}
pub fn pnm_2000a() {
    todo!();
}
pub fn pnm_2000b() {
    todo!();
}
pub fn pnm_2006a() {
    todo!();
}
pub fn pnm_1980() {
    todo!();
}
pub fn p_2006_equinox() {
    todo!();
}
pub fn polar_motion_matrix() {
    todo!();
}
pub fn pr_2000() {
    todo!();
}
pub fn prec_1976() {
    todo!();
}
pub fn s_2000() {
    todo!();
}
pub fn s_2000a() {
    todo!();
}
pub fn s_2000b() {
    todo!();
}
pub fn s_2006() {
    todo!();
}
pub fn s_2006a() {
    todo!();
}
pub fn sp_2000() {
    todo!();
}
pub fn xy_2006a() {
    todo!();
}
pub fn xys_2000a() {
    todo!();
}
pub fn xys_2000b() {
    todo!();
}
pub fn xys_2006a() {
    todo!();
}

mod tests {
    use super::*;
    // t_sofa.c t_bi00
    #[test]
    fn test_bias_2000() {
        todo!();
    }
    #[test]
    fn test_bias_2000_parity() {
        use rsofa::iauBi00;
        todo!();
    }

    // t_sofa.c t_bp00
    #[test]
    fn test_bias_precession_2000() {
        todo!();
    }
    #[test]
    fn test_bias_precession_2000_parity() {
        use rsofa::iauBp00;
        todo!();
    }

    // t_sofa.c t_bp06
    #[test]
    fn test_bias_precession_2006() {
        todo!();
    }
    #[test]
    fn test_bias_precession_2006_parity() {
        use rsofa::iauBp06;
        todo!();
    }

    // t_sofa.c t_bpn2xy
    #[test]
    fn test_npb_to_xy() {
        todo!();
    }
    #[test]
    fn test_npb_to_xy_parity() {
        use rsofa::iauBpn2xy;
        todo!();
    }

    // t_sofa.c t_c2i00a
    #[test]
    fn test_c2i_2000a() {
        todo!();
    }
    #[test]
    fn test_c2i_2000a_parity() {
        use rsofa::iauC2i00a;
        todo!();
    }

    // t_sofa.c t_c2i00b
    #[test]
    fn test_c2i_2000b() {
        todo!();
    }
    #[test]
    fn test_c2i_2000b_parity() {
        use rsofa::iauC2i00b;
        todo!();
    }

    // t_sofa.c t_c2i06a
    #[test]
    fn test_c2i_2006() {
        todo!();
    }
    #[test]
    fn test_c2i_2006_parity() {
        use rsofa::iauC2i06a;
        todo!();
    }

    // t_sofa.c t_c2ibpn
    #[test]
    fn test_c2i_npb() {
        todo!();
    }
    #[test]
    fn test_c2i_npb_parity() {
        use rsofa::iauC2ibpn;
        todo!();
    }

    // t_sofa.c t_c2ixy
    #[test]
    fn test_c2i_xy() {
        todo!();
    }
    #[test]
    fn test_c2i_xy_parity() {
        use rsofa::iauC2ixy;
        todo!();
    }

    // t_sofa.c t_c2ixys
    #[test]
    fn test_c2i_xys() {
        todo!();
    }
    #[test]
    fn test_c2i_xys_parity() {
        use rsofa::iauC2ixys;
        todo!();
    }

    // t_sofa.c t_c2t00a
    #[test]
    fn test_c2t_2000a() {
        todo!();
    }
    #[test]
    fn test_c2t_2000a_parity() {
        use rsofa::iauC2t00a;
        todo!();
    }

    // t_sofa.c t_c2t00b
    #[test]
    fn test_c2t_2000b() {
        todo!();
    }
    #[test]
    fn test_c2t_2000b_parity() {
        use rsofa::iauC2t00b;
        todo!();
    }

    // t_sofa.c t_c2t06a
    #[test]
    fn test_c2t_2006() {
        todo!();
    }
    #[test]
    fn test_c2t_2006_parity() {
        use rsofa::iauC2t06a;
        todo!();
    }

    // t_sofa.c t_c2tcio
    #[test]
    fn test_c2t_cio() {
        todo!();
    }
    #[test]
    fn test_c2t_cio_parity() {
        use rsofa::iauC2tcio;
        todo!();
    }

    // t_sofa.c t_c2teqx
    #[test]
    fn test_c2t_eqx() {
        todo!();
    }
    #[test]
    fn test_c2t_eqx_parity() {
        use rsofa::iauC2teqx;
        todo!();
    }

    // t_sofa.c t_c2tpe
    #[test]
    fn test_c2t_nutation() {
        todo!();
    }
    #[test]
    fn test_c2t_nutation_parity() {
        use rsofa::iauC2tpe;
        todo!();
    }

    // t_sofa.c t_c2txy
    #[test]
    fn test_xy_to_c2t() {
        todo!();
    }
    #[test]
    fn test_xy_to_c2t_parity() {
        use rsofa::iauC2txy;
        todo!();
    }

    // t_sofa.c t_eo06a
    #[test]
    fn test_eo_2006() {
        todo!();
    }
    #[test]
    fn test_eo_2006_parity() {
        use rsofa::iauEo06a;
        todo!();
    }

    // t_sofa.c t_eors
    #[test]
    fn test_eo_npb_s() {
        todo!();
    }
    #[test]
    fn test_eo_npb_s_parity() {
        use rsofa::iauEors;
        todo!();
    }

    // t_sofa.c t_fw2m
    #[test]
    fn test_fw_to_rmatrix() {
        todo!();
    }
    #[test]
    fn test_fw_to_rmatrix_parity() {
        use rsofa::iauFw2m;
        todo!();
    }

    // t_sofa.c t_fw2xy
    #[test]
    fn test_fw_to_xy() {
        todo!();
    }
    #[test]
    fn test_fw_to_xy_parity() {
        use rsofa::iauFw2xy;
        todo!();
    }

    // t_sofa.c t_ltp
    #[test]
    fn test_long_term_precession() {
        todo!();
    }
    #[test]
    fn test_long_term_precession_parity() {
        use rsofa::iauLtp;
        todo!();
    }

    // t_sofa.c t_ltpb
    #[test]
    fn test_long_term_precession_b() {
        todo!();
    }
    #[test]
    fn test_long_term_precession_b_parity() {
        use rsofa::iauLtpb;
        todo!();
    }

    // t_sofa.c t_ltpecl
    #[test]
    fn test_long_term_precession_e() {
        todo!();
    }
    #[test]
    fn test_long_term_precession_e_parity() {
        use rsofa::iauLtpecl;
        todo!();
    }

    // t_sofa.c t_ltpequ
    #[test]
    fn test_long_term_precession_eq() {
        todo!();
    }
    #[test]
    fn test_long_term_precession_eq_parity() {
        use rsofa::iauLtpequ;
        todo!();
    }

    // t_sofa.c t_num00a
    #[test]
    fn test_nutation_matrix_2000a() {
        todo!();
    }
    #[test]
    fn test_nutation_matrix_2000a_parity() {
        use rsofa::iauNum00a;
        todo!();
    }

    // t_sofa.c t_num00b
    #[test]
    fn test_nutation_matrix_2000b() {
        todo!();
    }
    #[test]
    fn test_nutation_matrix_2000b_parity() {
        use rsofa::iauNum00b;
        todo!();
    }

    // t_sofa.c t_num06a
    #[test]
    fn test_nutation_matrix_2006() {
        todo!();
    }
    #[test]
    fn test_nutation_matrix_2006_parity() {
        use rsofa::iauNum06a;
        todo!();
    }

    // t_sofa.c t_numat
    #[test]
    fn test_nutation_matrix() {
        todo!();
    }
    #[test]
    fn test_nutation_matrix_parity() {
        use rsofa::iauNumat;
        todo!();
    }

    // t_sofa.c t_nut00a
    #[test]
    fn test_nutation_2000a() {
        todo!();
    }
    #[test]
    fn test_nutation_2000a_parity() {
        use rsofa::iauNut00a;
        todo!();
    }

    // t_sofa.c t_nut00b
    #[test]
    fn test_nutation_2000b() {
        todo!();
    }
    #[test]
    fn test_nutation_2000b_parity() {
        use rsofa::iauNut00b;
        todo!();
    }

    // t_sofa.c t_nut06a
    #[test]
    fn test_nutation_2006() {
        todo!();
    }
    #[test]
    fn test_nutation_2006_parity() {
        use rsofa::iauNut06a;
        todo!();
    }

    // t_sofa.c t_nut80
    #[test]
    fn test_nutation_1980() {
        todo!();
    }
    #[test]
    fn test_nutation_1980_parity() {
        use rsofa::iauNut80;
        todo!();
    }

    // t_sofa.c t_nutm80
    #[test]
    fn test_nutation_matrix_1980() {
        todo!();
    }
    #[test]
    fn test_nutation_matrix_1980_parity() {
        use rsofa::iauNutm80;
        todo!();
    }

    // t_sofa.c t_obl06
    #[test]
    fn test_obliquity_2006() {
        todo!();
    }
    #[test]
    fn test_obliquity_2006_parity() {
        use rsofa::iauObl06;
        todo!();
    }

    // t_sofa.c t_obl80
    #[test]
    fn test_obliquity_1980() {
        todo!();
    }
    #[test]
    fn test_obliquity_1980_parity() {
        use rsofa::iauObl80;
        todo!();
    }

    // t_sofa.c t_pb06
    #[test]
    fn test_precession_bias_2006() {
        todo!();
    }
    #[test]
    fn test_precession_bias_2006_parity() {
        use rsofa::iauPb06;
        todo!();
    }

    // t_sofa.c t_pfw06
    #[test]
    fn test_precession_fw_2006() {
        todo!();
    }
    #[test]
    fn test_precession_fw_2006_parity() {
        use rsofa::iauPfw06;
        todo!();
    }

    // t_sofa.c t_pmat00
    #[test]
    fn test_precession_matrix_2000() {
        todo!();
    }
    #[test]
    fn test_precession_matrix_2000_parity() {
        use rsofa::iauPmat00;
        todo!();
    }

    // t_sofa.c t_pmat06
    #[test]
    fn test_precession_matrix_2006() {
        todo!();
    }
    #[test]
    fn test_precession_matrix_2006_parity() {
        use rsofa::iauPmat06;
        todo!();
    }

    // t_sofa.c t_pmat76
    #[test]
    fn test_precession_matrix_1976() {
        todo!();
    }
    #[test]
    fn test_precession_matrix_1976_parity() {
        use rsofa::iauPmat76;
        todo!();
    }

    // t_sofa.c t_pn00
    #[test]
    fn test_pn_2000() {
        todo!();
    }
    #[test]
    fn test_pn_2000_parity() {
        use rsofa::iauPn00;
        todo!();
    }

    // t_sofa.c t_pn00a
    #[test]
    fn test_pn_2000a() {
        todo!();
    }
    #[test]
    fn test_pn_2000a_parity() {
        use rsofa::iauPn00a;
        todo!();
    }

    // t_sofa.c t_pn00b
    #[test]
    fn test_pn_2000b() {
        todo!();
    }
    #[test]
    fn test_pn_2000b_parity() {
        use rsofa::iauPn00b;
        todo!();
    }

    // t_sofa.c t_pn06
    #[test]
    fn test_pn_2006() {
        todo!();
    }
    #[test]
    fn test_pn_2006_parity() {
        use rsofa::iauPn06;
        todo!();
    }

    // t_sofa.c t_pn06a
    #[test]
    fn test_pn_2006a() {
        todo!();
    }
    #[test]
    fn test_pn_2006a_parity() {
        use rsofa::iauPn06a;
        todo!();
    }

    // t_sofa.c t_pnm00a
    #[test]
    fn test_pnm_2000a() {
        todo!();
    }
    #[test]
    fn test_pnm_2000a_parity() {
        use rsofa::iauPnm00a;
        todo!();
    }

    // t_sofa.c t_pnm00b
    #[test]
    fn test_pnm_2000b() {
        todo!();
    }
    #[test]
    fn test_pnm_2000b_parity() {
        use rsofa::iauPnm00b;
        todo!();
    }

    // t_sofa.c t_pnm06a
    #[test]
    fn test_pnm_2006a() {
        todo!();
    }
    #[test]
    fn test_pnm_2006a_parity() {
        use rsofa::iauPnm06a;
        todo!();
    }

    // t_sofa.c t_pnm80
    #[test]
    fn test_pnm_1980() {
        todo!();
    }
    #[test]
    fn test_pnm_1980_parity() {
        use rsofa::iauPnm80;
        todo!();
    }

    // t_sofa.c t_p06e
    #[test]
    fn test_p_2006_equinox() {
        todo!();
    }
    #[test]
    fn test_p_2006_equinox_parity() {
        use rsofa::iauP06e;
        todo!();
    }

    // t_sofa.c t_pom00
    #[test]
    fn test_polar_motion_matrix() {
        todo!();
    }
    #[test]
    fn test_polar_motion_matrix_parity() {
        use rsofa::iauPom00;
        todo!();
    }

    // t_sofa.c t_pr00
    #[test]
    fn test_pr_2000() {
        todo!();
    }
    #[test]
    fn test_pr_2000_parity() {
        use rsofa::iauPr00;
        todo!();
    }

    // t_sofa.c t_prec76
    #[test]
    fn test_prec_1976() {
        todo!();
    }
    #[test]
    fn test_prec_1976_parity() {
        use rsofa::iauPrec76;
        todo!();
    }

    // t_sofa.c t_s00
    #[test]
    fn test_s_2000() {
        todo!();
    }
    #[test]
    fn test_s_2000_parity() {
        use rsofa::iauS00;
        todo!();
    }

    // t_sofa.c t_s00a
    #[test]
    fn test_s_2000a() {
        todo!();
    }
    #[test]
    fn test_s_2000a_parity() {
        use rsofa::iauS00a;
        todo!();
    }

    // t_sofa.c t_s00b
    #[test]
    fn test_s_2000b() {
        todo!();
    }
    #[test]
    fn test_s_2000b_parity() {
        use rsofa::iauS00b;
        todo!();
    }

    // t_sofa.c t_s06
    #[test]
    fn test_s_2006() {
        todo!();
    }
    #[test]
    fn test_s_2006_parity() {
        use rsofa::iauS06;
        todo!();
    }

    // t_sofa.c t_s06a
    #[test]
    fn test_s_2006a() {
        todo!();
    }
    #[test]
    fn test_s_2006a_parity() {
        use rsofa::iauS06a;
        todo!();
    }

    // t_sofa.c t_sp00
    #[test]
    fn test_sp_2000() {
        todo!();
    }
    #[test]
    fn test_sp_2000_parity() {
        use rsofa::iauSp00;
        todo!();
    }

    // t_sofa.c t_xy06
    #[test]
    fn test_xy_2006a() {
        todo!();
    }
    #[test]
    fn test_xy_2006a_parity() {
        use rsofa::iauXy06;
        todo!();
    }

    // t_sofa.c t_xys00a
    #[test]
    fn test_xys_2000a() {
        todo!();
    }
    #[test]
    fn test_xys_2000a_parity() {
        use rsofa::iauXys00a;
        todo!();
    }

    // t_sofa.c t_xys00b
    #[test]
    fn test_xys_2000b() {
        todo!();
    }
    #[test]
    fn test_xys_2000b_parity() {
        use rsofa::iauXys00b;
        todo!();
    }

    // t_sofa.c t_xys06a
    #[test]
    fn test_xys_2006a() {
        todo!();
    }
    #[test]
    fn test_xys_2006a_parity() {
        use rsofa::iauXys06a;
        todo!();
    }
}
