//! Ecliptic coordinates
//! - ec_to_ircs_2006         ECEQ06    ecliptic to ICRS, IAU 2006
//! - rm_ec_to_ircs_2006      ECM06     rotation matrix, ICRS to ecliptic, IAU 2006
//! - icrs_to_ec_2006         EQEC06    ICRS to ecliptic, IAU 2006
//! - ec_to_icrs_long         LTECEQ    ecliptic to ICRS, long term
//! - rm_ec_to_icrs_long      LTECM     rotation matrix, ICRS to ecliptic, long-term
//! - icrs_to_ec_long         LTEQEC    ICRS to ecliptic, long term
pub fn ec_to_ircs_2006() {
    todo!();
}
pub fn rm_ec_to_ircs_2006() {
    todo!();
}
pub fn icrs_to_ec_2006() {
    todo!();
}
pub fn ec_to_icrs_long() {
    todo!();
}
pub fn rm_ec_to_icrs_long() {
    todo!();
}
pub fn icrs_to_ec_long() {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;

    // t_sofa.c t_eceq06
    #[test]
    fn test_ec_to_ircs_2006() {
        todo!();
    }
    #[test]
    fn test_ec_to_ircs_2006_parity() {
        use rsofa::iauEceq06;
        todo!();
    }

    // t_sofa.c t_ecm06
    #[test]
    fn test_rm_ec_to_ircs_2006() {
        todo!();
    }
    #[test]
    fn test_rm_ec_to_ircs_2006_parity() {
        use rsofa::iauEcm06;
        todo!();
    }

    // t_sofa.c t_eqec06
    #[test]
    fn test_icrs_to_ec_2006() {
        todo!();
    }
    #[test]
    fn test_icrs_to_ec_2006_parity() {
        use rsofa::iauEqec06;
        todo!();
    }

    // t_sofa.c t_lteceq
    #[test]
    fn test_ec_to_icrs_long() {
        todo!();
    }
    #[test]
    fn test_ec_to_icrs_long_parity() {
        use rsofa::iauLteceq;
        todo!();
    }

    // t_sofa.c t_ltecm
    #[test]
    fn test_rm_ec_to_icrs_long() {
        todo!();
    }
    #[test]
    fn test_rm_ec_to_icrs_long_parity() {
        use rsofa::iauLtecm;
        todo!();
    }

    // t_sofa.c t_lteqec
    #[test]
    fn test_icrs_to_ec_long() {
        todo!();
    }
    #[test]
    fn test_icrs_to_ec_long_parity() {
        use rsofa::iauLteqec;
        todo!();
    }
}
