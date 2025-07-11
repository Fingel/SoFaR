//!Earth rotation angle and sidereal time
//! - ee_2000     EE00      equation of the equinoxes, IAU 2000
//! - ee_2000a    EE00A     equation of the equinoxes, IAU 2000A
//! - ee_2000b    EE00B     equation of the equinoxes, IAU 2000B
//! - ee_2006a    EE06A     equation of the equinoxes, IAU 2006/2000A
//! - ee_ct_2000  EECT00    equation of the equinoxes complementary terms, IAU 2000
//! - ee_1994     EQEQ94    equation of the equinoxes, IAU 1994
//! - era_2000    ERA00     Earth rotation angle, IAU 2000
//! - gmst_2000   GMST00    Greenwich mean sidereal time, IAU 2000
//! - gmst_2006   GMST06    Greenwich mean sidereal time, IAU 2006
//! - gmst_1982   GMST82    Greenwich mean sidereal time, IAU 1982
//! - gst_2000a   GST00A    Greenwich apparent sidereal time, IAU 2000A
//! - gst_2000b   GST00B    Greenwich apparent sidereal time, IAU 2000B
//! - gst_2006    GST06     Greenwich apparent ST, IAU 2006, given NPB matrix
//! - gst_2006a   GST06A    Greenwich apparent sidereal time, IAU 2006/2000A
//! - gst_1994    GST94     Greenwich apparent sidereal time, IAU 1994
pub fn ee_2000() {
    todo!();
}
pub fn ee_2000a() {
    todo!();
}
pub fn ee_2000b() {
    todo!();
}
pub fn ee_2006a() {
    todo!();
}
pub fn ee_ct_2000() {
    todo!();
}
pub fn ee_1994() {
    todo!();
}
pub fn era_2000() {
    todo!();
}
pub fn gmst_2000() {
    todo!();
}
pub fn gmst_2006() {
    todo!();
}
pub fn gmst_1982() {
    todo!();
}
pub fn gst_2000a() {
    todo!();
}
pub fn gst_2000b() {
    todo!();
}
pub fn gst_2006() {
    todo!();
}
pub fn gst_2006a() {
    todo!();
}
pub fn gst_1994() {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;

    // t_sofa.c t_ee00
    #[test]
    fn test_ee_2000() {
        todo!();
    }
    #[test]
    fn test_ee_2000_parity() {
        use rsofa::iauEe00;
        todo!();
    }

    // t_sofa.c t_ee00a
    #[test]
    fn test_ee_2000a() {
        todo!();
    }
    #[test]
    fn test_ee_2000a_parity() {
        use rsofa::iauEe00a;
        todo!();
    }

    // t_sofa.c t_ee00b
    #[test]
    fn test_ee_2000b() {
        todo!();
    }
    #[test]
    fn test_ee_2000b_parity() {
        use rsofa::iauEe00b;
        todo!();
    }

    // t_sofa.c t_ee06a
    #[test]
    fn test_ee_2006a() {
        todo!();
    }
    #[test]
    fn test_ee_2006a_parity() {
        use rsofa::iauEe06a;
        todo!();
    }

    // t_sofa.c t_eect00
    #[test]
    fn test_ee_ct_2000() {
        todo!();
    }
    #[test]
    fn test_ee_ct_200_parity() {
        use rsofa::iauEect00;
        todo!();
    }

    // t_sofa.c t_eqeq94
    #[test]
    fn test_ee_1994() {
        todo!();
    }
    #[test]
    fn test_ee_1994_parity() {
        use rsofa::iauEqeq94;
        todo!();
    }

    // t_sofa.c t_era00
    #[test]
    fn test_era_2000() {
        todo!();
    }
    #[test]
    fn test_era_2000_parity() {
        use rsofa::iauEra00;
        todo!();
    }

    // t_sofa.c t_gmst00
    #[test]
    fn test_gmst_2000() {
        todo!();
    }
    #[test]
    fn test_gmst_2000_parity() {
        use rsofa::iauGmst00;
        todo!();
    }

    // t_sofa.c t_gmst06
    #[test]
    fn test_gmst_2006() {
        todo!();
    }
    #[test]
    fn test_gmst_2006_parity() {
        use rsofa::iauGmst06;
        todo!();
    }

    // t_sofa.c t_gmst82
    #[test]
    fn test_gmst_1982() {
        todo!();
    }
    #[test]
    fn test_gmst_1982_parity() {
        use rsofa::iauGmst82;
        todo!();
    }

    // t_sofa.c t_gst00a
    #[test]
    fn test_gst_2000a() {
        todo!();
    }
    #[test]
    fn test_gst_2000a_parity() {
        use rsofa::iauGst00a;
        todo!();
    }

    // t_sofa.c t_gst00b
    #[test]
    fn test_gst_2000b() {
        todo!();
    }
    #[test]
    fn test_gst_2000b_parity() {
        use rsofa::iauGst00b;
        todo!();
    }

    // t_sofa.c t_gst06
    #[test]
    fn test_gst_2006() {
        todo!();
    }
    #[test]
    fn test_gst_2006_parity() {
        use rsofa::iauGst06;
        todo!();
    }

    // t_sofa.c t_gst06a
    #[test]
    fn test_gst_2006a() {
        todo!();
    }
    #[test]
    fn test_gst_2006a_parity() {
        use rsofa::iauGst06a;
        todo!();
    }

    // t_sofa.c t_gst94
    #[test]
    fn test_gst_1994() {
        todo!();
    }
    #[test]
    fn test_gst_1994_parity() {
        use rsofa::iauGst94;
        todo!();
    }
}
