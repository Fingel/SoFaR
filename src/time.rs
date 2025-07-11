//!Time scales
//! - format_jd         D2DTF     format 2-part JD for output
//! - delta_at          DAT       Delta(AT) (=TAI-UTC) for a given UTC date
//! - tdb_tt_diff       DTDB      TDB-TT
//! - datetime_to_jd    DTF2D     encode time and date fields into 2-part JD
//! - tai_to_tt         TAITT     TAI to TT
//! - tai_to_ut1        TAIUT1    TAI to UT1
//! - tai_to_utc        TAIUTC    TAI to UTC
//! - tcb_to_tdb        TCBTDB    TCB to TDB
//! - tcg_to_tt         TCGTT     TCG to TT
//! - tdb_to_tcb        TDBTCB    TDB to TCB
//! - tdb_to_tt         TDBTT     TDB to TT
//! - tt_to_tai         TTTAI     TT to TAI
//! - tt_to_tcg         TTTCG     TT to TCG
//! - tt_to_tdb         TTTDB     TT to TDB
//! - tt_to_ut1         TTUT1     TT to UT1
//! - ut1_to_tai        UT1TAI    UT1 to TAI
//! - ut1_to_tt         UT1TT     UT1 to TT
//! - ut1_to_utc        UT1UTC    UT1 to UTC
//! - utc_to_tai        UTCTAI    UTC to TAI
//! - utc_to_ut1        UTCUT1    UTC to UT1
pub fn format_jd() {
    todo!();
}
pub fn delta_at() {
    todo!();
}
pub fn tdb_tt_diff() {
    todo!();
}
pub fn datetime_to_jd() {
    todo!();
}
pub fn tai_to_tt() {
    todo!();
}
pub fn tai_to_ut1() {
    todo!();
}
pub fn tai_to_utc() {
    todo!();
}
pub fn tcb_to_tdb() {
    todo!();
}
pub fn tcg_to_tt() {
    todo!();
}
pub fn tdb_to_tcb() {
    todo!();
}
pub fn tdb_to_tt() {
    todo!();
}
pub fn tt_to_tai() {
    todo!();
}
pub fn tt_to_tcg() {
    todo!();
}
pub fn tt_to_tdb() {
    todo!();
}
pub fn tt_to_ut1() {
    todo!();
}
pub fn ut1_to_tai() {
    todo!();
}
pub fn ut1_to_tt() {
    todo!();
}
pub fn ut1_to_utc() {
    todo!();
}
pub fn utc_to_tai() {
    todo!();
}
pub fn utc_to_ut1() {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;
    // t_sofa.c t_d2dtf
    #[test]
    fn test_format_jd() {
        todo!();
    }
    #[test]
    fn test_format_jd_parity() {
        use rsofa::iauD2dtf;
        todo!();
    }

    // t_sofa.c t_dat
    #[test]
    fn test_delta_at() {
        todo!();
    }
    #[test]
    fn test_delta_at_parity() {
        use rsofa::iauDat;
        todo!();
    }

    // t_sofa.c t_dtdb
    #[test]
    fn test_tdb_tt_diff() {
        todo!();
    }
    #[test]
    fn test_tdb_tt_diff_parity() {
        use rsofa::iauDtdb;
        todo!();
    }

    // t_sofa.c t_dtf2d
    #[test]
    fn test_datetime_to_jd() {
        todo!();
    }
    #[test]
    fn test_datetime_to_jd_parity() {
        use rsofa::iauDtf2d;
        todo!();
    }

    // t_sofa.c t_taitt
    #[test]
    fn test_tai_to_tt() {
        todo!();
    }
    #[test]
    fn test_tai_to_tt_parity() {
        use rsofa::iauTaitt;
        todo!();
    }

    // t_sofa.c t_taiut1
    #[test]
    fn test_tai_to_ut1() {
        todo!();
    }
    #[test]
    fn test_tai_to_ut1_parity() {
        use rsofa::iauTaiut1;
        todo!();
    }

    // t_sofa.c t_taiutc
    #[test]
    fn test_tai_to_utc() {
        todo!();
    }
    #[test]
    fn test_tai_to_utc_parity() {
        use rsofa::iauTaiutc;
        todo!();
    }

    // t_sofa.c t_tcbtdb
    #[test]
    fn test_tcb_to_tdb() {
        todo!();
    }
    #[test]
    fn test_tcb_to_tdb_parity() {
        use rsofa::iauTcbtdb;
        todo!();
    }

    // t_sofa.c t_tcgtt
    #[test]
    fn test_tcg_to_tt() {
        todo!();
    }
    #[test]
    fn test_tcg_to_tt_parity() {
        use rsofa::iauTcgtt;
        todo!();
    }

    // t_sofa.c t_tdbtcb
    #[test]
    fn test_tdb_to_tcb() {
        todo!();
    }
    #[test]
    fn test_tdb_to_tcb_parity() {
        use rsofa::iauTdbtcb;
        todo!();
    }

    // t_sofa.c t_tdbtt
    #[test]
    fn test_tdb_to_tt() {
        todo!();
    }
    #[test]
    fn test_tdb_to_tt_parity() {
        use rsofa::iauTdbtt;
        todo!();
    }

    // t_sofa.c t_tttai
    #[test]
    fn test_tt_to_tai() {
        todo!();
    }
    #[test]
    fn test_tt_to_tai_parity() {
        use rsofa::iauTttai;
        todo!();
    }

    // t_sofa.c t_tttcg
    #[test]
    fn test_tt_to_tcg() {
        todo!();
    }
    #[test]
    fn test_tt_to_tcg_parity() {
        use rsofa::iauTttcg;
        todo!();
    }

    // t_sofa.c t_tttdb
    #[test]
    fn test_tt_to_tdb() {
        todo!();
    }
    #[test]
    fn test_tt_to_tdb_parity() {
        use rsofa::iauTttdb;
        todo!();
    }

    // t_sofa.c t_ttut1
    #[test]
    fn test_tt_to_ut1() {
        todo!();
    }
    #[test]
    fn test_tt_to_ut1_parity() {
        use rsofa::iauTtut1;
        todo!();
    }

    // t_sofa.c t_ut1tai
    #[test]
    fn test_ut1_to_tai() {
        todo!();
    }
    #[test]
    fn test_ut1_to_tai_parity() {
        use rsofa::iauUt1tai;
        todo!();
    }

    // t_sofa.c t_ut1tt
    #[test]
    fn test_ut1_to_tt() {
        todo!();
    }
    #[test]
    fn test_ut1_to_tt_parity() {
        use rsofa::iauUt1tt;
        todo!();
    }

    // t_sofa.c t_ut1utc
    #[test]
    fn test_ut1_to_utc() {
        todo!();
    }
    #[test]
    fn test_ut1_to_utc_parity() {
        use rsofa::iauUt1utc;
        todo!();
    }

    // t_sofa.c t_utctai
    #[test]
    fn test_utc_to_tai() {
        todo!();
    }
    #[test]
    fn test_utc_to_tai_parity() {
        use rsofa::iauUtctai;
        todo!();
    }

    // t_sofa.c t_utcut1
    #[test]
    fn test_utc_to_ut1() {
        todo!();
    }
    #[test]
    fn test_utc_to_ut1_parity() {
        use rsofa::iauUtcut1;
        todo!();
    }
}
