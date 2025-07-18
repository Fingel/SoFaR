//! Star catalog conversions
//! fk5_to_hip              FK52H     transform FK5 star data into the Hipparcos system
//! fk5_to_hipp_rot         FK5HIP    FK5 to Hipparcos rotation and spin
//! fk5_to_hipp_zero        FK5HZ     FK5 to Hipparcos assuming zero Hipparcos proper motion
//! hipp_to_fk5             H2FK5     transform Hipparcos star data into the FK5 system
//! hipp_to_fk5_zero        HFK5Z     Hipparcos to FK5 assuming zero Hipparcos proper motion
//! fk4_to_fk5              FK425     transform FK4 star data into FK5
//! fk4_to_fk5_zero         FK45Z     FK4 to FK5 assuming zero FK5 proper motion
//! fk5_to_fk4              FK524     transform FK5 star data into FK4
//! fk5_to_fk4_zero         FK54Z     FK5 to FK4 assuming zero FK5 proper motion

pub fn fk5_to_hip() {
    todo!();
}
pub fn fk5_to_hipp_rot() {
    todo!();
}
pub fn fk5_to_hipp_zero() {
    todo!();
}
pub fn hipp_to_fk5() {
    todo!();
}
pub fn hipp_to_fk5_zero() {
    todo!();
}
pub fn fk4_to_fk5() {
    todo!();
}
pub fn fk4_to_fk5_zero() {
    todo!();
}
pub fn fk5_to_fk4() {
    todo!();
}
pub fn fk5_to_fk4_zero() {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;

    //t_sofa.c t_fk52h
    #[test]
    fn test_fk5_to_hip() {
        todo!();
    }

    #[test]
    fn test_fk5_to_hip_parity() {
        use rsofa::iauFk52h;
        todo!();
    }

    //t_sofa.c t_fk5hip
    #[test]
    fn test_fk5_to_hipp_rot() {
        todo!();
    }

    #[test]
    fn test_fk5_to_hipp_rot_parity() {
        use rsofa::iauFk5hip;
        todo!();
    }

    //t_sofa.c t_fk5hz
    #[test]
    fn test_fk5_to_hipp_zero() {
        todo!();
    }

    #[test]
    fn test_fk5_to_hipp_zero_parity() {
        use rsofa::iauFk5hz;
        todo!();
    }

    //t_sofa.c t_h2fk5
    #[test]
    fn test_hipp_to_fk5() {
        todo!();
    }

    #[test]
    fn test_hipp_to_fk5_parity() {
        use rsofa::iauH2fk5;
        todo!();
    }

    //t_sofa.c t_hfk5z
    #[test]
    fn test_hipp_to_fk5_zero() {
        todo!();
    }

    #[test]
    fn test_hipp_to_fk5_zero_parity() {
        use rsofa::iauHfk5z;
        todo!();
    }

    //t_sofa.c t_fk425
    #[test]
    fn test_fk4_to_fk5() {
        todo!();
    }

    #[test]
    fn test_fk4_to_fk5_parity() {
        use rsofa::iauFk425;
        todo!();
    }

    //t_sofa.c t_fk45z
    #[test]
    fn test_fk4_to_fk5_zero() {
        todo!();
    }

    #[test]
    fn test_fk4_to_fk5_zero_parity() {
        use rsofa::iauFk45z;
        todo!();
    }

    //t_sofa.c t_fk524
    #[test]
    fn test_fk5_to_fk4() {
        todo!();
    }

    #[test]
    fn test_fk5_to_fk4_parity() {
        use rsofa::iauFk524;
        todo!();
    }

    //t_sofa.c t_fk54z
    #[test]
    fn test_fk5_to_fk4_zero() {
        todo!();
    }

    #[test]
    fn test_fk5_to_fk4_zero_parity() {
        use rsofa::iauFk54z;
        todo!();
    }
}
