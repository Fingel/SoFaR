//!  Galactic coordinates
//! - galactic_to_icrs      G2ICRS    transform IAU 1958 galactic coordinates to ICRS
//! - icrs_to_galactic      ICRS2G    transform ICRS coordinates to IAU 1958 Galactic

pub fn galactic_to_icrs() {
    todo!();
}
pub fn icrs_to_galactic() {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;

    // t_sofa.c t_g2icrs
    #[test]
    fn test_galactic_to_icrs() {
        todo!()
    }
    #[test]
    fn test_galactic_to_icrs_parity() {
        use rsofa::iauG2icrs;
        todo!()
    }

    // t_sofa.c t_icrs2g
    #[test]
    fn test_icrs_to_galactic() {
        todo!()
    }
    #[test]
    fn test_icrs_to_galactic_parity() {
        use rsofa::iauIcrs2g;
        todo!()
    }
}
