//!Ephemerides (limited precision)
//!
//! - earth_pv   EPV00     Earth position and velocity
//! - moon_pv    MOON98    Moon position and velocity
//! - planet_pv  PLAN94    major-planet position and velocity
pub fn earth_pv() {
    todo!();
}
pub fn moon_pv() {
    todo!();
}
pub fn planet_pv() {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;

    //t_sofa.c t_epv00
    #[test]
    fn test_earth_pv() {
        todo!();
    }
    #[test]
    fn test_earth_pv_parity() {
        use rsofa::iauEpv00;
        todo!();
    }

    //t_sofa.c t_moon98
    #[test]
    fn test_moon_pv() {
        todo!();
    }
    #[test]
    fn test_moon_pv_parity() {
        use rsofa::iauMoon98;
        todo!();
    }

    //t_sofa.c t_plan94
    #[test]
    fn test_planet_pv() {
        todo!();
    }
    #[test]
    fn test_planet_pv_parity() {
        use rsofa::iauPlan94;
        todo!();
    }
}
