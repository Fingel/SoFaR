//! Operations on angles

///   To sexagesimal
///   - radians_to_hms (A2TF)      decompose radians into hours, minutes, seconds
///   - A2AF      decompose radians into degrees, arcminutes, arcseconds
///   - days_to_hms (D2TF)      decompose days into hours, minutes, seconds
pub mod to_sexagesimal {
    use crate::constants::*;
    use crate::{dint, dnint};
    ///  Decompose radians into hours, minutes, seconds, fraction.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     ndp     int     resolution (Note 1)
    ///     angle   double  angle in radians
    ///
    ///  Returned:
    ///     sign    char*   '+' or '-'
    ///     ihmsf   int[4]  hours, minutes, seconds, fraction
    ///
    ///  Notes:
    ///
    ///  1) The argument ndp is interpreted as follows:
    ///
    ///  ndp         resolution
    ///   :      ...0000 00 00
    ///  -7         1000 00 00
    ///  -6          100 00 00
    ///  -5           10 00 00
    ///  -4            1 00 00
    ///  -3            0 10 00
    ///  -2            0 01 00
    ///  -1            0 00 10
    ///   0            0 00 01
    ///   1            0 00 00.1
    ///   2            0 00 00.01
    ///   3            0 00 00.001
    ///   :            0 00 00.000...
    ///
    ///  2) The largest positive useful value for ndp is determined by the
    ///     size of angle, the format of doubles on the target platform, and
    ///     the risk of overflowing ihmsf[3].  On a typical platform, for
    ///     angle up to 2pi, the available floating-point precision might
    ///     correspond to ndp=12.  However, the practical limit is typically
    ///     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
    ///     only 16 bits.
    ///
    ///  3) The absolute value of angle may exceed 2pi.  In cases where it
    ///     does not, it is up to the caller to test for and handle the
    ///     case where angle is very nearly 2pi and rounds up to 24 hours,
    ///     by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
    ///
    ///  Called:
    ///     iauD2tf      decompose days to hms
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn radians_to_hms(ndp: i32, angle: f64) -> (char, [i32; 4]) {
        // Scale then use days to h,m,s function.
        days_to_hms(ndp, angle / D2PI)
    }

    //TODO: optimize
    ///  Decompose days to hours, minutes, seconds, fraction.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     ndp     int     resolution (Note 1)
    ///     days    double  interval in days
    ///
    ///  Returned:
    ///     sign    char*   '+' or '-'
    ///     ihmsf   int[4]  hours, minutes, seconds, fraction
    ///
    ///  Notes:
    ///
    ///  1) The argument ndp is interpreted as follows:
    ///
    ///  ndp         resolution
    ///   :      ...0000 00 00
    ///  -7         1000 00 00
    ///  -6          100 00 00
    ///  -5           10 00 00
    ///  -4            1 00 00
    ///  -3            0 10 00
    ///  -2            0 01 00
    ///  -1            0 00 10
    ///   0            0 00 01
    ///   1            0 00 00.1
    ///   2            0 00 00.01
    ///   3            0 00 00.001
    ///   :            0 00 00.000...
    ///
    ///  2) The largest positive useful value for ndp is determined by the
    ///     size of days, the format of double on the target platform, and
    ///     the risk of overflowing ihmsf[3].  On a typical platform, for
    ///     days up to 1.0, the available floating-point precision might
    ///     correspond to ndp=12.  However, the practical limit is typically
    ///     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
    ///     only 16 bits.
    ///
    ///  3) The absolute value of days may exceed 1.0.  In cases where it
    ///     does not, it is up to the caller to test for and handle the
    ///     case where days is very nearly 1.0 and rounds up to 24 hours,
    ///     by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    ///
    ///  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
    pub fn days_to_hms(ndp: i32, days: f64) -> (char, [i32; 4]) {
        let mut nrs: i32;
        let mut rs: f64;

        // Handle sign.
        let sign = if days >= 0.0 { '+' } else { '-' };

        // Interval in seconds.
        let mut a = DAYSEC * days.abs();

        // Pre-round  if resolution coarser than 1s (then pretend ndp=1).
        if ndp < 0 {
            nrs = 1;
            for n in 1..=-ndp {
                nrs *= if n == 2 || n == 4 { 6 } else { 10 };
            }
            rs = nrs as f64;
            let w = a / rs;
            a = rs * dnint(w);
        }

        // Express the unit of each field in resolution units.
        nrs = 1;
        for _ in 1..=ndp {
            //TODO: Optimize
            nrs *= 10;
        }
        rs = nrs as f64;
        let rm = rs * 60.0;
        let rh = rm * 60.0;

        // Round the interval and express in resolution units.
        a = dnint(rs * a);

        // Break into fields
        //TODO: no need for mut
        let mut ah = a / rh;
        ah = dint(ah);
        a -= ah * rh;
        let mut am = a / rm;
        am = dint(am);
        a -= am * rm;
        let mut r#as = a / rs;
        r#as = dint(r#as);
        let af = a - r#as * rs;

        // Return results
        (sign, [ah as i32, am as i32, r#as as i32, af as i32])
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// t_sofa.c t_a2tf
        #[test]
        fn test_a2tf() {
            let result = radians_to_hms(4, -3.01234);
            assert_eq!(result, ('-', [11, 30, 22, 6484]));
        }

        /// t_sofa.c t_d2tf
        #[test]
        fn test_d2tf() {
            let result = days_to_hms(4, -0.987654321);
            assert_eq!(result, ('-', [23, 42, 13, 3333]));
        }
    }
}
