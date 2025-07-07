//! Operations on angles

// Wrap
// - angle_normalize_positive   ANP       normalize radians to range 0 to 2pi
// - angle_normalize_pm         ANPM      normalize radians to range -pi to +pi
pub mod wrap {

    use crate::constants::{D2PI, DPI};
    ///  Normalize angle into the range 0 <= a < 2pi.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double     angle (radians)
    ///
    ///  Returned (function value):
    ///              double     angle in range 0-2pi
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn angle_normalize_positive(a: f64) -> f64 {
        let w = a % D2PI;
        if w < 0.0 { w + D2PI } else { w }
    }

    ///  Normalize angle into the range -pi <= a < +pi.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  vector/matrix support function.
    ///
    ///  Given:
    ///     a        double     angle (radians)
    ///
    ///  Returned (function value):
    ///              double     angle in range +/-pi
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn angle_normalize_pm(a: f64) -> f64 {
        let w = a % D2PI;
        if w.abs() >= DPI {
            w - D2PI.copysign(a)
        } else {
            w
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use assert_approx_eq::assert_approx_eq;

        // t_sofa.c t_anp
        #[test]
        fn test_anp() {
            assert_approx_eq!(angle_normalize_positive(-0.1), 6.183185307179586477, 1e-12);
        }

        #[test]
        fn test_anp_parity() {
            use rsofa::iauAnp;
            let r = angle_normalize_positive(-0.1);
            let r_iau = unsafe { iauAnp(-0.1) };
            assert_eq!(r, r_iau);
        }

        // t_sofa.c t_anpm
        #[test]
        fn test_anpm() {
            assert_approx_eq!(angle_normalize_pm(-4.0), 2.283185307179586477, 1e-12);
        }

        #[test]
        fn test_anpm_parity() {
            use rsofa::iauAnpm;
            let r = angle_normalize_pm(4.0);
            let r_iau = unsafe { iauAnpm(4.0) };
            assert_eq!(r, r_iau);
        }
    }
}

///   To sexagesimal
///   - radians_to_hms  (A2TF)      decompose radians into hours, minutes, seconds
///   - radians_to_dms  (A2AF)      decompose radians into degrees, arcminutes, arcseconds
///   - days_to_hms     (D2TF)      decompose days into hours, minutes, seconds
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
    pub fn radians_to_hms(ndp: i32, angle: f64) -> (char, [i32; 4]) {
        // Scale then use days to h,m,s function.
        days_to_hms(ndp, angle / D2PI)
    }

    ///Decompose radians into degrees, arcminutes, arcseconds, fraction.
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
    ///     idmsf   int[4]  degrees, arcminutes, arcseconds, fraction
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
    ///     the risk of overflowing idmsf[3].  On a typical platform, for
    ///     angle up to 2pi, the available floating-point precision might
    ///     correspond to ndp=12.  However, the practical limit is typically
    ///     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
    ///     only 16 bits.
    ///
    ///  3) The absolute value of angle may exceed 2pi.  In cases where it
    ///     does not, it is up to the caller to test for and handle the
    ///     case where angle is very nearly 2pi and rounds up to 360 degrees,
    ///     by testing for idmsf[0]=360 and setting idmsf[0-3] to zero.
    ///
    ///  Called:
    ///     iauD2tf      decompose days to hms
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn radians_to_dms(ndp: i32, angle: f64) -> (char, [i32; 4]) {
        let f = 15.0 / D2PI;
        days_to_hms(ndp, angle * f)
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
        use std::os::raw::c_char;

        /// t_sofa.c t_a2tf
        #[test]
        fn test_a2tf() {
            let (sign, hmsf) = radians_to_hms(4, -3.01234);
            assert_eq!(sign, '-');
            assert_eq!(hmsf, [11, 30, 22, 6484]);
        }

        #[test]
        fn test_a2tf_parity() {
            use rsofa::iauA2tf;
            let (sign, hmsf) = radians_to_hms(4, -3.01234);
            let mut sign_iau: c_char = 0;
            let mut hmsf_iau = [0; 4];

            unsafe {
                iauA2tf(4, -3.01234, &mut sign_iau, hmsf_iau.as_mut_ptr());
            }
            assert_eq!(sign, sign_iau as u8 as char);
            assert_eq!(hmsf, hmsf_iau);
        }

        /// t_sofa.c t_a2af
        #[test]
        fn test_a2af() {
            let (sign, dmsf) = radians_to_dms(4, 2.345);
            assert_eq!(sign, '+');
            assert_eq!(dmsf, [134, 21, 30, 9706]);
        }

        #[test]
        fn test_a2af_parity() {
            use rsofa::iauA2af;
            let (sign, hmsf) = radians_to_dms(4, 2.345);
            let mut sign_iau: c_char = 0;
            let mut hmsf_iau = [0; 4];
            unsafe {
                iauA2af(4, 2.345, &mut sign_iau, hmsf_iau.as_mut_ptr());
            }
            assert_eq!(sign, sign_iau as u8 as char);
            assert_eq!(hmsf, hmsf_iau);
        }

        /// t_sofa.c t_d2tf
        #[test]
        fn test_d2tf() {
            let (sign, hmsf) = days_to_hms(4, -0.987654321);
            assert_eq!(sign, '-');
            assert_eq!(hmsf, [23, 42, 13, 3333]);
        }

        #[test]
        fn test_d2tf_parity() {
            use rsofa::iauD2tf;
            let (sign, hmsf) = days_to_hms(4, -0.987654321);
            let mut sign_iau: c_char = 0;
            let mut hmsf_iau = [0; 4];

            unsafe {
                iauD2tf(4, -0.987654321, &mut sign_iau, hmsf_iau.as_mut_ptr());
            }
            assert_eq!(sign, sign_iau as u8 as char);
            assert_eq!(hmsf, hmsf_iau);
        }
    }
}

///From sexagesimal
///
/// - dms_to_radians    AF2A      degrees, arcminutes, arcseconds to radians
/// - hms_to_radians    TF2A      hours, minutes, seconds to radians
/// - hms_to_days       TF2D      hours, minutes, seconds to days
pub mod from_sexagesimal {
    use crate::constants::{DAS2R, DAYSEC, DS2R};
    #[derive(Debug, PartialEq)]
    pub enum DmsStatus {
        Ok,
        /// ideg outside 0-359
        DegOutOfRange,
        /// iamin outside 0-59
        AminOutOfRange,
        /// asec outside 0-59.999...
        AsecOutOfRange,
    }

    impl From<i32> for DmsStatus {
        fn from(status_code: i32) -> Self {
            match status_code {
                0 => DmsStatus::Ok,
                1 => DmsStatus::DegOutOfRange,
                2 => DmsStatus::AminOutOfRange,
                3 => DmsStatus::AsecOutOfRange,
                _ => DmsStatus::Ok, // Default to Ok for unknown values
            }
        }
    }
    ///  Convert degrees, arcminutes, arcseconds to radians.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  support function.
    ///
    ///  Given:
    ///     s         char    sign:  '-' = negative, otherwise positive
    ///     ideg      int     degrees
    ///     iamin     int     arcminutes
    ///     asec      double  arcseconds
    ///
    ///  Returned:
    ///     rad       double  angle in radians
    ///
    ///  Returned (function value):
    ///               int     status:  0 = OK
    ///                                1 = ideg outside range 0-359
    ///                                2 = iamin outside range 0-59
    ///                                3 = asec outside range 0-59.999...
    ///
    ///  Notes:
    ///
    ///  1)  The result is computed even if any of the range checks fail.
    ///
    ///  2)  Negative ideg, iamin and/or asec produce a warning status, but
    ///      the absolute value is used in the conversion.
    ///
    ///  3)  If there are multiple errors, the status value reflects only the
    ///      first, the smallest taking precedence.
    ///
    ///  This revision:  2021 May 11
    ///
    ///  SOFA release 2023-10-11
    pub fn dms_to_radians(s: char, ideg: i32, iamin: i32, asec: f64) -> (f64, DmsStatus) {
        let sign = match s {
            '-' => -1.0,
            _ => 1.0,
        };

        // Compute the interval.
        let rad =
            sign * (60.0 * (60.0 * ideg.abs() as f64 + iamin.abs() as f64) + asec.abs()) * DAS2R;

        // Validate arguments and return status.
        if !(0..=359).contains(&ideg) {
            (rad, DmsStatus::DegOutOfRange)
        } else if !(0..=59).contains(&iamin) {
            (rad, DmsStatus::AminOutOfRange)
        } else if !(0.0..60.0).contains(&asec) {
            (rad, DmsStatus::AsecOutOfRange)
        } else {
            (rad, DmsStatus::Ok)
        }
    }

    #[derive(Debug, PartialEq)]
    pub enum HmsStatus {
        Ok,
        /// ihour outside range 0-23
        IhourOutOfRange,
        /// imin outside range 0-59
        IminOutOfRange,
        /// sec outside range 0-59.999...
        SecOutOfRange,
    }

    impl From<i32> for HmsStatus {
        fn from(status_code: i32) -> Self {
            match status_code {
                0 => HmsStatus::Ok,
                1 => HmsStatus::IhourOutOfRange,
                2 => HmsStatus::IminOutOfRange,
                3 => HmsStatus::SecOutOfRange,
                _ => HmsStatus::Ok, // Default to Ok for unknown values
            }
        }
    }
    ///  Convert hours, minutes, seconds to radians.
    ///
    ///  This function is part of the International Astronomical Union's
    ///  SOFA (Standards of Fundamental Astronomy) software collection.
    ///
    ///  Status:  support function.
    ///
    ///  Given:
    ///     s         char    sign:  '-' = negative, otherwise positive
    ///     ihour     int     hours
    ///     imin      int     minutes
    ///     sec       double  seconds
    ///
    ///  Returned:
    ///     rad       double  angle in radians
    ///
    ///  Returned (function value):
    ///               int     status:  0 = OK
    ///                                1 = ihour outside range 0-23
    ///                                2 = imin outside range 0-59
    ///                                3 = sec outside range 0-59.999...
    ///
    ///  Notes:
    ///
    ///  1)  The result is computed even if any of the range checks fail.
    ///
    ///  2)  Negative ihour, imin and/or sec produce a warning status, but
    ///      the absolute value is used in the conversion.
    ///
    ///  3)  If there are multiple errors, the status value reflects only the
    ///      first, the smallest taking precedence.
    ///
    ///  This revision:  2021 May 11
    pub fn hms_to_radians(s: char, ihour: i32, imin: i32, sec: f64) -> (f64, HmsStatus) {
        let sign = match s {
            '-' => -1.0,
            _ => 1.0,
        };

        // Compute the interval.
        let rad =
            sign * (60.0 * (60.0 * ihour.abs() as f64 + imin.abs() as f64) + sec.abs()) * DS2R;

        // Validate arguments and return status.
        if !(0..=23).contains(&ihour) {
            (rad, HmsStatus::IhourOutOfRange)
        } else if !(0..=59).contains(&imin) {
            (rad, HmsStatus::IminOutOfRange)
        } else if !(0.0..60.0).contains(&sec) {
            (rad, HmsStatus::SecOutOfRange)
        } else {
            (rad, HmsStatus::Ok)
        }
    }

    //  Convert hours, minutes, seconds to days.
    //
    //  This function is part of the International Astronomical Union's
    //  SOFA (Standards of Fundamental Astronomy) software collection.
    //
    //  Status:  support function.
    //
    //  Given:
    //     s         char    sign:  '-' = negative, otherwise positive
    //     ihour     int     hours
    //     imin      int     minutes
    //     sec       double  seconds
    //
    //  Returned:
    //     days      double  interval in days
    //
    //  Returned (function value):
    //               int     status:  0 = OK
    //                                1 = ihour outside range 0-23
    //                                2 = imin outside range 0-59
    //                                3 = sec outside range 0-59.999...
    //
    //  Notes:
    //
    //  1)  The result is computed even if any of the range checks fail.
    //
    //  2)  Negative ihour, imin and/or sec produce a warning status, but
    //      the absolute value is used in the conversion.
    //
    //  3)  If there are multiple errors, the status value reflects only the
    //      first, the smallest taking precedence.
    //
    //  This revision:  2021 May 11
    //
    //  SOFA release 2023-10-11
    pub fn hms_to_days(s: char, ihour: i32, imin: i32, sec: f64) -> (f64, HmsStatus) {
        let sign = match s {
            '-' => -1.0,
            _ => 1.0,
        };
        // Compute the interval.
        let days =
            sign * (60.0 * (60.0 * ihour.abs() as f64 + imin.abs() as f64) + sec.abs()) / DAYSEC;

        // Validate arguments and return status.
        if !(0..=23).contains(&ihour) {
            (days, HmsStatus::IhourOutOfRange)
        } else if !(0..=59).contains(&imin) {
            (days, HmsStatus::IminOutOfRange)
        } else if !(0.0..60.0).contains(&sec) {
            (days, HmsStatus::SecOutOfRange)
        } else {
            (days, HmsStatus::Ok)
        }
    }

    #[cfg(test)]
    mod tests {
        use std::os::raw::c_char;

        use assert_approx_eq::assert_approx_eq;

        use super::*;

        // t_sofa.c t_af2a
        #[test]
        fn test_af2a() {
            let (a, status) = dms_to_radians('-', 45, 13, 27.2);
            assert_approx_eq!(a, -0.7893115794313644842, 1e-12);
            assert_eq!(status, DmsStatus::Ok);
        }

        #[test]
        fn test_af2a_status() {
            let (_, status) = dms_to_radians('-', 45, 13, 60.0);
            assert_eq!(status, DmsStatus::AsecOutOfRange);
        }

        #[test]
        fn test_af2a_parity() {
            use rsofa::iauAf2a;
            let (a, status) = dms_to_radians('-', 45, 13, 27.2);

            let mut a_iau = 0.0;
            let status_iau = unsafe { iauAf2a('-' as c_char, 45, 13, 27.2, &mut a_iau) };
            assert_eq!(a, a_iau);
            assert_eq!(status, DmsStatus::from(status_iau));
        }

        // t_sofa.c t_tf2a
        #[test]
        fn test_tf2a() {
            let (a, status) = hms_to_radians('+', 4, 58, 20.2);
            assert_approx_eq!(a, 1.301739278189537429, 1e-12);
            assert_eq!(status, HmsStatus::Ok);
        }

        #[test]
        fn test_tf2a_parity() {
            use rsofa::iauTf2a;
            let (a, status) = hms_to_radians('+', 4, 58, 20.2);

            let mut a_iau = 0.0;
            let status_iau = unsafe { iauTf2a('+' as c_char, 4, 58, 20.2, &mut a_iau) };
            assert_eq!(a, a_iau);
            assert_eq!(status, HmsStatus::from(status_iau));
        }

        #[test]
        fn test_tf2a_status() {
            let (_, status) = hms_to_radians('+', 4, 120, 20.2);
            assert_eq!(status, HmsStatus::IminOutOfRange);
        }

        // t_sofa.c t_tf2d
        #[test]
        fn test_tf2d() {
            let (d, status) = hms_to_days(' ', 23, 55, 10.9);
            assert_approx_eq!(d, 0.9966539351851851852, 1e-12);
            assert_eq!(status, HmsStatus::Ok);
        }

        #[test]
        fn test_tf2d_parity() {
            use rsofa::iauTf2d;
            let (d, status) = hms_to_days(' ', 23, 55, 10.9);

            let mut d_iau = 0.0;
            let status_iau = unsafe { iauTf2d(' ' as c_char, 23, 55, 10.9, &mut d_iau) };
            assert_eq!(d, d_iau);
            assert_eq!(status, HmsStatus::from(status_iau));
        }
    }
}
