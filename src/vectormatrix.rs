use crate::{DAYSEC, dint, dnint};

///  - - - - - - - -
///   i a u D 2 t f
///  - - - - - - - -
///
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
pub fn iau_d2tf(ndp: i32, days: f64) -> (char, [i32; 4]) {
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
        nrs *= 10;
    }
    rs = nrs as f64;
    let rm = rs * 60.0;
    let rh = rm * 60.0;

    // Round the interval and express in resolution units.
    a = dnint(rs * a);

    // Break into fields
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

pub use iau_d2tf as days_to_hours_minutes_seconds_fraction;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_d2tf() {
        let result = iau_d2tf(4, -0.987654321);
        assert_eq!(result, ('-', [23, 42, 13, 3333]));
    }
}
