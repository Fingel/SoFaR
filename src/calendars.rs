//! Calendars
//! - cal_to_jd                 CAL2JD    Gregorian calendar to Julian Day number
//! - jd_to_besselian_e         EPB       Julian Date to Besselian Epoch
//! - besselian_e_to_jd         EPB2JD    Besselian Epoch to Julian Date
//! - jd_to_je                  EPJ       Julian Date to Julian Epoch
//! - je_to_jd                  EPJ2JD    Julian Epoch to Julian Date
//! - jd_to_cal                 JD2CAL    Julian Date to Gregorian year, month, day, fraction
//! - jd_to_cal_fmt             JDCALF    Julian Date to Gregorian date for formatted output
use crate::{Warned, Warning, constants::*, dnint, warning};
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum DateError {
    #[error("Bad year")]
    InvalidYear,
    #[error("Bad month")]
    InvalidMonth,
}
///  Gregorian Calendar to Julian Date.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  support function.
///
///  Given:
///     iy,im,id  int     year, month, day in Gregorian calendar (Note 1)
///
///  Returned:
///     djm0      double  MJD zero-point: always 2400000.5
///     djm       double  Modified Julian Date for 0 hrs
///
///  Returned (function value):
///               int     status:
///                           0 = OK
///                          -1 = bad year   (Note 3: JD not computed)
///                          -2 = bad month  (JD not computed)
///                          -3 = bad day    (JD computed)
///
///  Notes:
///
///  1) The algorithm used is valid from -4800 March 1, but this
///     implementation rejects dates before -4799 January 1.
///
///  2) The Julian Date is returned in two pieces, in the usual SOFA
///     manner, which is designed to preserve time resolution.  The
///     Julian Date is available as a single number by adding djm0 and
///     djm.
///
///  3) In early eras the conversion is from the "Proleptic Gregorian
///     Calendar";  no account is taken of the date(s) of adoption of
///     the Gregorian Calendar, nor is the AD/BC numbering convention
///     observed.
///
///  Reference:
///
///     Explanatory Supplement to the Astronomical Almanac,
///     P. Kenneth Seidelmann (ed), University Science Books (1992),
///     Section 12.92 (p604).
///
///  This revision:  2021 May 11
///
///  SOFA release 2023-10-11
pub fn cal_to_jd(iy: i32, im: i32, id: i32) -> Result<Warned<(f64, f64)>, DateError> {
    // Earliest year allowed (4800BC)
    const IYMIN: i32 = -4799;

    // Month lengths in days
    const MTAB: [i32; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    // Validate year and month.
    if iy < IYMIN {
        return Err(DateError::InvalidYear);
    }
    if !(1..=12).contains(&im) {
        return Err(DateError::InvalidMonth);
    }

    // If February in a leap year, 1, otherwise 0.
    let ly = if (im == 2) && (iy % 4 == 0) && (iy % 100 != 0 || (iy % 400 == 0)) {
        1
    } else {
        0
    };

    let my = (im - 14) / 12;
    let iypmy = (iy + my) as i64;
    let djm0 = DJM0;
    let djm = ((1461_i64 * (iypmy + 4800_i64)) / 4_i64
        + (367_i64 * (im - 2 - 12 * my) as i64) / 12_i64
        - (3_i64 * ((iypmy + 4900_i64) / 100_i64)) / 4_i64
        + id as i64
        - 2432076_i64) as f64;

    // Validate day, taking into account leap years.
    if (id < 1) || (id > (MTAB[im as usize - 1] + ly)) {
        Ok(Warned::with_warning((djm0, djm), warning!(-3, "bad day")))
    } else {
        Ok(Warned::new((djm0, djm)))
    }
}

///  Julian Date to Besselian Epoch.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  support function.
///
///  Given:
///     dj1,dj2    double     Julian Date (Notes 3,4)
///
///  Returned (function value):
///                double     Besselian Epoch.
///
///  Notes:
///
///  1) Besselian Epoch is a method of expressing a moment in time as a
///     year plus fraction.  It was superseded by Julian Year (see the
///     function iauEpj).
///
///  2) The start of a Besselian year is when the right ascension of
///     the fictitious mean Sun is 18h 40m, and the unit is the tropical
///     year.  The conventional definition (see Lieske 1979) is that
///     Besselian Epoch B1900.0 is JD 2415020.31352 and the length of the
///     year is 365.242198781 days.
///
///  3) The time scale for the JD, originally Ephemeris Time, is TDB,
///     which for all practical purposes in the present context is
///     indistinguishable from TT.
///
///  4) The Julian Date is supplied in two pieces, in the usual SOFA
///     manner, which is designed to preserve time resolution.  The
///     Julian Date is available as a single number by adding dj1 and
///     dj2.  The maximum resolution is achieved if dj1 is 2451545.0
///     (J2000.0).
///
///  Reference:
///
///     Lieske, J.H., 1979. Astron.Astrophys., 73, 282.
///
///  This revision:  2023 May 5
///
///  SOFA release 2023-10-11
pub fn jd_to_besselian_e(dj1: f64, dj2: f64) -> f64 {
    const D1900: f64 = 36524.68648;

    1900.0 + ((dj1 - DJ00) + (dj2 + D1900)) / DTY
}

///  Besselian Epoch to Julian Date.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  support function.
///
///  Given:
///     epb      double    Besselian Epoch (e.g. 1957.3)
///
///  Returned:
///     djm0     double    MJD zero-point: always 2400000.5
///     djm      double    Modified Julian Date
///
///  Note:
///
///     The Julian Date is returned in two pieces, in the usual SOFA
///     manner, which is designed to preserve time resolution.  The
///     Julian Date is available as a single number by adding djm0 and
///     djm.
///
///  Reference:
///
///     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
///
///  This revision:  2021 May 11
///
///  SOFA release 2023-10-11
pub fn besselian_e_to_jd(epb: f64) -> (f64, f64) {
    let djm = 15019.81352 + (epb - 1900.0) * DTY;
    (DJM0, djm)
}

///  Julian Date to Julian Epoch.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  support function.
///
///  Given:
///     dj1,dj2    double     Julian Date (Note 4)
///
///  Returned (function value):
///                double     Julian Epoch
///
///  Notes:
///
///  1) Julian Epoch is a method of expressing a moment in time as a
///     year plus fraction.
///
///  2) Julian Epoch J2000.0 is 2000 Jan 1.5, and the length of the year
///     is 365.25 days.
///
///  3) For historical reasons, the time scale formally associated with
///     Julian Epoch is TDB (or TT, near enough).  However, Julian Epoch
///     can be used more generally as a calendrical convention to
///     represent other time scales such as TAI and TCB.  This is
///     analogous to Julian Date, which was originally defined
///     specifically as a way of representing Universal Times but is now
///     routinely used for any of the regular time scales.
///
///  4) The Julian Date is supplied in two pieces, in the usual SOFA
///     manner, which is designed to preserve time resolution.  The
///     Julian Date is available as a single number by adding dj1 and
///     dj2.  The maximum resolution is achieved if dj1 is 2451545.0
///     (J2000.0).
///
///  Reference:
///
///     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
///
///  This revision:  2022 May 6
///
///  SOFA release 2023-10-11
pub fn jd_to_je(dj1: f64, dj2: f64) -> f64 {
    2000.0 + ((dj1 - DJ00) + dj2) / DJY
}

///  Julian Epoch to Julian Date.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  support function.
///
///  Given:
///     epj      double    Julian Epoch (e.g. 1996.8)
///
///  Returned:
///     djm0     double    MJD zero-point: always 2400000.5
///     djm      double    Modified Julian Date
///
///  Note:
///
///     The Julian Date is returned in two pieces, in the usual SOFA
///     manner, which is designed to preserve time resolution.  The
///     Julian Date is available as a single number by adding djm0 and
///     djm.
///
///  Reference:
///
///     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
///
///  This revision:  2021 May 11
///
///  SOFA release 2023-10-11
pub fn je_to_jd(epj: f64) -> (f64, f64) {
    (DJM0, DJM00 + (epj - 2000.0) * 365.25)
}

#[derive(Error, Debug, PartialEq)]
pub enum JDDateError {
    #[error("Unacceptable date (Note 1)")]
    UnacceptableDate,
}
///  Julian Date to Gregorian year, month, day, and fraction of a day.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  support function.
///
///  Given:
///     dj1,dj2   double   Julian Date (Notes 1, 2)
///
///  Returned (arguments):
///     iy        int      year
///     im        int      month
///     id        int      day
///     fd        double   fraction of day
///
///  Returned (function value):
///               int      status:
///                           0 = OK
///                          -1 = unacceptable date (Note 1)
///
///  Notes:
///
///  1) The earliest valid date is -68569.5 (-4900 March 1).  The
///     largest value accepted is 1e9.
///
///  2) The Julian Date is apportioned in any convenient way between
///     the arguments dj1 and dj2.  For example, JD=2450123.7 could
///     be expressed in any of these ways, among others:
///
///        dj1             dj2
///
///     2450123.7           0.0       (JD method)
///     2451545.0       -1421.3       (J2000 method)
///     2400000.5       50123.2       (MJD method)
///     2450123.5           0.2       (date & time method)
///
///     Separating integer and fraction uses the "compensated summation"
///     algorithm of Kahan-Neumaier to preserve as much precision as
///     possible irrespective of the jd1+jd2 apportionment.
///
///  3) In early eras the conversion is from the "proleptic Gregorian
///     calendar";  no account is taken of the date(s) of adoption of
///     the Gregorian calendar, nor is the AD/BC numbering convention
///     observed.
///
///  References:
///
///     Explanatory Supplement to the Astronomical Almanac,
///     P. Kenneth Seidelmann (ed), University Science Books (1992),
///     Section 12.92 (p604).
///
///     Klein, A., A Generalized Kahan-Babuska-Summation-Algorithm.
///     Computing, 76, 279-293 (2006), Section 3.
///
///  This revision:  2021 May 11
///
///  SOFA release 2023-10-11
pub fn jd_to_cal(dj1: f64, dj2: f64) -> Result<(i32, i32, i32, f64), JDDateError> {
    // Minimum and maximum allowed JD
    const DJMIN: f64 = -68569.5;
    const DJMAX: f64 = 1e9;

    // Verify date is acceptable.
    let dj = dj1 + dj2;
    if !(DJMIN..=DJMAX).contains(&dj) {
        return Err(JDDateError::UnacceptableDate);
    }

    // Separate day and fraction (where -0.5 <= fraction < 0.5).
    let mut d = dnint(dj1);
    let f1 = dj1 - d;
    let mut jd = d as i64;
    d = dnint(dj2);
    let f2 = dj2 - d;
    jd += d as i64;

    // Compute f1+f2+0.5 using compensated summation (Klein 2006).
    let mut s = 0.5;
    let mut cs = 0.0;
    let v: [f64; 2] = [f1, f2];
    for x in v {
        let t = s + x;
        cs += if s.abs() >= x.abs() {
            (s - t) + x
        } else {
            (x - t) + s
        };
        s = t;
        if s >= 1.0 {
            jd += 1;
            s -= 1.0;
        }
    }
    let mut f = s + cs;
    cs = f - s;

    // Deal with negative f.
    if f < 0.0 {
        // Compensated summation: assume that |s| <= 1.0.
        f = s + 1.0;
        cs += (1.0 - f) + s;
        s = f;
        f = s + cs;
        cs = f - s;
        jd -= 1;
    }

    // Deal with f that is 1.0 or more (when rounded to double).
    if (f - 1.0) >= -f64::EPSILON / 4.0 {
        // Compensated summation: assume that |s| <= 1.0.
        let t = s - 1.0;
        cs += (s - t) - 1.0;
        s = t;
        f = s + cs;
        if -f64::EPSILON / 2.0 < f {
            jd += 1;
            f = f.max(0.0);
        }
    }

    // Express day in Gregorian calendar.
    let mut l = jd + 68569_i64;
    let n = (4_i64 * l) / 146097_i64;
    l -= (146097_i64 * n + 3_i64) / 4_i64;
    let i = (4000_i64 * (l + 1_i64)) / 1461001_i64;
    l -= (1461_i64 * i) / 4_i64 - 31_i64;
    let k = (80_i64 * l) / 2447_i64;
    let id = l - (2447_i64 * k) / 80_i64;
    l = k / 11_i64;
    let im = k + 2_i64 - 12_i64 * l;
    let iy = 100_i64 * (n - 49_i64) + i + l;
    let fd = f;

    Ok((iy as i32, im as i32, id as i32, fd))
}

///  Julian Date to Gregorian Calendar, expressed in a form convenient
///  for formatting messages:  rounded to a specified precision.
///
///  This function is part of the International Astronomical Union's
///  SOFA (Standards of Fundamental Astronomy) software collection.
///
///  Status:  support function.
///
///  Given:
///     ndp       int      number of decimal places of days in fraction
///     dj1,dj2   double   dj1+dj2 = Julian Date (Note 1)
///
///  Returned:
///     iymdf     int[4]   year, month, day, fraction in Gregorian
///                        calendar
///
///  Returned (function value):
///               int      status:
///                          -1 = date out of range
///                           0 = OK
///                          +1 = ndp not 0-9 (interpreted as 0)
///
///  Notes:
///
///  1) The Julian Date is apportioned in any convenient way between
///     the arguments dj1 and dj2.  For example, JD=2450123.7 could
///     be expressed in any of these ways, among others:
///
///    dj1            dj2
///
///  2450123.7           0.0       (JD method)
///  2451545.0       -1421.3       (J2000 method)
///  2400000.5       50123.2       (MJD method)
///  2450123.5           0.2       (date & time method)
///
///  2) In early eras the conversion is from the "Proleptic Gregorian
///     Calendar";  no account is taken of the date(s) of adoption of
///     the Gregorian Calendar, nor is the AD/BC numbering convention
///     observed.
///
///  3) See also the function iauJd2cal.
///
///  4) The number of decimal places ndp should be 4 or less if internal
///     overflows are to be avoided on platforms which use 16-bit
///     integers.
///
///  Called:
///     iauJd2cal    JD to Gregorian calendar
///
///  Reference:
///
///     Explanatory Supplement to the Astronomical Almanac,
///     P. Kenneth Seidelmann (ed), University Science Books (1992),
///     Section 12.92 (p604).
///
///  This revision:  2023 January 16
///
///  SOFA release 2023-10-11
pub fn jd_to_cal_fmt(ndp: i32, dj1: f64, dj2: f64) -> Result<Warned<[i32; 4]>, JDDateError> {
    let mut warning: Option<Warning> = None;

    // Denominator of fraction (e.g. 100 for 2 decimal places).
    let denom = if (0..=9).contains(&ndp) {
        10_i32.pow(ndp as u32) as f64
    } else {
        warning = Some(warning!(1, "ndp not 0-9 (interpreted as 0)"));
        1.0
    };

    // Copy the date, big then small.
    let (mut d1, d2) = if dj1.abs() >= dj2.abs() {
        (dj1, dj2)
    } else {
        (dj2, dj1)
    };

    // Realign to midnight (without rounding error).
    d1 -= 0.5;

    // Separate day and fraction (as precisely as possible).
    let mut d = dnint(d1);
    let f1 = d1 - d;
    let mut djd = d;
    d = dnint(d2);
    let f2 = d2 - d;
    djd += d;
    d = dnint(f1 + f2);
    let mut f = (f1 - d) + f2;
    if f < 0.0 {
        f += 1.0;
        d -= 1.0;
    }
    djd += d;

    // Round the total fraction to the specified number of places.
    let rf = dnint(f * denom) / denom;

    // Re-align to noon.
    djd += 0.5;

    // Convert to Gregorian calendar.
    let (y, m, d, f) = jd_to_cal(djd, rf)?;
    let r_f = dnint(f * denom) as i32;

    Ok(Warned {
        value: [y, m, d, r_f],
        warning,
    })
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    // t_sofa.c t_cal2jd
    #[test]
    fn test_cal_to_jd() {
        let j = cal_to_jd(2003, 6, 1).unwrap();
        assert_eq!(j.value.0, 2400000.5);
        assert_eq!(j.value.1, 52791.0);
    }

    #[test]
    fn test_cal_to_jd_leap_year() {
        // leap year
        let j = cal_to_jd(2400, 2, 29).unwrap();
        assert_eq!(j.ok_code(), 0);
        assert_eq!(j.value.0, 2400000.5);
        assert_eq!(j.value.1, 197700.0);

        // not a leap year
        let j_bad_day = cal_to_jd(2300, 2, 29).unwrap();
        assert!(j_bad_day.value_safe().is_err());
        assert_eq!(j_bad_day.value_safe().err().unwrap().code, -3);
    }

    #[test]
    fn test_cal_to_jd_invalid_month() {
        let j_bad_month = cal_to_jd(2003, 13, 1);
        assert_eq!(j_bad_month.err(), Some(DateError::InvalidMonth));
    }

    #[test]
    fn test_cal_to_jd_parity() {
        use rsofa::iauCal2jd;
        let j = cal_to_jd(2003, 6, 1).unwrap();
        let mut djm0 = 0.0;
        let mut djm = 0.0;
        let status_iau = unsafe { iauCal2jd(2003, 6, 1, &mut djm0, &mut djm) };
        assert_eq!(status_iau, j.ok_code());
        assert_eq!(djm0, j.value.0);
        assert_eq!(djm, j.value.1);
    }
    // t_sofa.c t_epb;
    #[test]
    fn test_jd_to_besselian_e() {
        let epb = jd_to_besselian_e(2415019.8135, 30103.18648);
        assert_approx_eq!(epb, 1982.418424159278580, 1e-12);
    }

    #[test]
    fn test_jd_to_besselian_e_parity() {
        use rsofa::iauEpb;
        let epb = jd_to_besselian_e(2415019.8135, 30103.18648);
        let epb_iau = unsafe { iauEpb(2415019.8135, 30103.18648) };
        assert_eq!(epb, epb_iau);
    }

    // t_sofa.c t_epb2jd
    #[test]
    fn test_besselian_e_to_jd() {
        let epb = 1957.3;
        let (djm0, djm) = besselian_e_to_jd(epb);
        assert_approx_eq!(djm0, 2400000.5, 1e-9);
        assert_approx_eq!(djm, 35948.1915101513, 1e-9);
    }
    #[test]
    fn test_besselian_e_to_jd_parity() {
        use rsofa::iauEpb2jd;
        let epb = 1957.3;
        let (djm0, djm) = besselian_e_to_jd(epb);

        let mut djm0_iau = 0.0;
        let mut djm_iau = 0.0;
        unsafe { iauEpb2jd(epb, &mut djm0_iau, &mut djm_iau) };
        assert_eq!(djm0, djm0_iau);
        assert_eq!(djm, djm_iau);
    }

    // t_sofa.c t_epj
    #[test]
    fn test_jd_to_je() {
        let epj = jd_to_je(2451545.0, -7392.5);
        assert_approx_eq!(epj, 1979.760438056125941, 1e-12);
    }

    #[test]
    fn test_jd_to_je_parity() {
        use rsofa::iauEpj;
        let epj = jd_to_je(2451545.0, -7392.5);
        let epj_iau = unsafe { iauEpj(2451545.0, -7392.5) };
        assert_eq!(epj, epj_iau);
    }

    // t_sofa.c t_epj2jd
    #[test]
    fn test_je_to_jd() {
        let epj = 1996.8;
        let (djm0, djm) = je_to_jd(epj);
        assert_approx_eq!(djm0, 2400000.5, 1e-9);
        assert_approx_eq!(djm, 50375.7, 1e-9);
    }
    #[test]
    fn test_je_to_jd_parity() {
        use rsofa::iauEpj2jd;
        let epj = 1996.8;
        let (djm0, djm) = je_to_jd(epj);
        let mut djm0_iau = 0.0;
        let mut djm_iau = 0.0;
        unsafe { iauEpj2jd(epj, &mut djm0_iau, &mut djm_iau) };
        assert_eq!(djm0, djm0_iau);
        assert_eq!(djm, djm_iau);
    }

    // t_sofa.c t_jd2cal
    #[test]
    fn test_jd_to_cal() {
        let dj1 = 2400000.5;
        let dj2 = 50123.9999;

        let (iy, im, id, fd) = jd_to_cal(dj1, dj2).unwrap();

        assert_eq!(iy, 1996);
        assert_eq!(im, 2);
        assert_eq!(id, 10);
        assert_approx_eq!(fd, 0.9999, 1e-7);
    }

    #[test]
    fn test_jd_to_cal_parity() {
        use rsofa::iauJd2cal;
        let dj1 = 2400000.5;
        let dj2 = 50123.9999;

        let (iy, im, id, fd) = jd_to_cal(dj1, dj2).unwrap();

        let mut iy_iau = 0;
        let mut im_iau = 0;
        let mut id_iau = 0;
        let mut fd_iau = 0.0;
        unsafe { iauJd2cal(dj1, dj2, &mut iy_iau, &mut im_iau, &mut id_iau, &mut fd_iau) };
        assert_eq!(iy, iy_iau);
        assert_eq!(im, im_iau);
        assert_eq!(id, id_iau);
        assert_eq!(fd, fd_iau);
    }

    // t_sofa.c t_jdcalf
    #[test]
    fn test_jd_to_cal_fmt() {
        let dj1 = 2400000.5;
        let dj2 = 50123.9999;

        let iymdf = jd_to_cal_fmt(4, dj1, dj2).unwrap().value;

        assert_eq!(iymdf[0], 1996);
        assert_eq!(iymdf[1], 2);
        assert_eq!(iymdf[2], 10);
        assert_eq!(iymdf[3], 9999);
    }
    #[test]
    fn test_jd_to_cal_fmt_parity() {
        use rsofa::iauJdcalf;
        let dj1 = 2400000.5;
        let dj2 = 50123.9999;

        let iymdf = jd_to_cal_fmt(4, dj1, dj2).unwrap().value;

        let mut iymdf_iau = [0; 4];
        unsafe {
            iauJdcalf(4, dj1, dj2, iymdf_iau.as_mut_ptr());
        }

        assert_eq!(iymdf, iymdf_iau);
    }
}
