SoFar

# Implementation Status

## Main Routines

### Calendars

|  | sofar | iau | description |
| --- | --- | --- | --- |
| ✅ | cal_to_jd |                 CAL2JD |    Gregorian calendar to Julian Day number |
| ✅ | jd_to_besselian_e |         EPB |       Julian Date to Besselian Epoch |
| ✅ | besselian_e_to_jd |         EPB2JD |    Besselian Epoch to Julian Date |
| ✅ | jd_to_je |                  EPJ |       Julian Date to Julian Epoch |
| ✅ | je_to_jd |                  EPJ2JD |    Julian Epoch to Julian Date |
| ✅ | jd_to_cal |                 JD2CAL |    Julian Date to Gregorian year, month, day, fraction |
| ✅ | jd_to_cal_fmt |             JDCALF |    Julian Date to Gregorian date for formatted output |


### Astrometry
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ❌ | stellar_aberration|             AB |        apply stellar aberration |
| ❌ | prep_i2g_geocentric|            APCG |      prepare for ICRS <-> GCRS, geocentric, special |
| ❌ | prep_i2g_geocentric_13|         APCG13 |    prepare for ICRS <-> GCRS, geocentric |
| ❌ | prep_i2c_terrestrial|           APCI |      prepare for ICRS <-> CIRS, terrestrial, special |
| ❌ | prep_i2c_terrestrial_13|        APCI13 |    prepare for ICRS <-> CIRS, terrestrial |
| ❌ | prep_i2obs_terrestrial|         APCO |      prepare for ICRS <-> observed, terrestrial, special |
| ❌ | prep_i2obs_terrestrial_13|      APCO13 |    prepare for ICRS <-> observed, terrestrial |
| ❌ | prep_i2c_space|                 APCS |      prepare for ICRS <-> CIRS, space, special |
| ❌ | prep_i2c_space_13|              APCS13 |    prepare for ICRS <-> CIRS, space |
| ❌ | earth_rotation|                 APER |      insert ERA into context |
| ❌ | earth_rotation_13|              APER13 |    update context for Earth rotation |
| ❌ | prep_c2o_terrestrial|           APIO |      prepare for CIRS <-> observed, terrestrial, special |
| ❌ | prep_c2o_terrestrial_13|        APIO13 |    prepare for CIRS <-> observed, terrestrial |
| ❌ | catalog_to_astrometric|         ATCC13 |    catalog -> astrometric |
| ❌ | catalog_to_astrometric_quick|   ATCCQ |     quick catalog -> astrometric |
| ❌ | catalog_to_cirs|                ATCI13 |    catalog -> CIRS |
| ❌ | icrs_to_cirs_quick|             ATCIQ |     quick ICRS -> CIRS |
| ❌ | icrs_to_cirs_quick_def|         ATCIQN |    quick ICRS -> CIRS, multiple deflections |
| ❌ | astro_icrs_to_cirs_quick|       ATCIQZ |    quick astrometric ICRS -> CIRS |
| ❌ | icrs_to_observed|               ATCO13 |    ICRS -> observed |
| ❌ | cirs_to_icrs|                   ATIC13 |    CIRS -> ICRS |
| ❌ | cirs_to_icrs_quick|             ATICQ |     quick CIRS -> ICRS |
| ❌ | cirs_to_icrs_quick_def|         ATICQN |    quick CIRS -> ICRS, multiple deflections |
| ❌ | cirs_to_observed|               ATIO13 |    CIRS -> observed |
| ❌ | cirs_to_observed_quick|         ATIOQ |     quick CIRS -> observed |
| ❌ | observed_to_astro_icrs|         ATOC13 |    observed -> astrometric ICRS |
| ❌ | observed_to_cirs|               ATOI13 |    observed -> CIRS |
| ❌ | observed_to_cirs_quick|         ATOIQ |     quick observed -> CIRS |
| ❌ | light_deflection|               LD |        light deflection by a single solar-system body |
| ❌ | light_deflection_n|             LDN |       light deflection by multiple solar-system bodies |
| ❌ | light_delfection_sun|           LDSUN |     light deflection by the Sun |
| ❌ | appy_pm_px|                     PMPX |      apply proper motion and parallax |
| ❌ | apply_pm_safe|                  PMSAFE |    apply proper motion, with zero-parallax precautions |
| ❌ | observatory_pv|                 PVTOB |     observatory position and velocity |
| ❌ | pv_to_star_catalog|             PVSTAR |    space motion pv-vector to star catalog data |
| ❌ | refraction_constants|           REFCO |     refraction constants |
| ❌ | star_pm|                        STARPM |    apply proper motion |
| ❌ | star_pv|                        STARPV |    star catalog data to space motion pv-vector |

### Time
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ❌ |format_jd|         D2DTF|     format 2-part JD for output|
| ❌ |delta_at|          DAT|       Delta(AT) (=TAI-UTC) for a given UTC date|
| ❌ |tdb_tt_diff|       DTDB|      TDB-TT|
| ❌ |datetime_to_jd|    DTF2D|     encode time and date fields into 2-part JD|
| ❌ |tai_to_tt|         TAITT|     TAI to TT|
| ❌ |tai_to_ut1|        TAIUT1|    TAI to UT1|
| ❌ |tai_to_utc|        TAIUTC|    TAI to UTC|
| ❌ |tcb_to_tdb|        TCBTDB|    TCB to TDB|
| ❌ |tcg_to_tt|         TCGTT|     TCG to TT|
| ❌ |tdb_to_tcb|        TDBTCB|    TDB to TCB|
| ❌ |tdb_to_tt|         TDBTT|     TDB to TT|
| ❌ |tt_to_tai|         TTTAI|     TT to TAI|
| ❌ |tt_to_tcg|         TTTCG|     TT to TCG|
| ❌ |tt_to_tdb|         TTTDB|     TT to TDB|
| ❌ |tt_to_ut1|         TTUT1|     TT to UT1|
| ❌ |ut1_to_tai|        UT1TAI|    UT1 to TAI|
| ❌ |ut1_to_tt|         UT1TT|     UT1 to TT|
| ❌ |ut1_to_utc|        UT1UTC|    UT1 to UTC|
| ❌ |utc_to_tai|        UTCTAI|    UTC to TAI|
| ❌ |utc_to_ut1|        UTCUT1|    UTC to UT1|

### Earth Rotation
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ❌ |ee_2000|     EE00|      equation of the equinoxes, IAU 2000|
| ❌ |ee_2000a|    EE00A|     equation of the equinoxes, IAU 2000A|
| ❌ |ee_2000b|    EE00B|     equation of the equinoxes, IAU 2000B|
| ❌ |ee_2006a|    EE06A|     equation of the equinoxes, IAU 2006/2000A|
| ❌ |ee_ct_2000|  EECT00|    equation of the equinoxes complementary terms, IAU 2000|
| ❌ |ee_1994|     EQEQ94|    equation of the equinoxes, IAU 1994|
| ❌ |era_2000|    ERA00|     Earth rotation angle, IAU 2000|
| ❌ |gmst_2000|   GMST00|    Greenwich mean sidereal time, IAU 2000|
| ❌ |gmst_2006|   GMST06|    Greenwich mean sidereal time, IAU 2006|
| ❌ |gmst_1982|   GMST82|    Greenwich mean sidereal time, IAU 1982|
| ❌ |gst_2000a|   GST00A|    Greenwich apparent sidereal time, IAU 2000A|
| ❌ |gst_2000b|   GST00B|    Greenwich apparent sidereal time, IAU 2000B|
| ❌ |gst_2006|    GST06|     Greenwich apparent ST, IAU 2006, given NPB matrix|
| ❌ |gst_2006a|   GST06A|    Greenwich apparent sidereal time, IAU 2006/2000A|
| ❌ |gst_1994|    GST94|     Greenwich apparent sidereal time, IAU 1994|

### Ephemerides
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ❌ |earth_pv|   EPV00|     Earth position and velocity|
| ❌ |moon_pv|    MOON98|    Moon position and velocity|
| ❌ |planet_pv|  PLAN94|    major-planet position and velocity|

## Precession, nutation, polar motion
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ❌ |bias_2000|                 BI00|      frame bias components, IAU 2000|
| ❌ |bias_precession_2000|      BP00|      frame bias and precession matrices, IAU 2000|
| ❌ |bias_precession_2006|      BP06|      frame bias and precession matrices, IAU 2006|
| ❌ |npb_to_xy|                 BPN2XY|    extract CIP X,Y coordinates from NPB matrix|
| ❌ |c2i_2000a|                 C2I00A|    celestial-to-intermediate matrix, IAU 2000A|
| ❌ |c2i_2000b|                 C2I00B|    celestial-to-intermediate matrix, IAU 2000B|
| ❌ |c2i_2006|                  C2I06A|    celestial-to-intermediate matrix, IAU 2006/2000A|
| ❌ |c2i_npb|                   C2IBPN|    celestial-to-intermediate matrix, given NPB matrix, IAU 2000|
| ❌ |c2i_xy|                    C2IXY|     celestial-to-intermediate matrix, given X,Y, IAU 2000|
| ❌ |c2i_xys|                   C2IXYS|    celestial-to-intermediate matrix, given X,Y and s|
| ❌ |c2t_2000a|                 C2T00A|    celestial-to-terrestrial matrix, IAU 2000A|
| ❌ |c2t_2000b|                 C2T00B|    celestial-to-terrestrial matrix, IAU 2000B|
| ❌ |c2t_2006|                  C2T06A|    celestial-to-terrestrial matrix, IAU 2006/2000A|
| ❌ |c2t_cio|                   C2TCIO|    form CIO-based celestial-to-terrestrial matrix|
| ❌ |c2t_eqx|                   C2TEQX|    form equinox-based celestial-to-terrestrial matrix|
| ❌ |c2t_nutation|              C2TPE|     celestial-to-terrestrial matrix given nutation, IAU 2000|
| ❌ |xy_to_c2t|                 C2TXY|     celestial-to-terrestrial matrix given CIP, IAU 2000|
| ❌ |eo_2006|                   EO06A|     equation of the origins, IAU 2006/2000A|
| ❌ |eo_npb_s|                  EORS|      equation of the origins, given NPB matrix and s|
| ❌ |fw_to_rmatrix|             FW2M|      Fukushima-Williams angles to r-matrix|
| ❌ |fw_to_xy|                  FW2XY|     Fukushima-Williams angles to X,Y|
| ❌ |long_term_precession|      LTP|       long-term precession matrix|
| ❌ |long_term_precession_b|    LTPB|      long-term precession matrix, including ICRS frame bias|
| ❌ |long_term_precession_ecl|  LTPECL|    long-term precession of the ecliptic|
| ❌ |long_term_precession_eq|   LTPEQU|    long-term precession of the equator|
| ❌ |nutation_matrix_2000a|     NUM00A|    nutation matrix, IAU 2000A|
| ❌ |nutation_matrix_2000b|     NUM00B|    nutation matrix, IAU 2000B|
| ❌ |nutation_matrix_2006|      NUM06A|    nutation matrix, IAU 2006/2000A|
| ❌ |nutation_matrix|           NUMAT|     form nutation matrix|
| ❌ |nutation_2000a|            NUT00A|    nutation, IAU 2000A|
| ❌ |nutation_2000b|            NUT00B|    nutation, IAU 2000B|
| ❌ |nutation_2006|             NUT06A|    nutation, IAU 2006/2000A|
| ❌ |nutation_1980|             NUT80|     nutation, IAU 1980|
| ❌ |nutation_matrix_1980|      NUTM80|    nutation matrix, IAU 1980|
| ❌ |obliquity_2006|            OBL06|     mean obliquity, IAU 2006|
| ❌ |obliquity_1980|            OBL80|     mean obliquity, IAU 1980|
| ❌ |precession_bias_2006|      PB06|      zeta,z,theta precession angles, IAU 2006, including bias|
| ❌ |precession_fw_2006|        PFW06|     bias-precession Fukushima-Williams angles, IAU 2006|
| ❌ |precession_matrix_2000|    PMAT00|    precession matrix (including frame bias), IAU 2000|
| ❌ |precession_matrix_2006|    PMAT06|    PB matrix, IAU 2006|
| ❌ |precession_matrix_1976|    PMAT76|    precession matrix, IAU 1976|
| ❌ |pn_2000|                   PN00|      bias/precession/nutation results, IAU 2000|
| ❌ |pn_2000a|                  PN00A|     bias/precession/nutation, IAU 2000A|
| ❌ |pn_2000b|                  PN00B|     bias/precession/nutation, IAU 2000B|
| ❌ |pn_2006|                   PN06|      bias/precession/nutation results, IAU 2006|
| ❌ |pn_2006a|                  PN06A|     bias/precession/nutation results, IAU 2006/2000A|
| ❌ |pnm_2000a|                 PNM00A|    classical NPB matrix, IAU 2000A|
| ❌ |pnm_2000b|                 PNM00B|    classical NPB matrix, IAU 2000B|
| ❌ |pnm_2006a|                 PNM06A|    classical NPB matrix, IAU 2006/2000A|
| ❌ |pnm_1980|                  PNM80|     precession/nutation matrix, IAU 1976/1980|
| ❌ |p_2006_equinox|            P06E|      precession angles, IAU 2006, equinox based|
| ❌ |polar_motion_matrix|       POM00|     polar motion matrix|
| ❌ |pr_2000|                   PR00|      IAU 2000 precession adjustments|
| ❌ |prec_1976|                 PREC76|    accumulated precession angles, IAU 1976|
| ❌ |s_2000|                    S00|       the CIO locator s, given X,Y, IAU 2000A|
| ❌ |s_2000a|                   S00A|      the CIO locator s, IAU 2000A|
| ❌ |s_2000b|                   S00B|      the CIO locator s, IAU 2000B|
| ❌ |s_2006|                    S06|       the CIO locator s, given X,Y, IAU 2006|
| ❌ |s_2006a|                   S06A|      the CIO locator s, IAU 2006/2000A|
| ❌ |sp_2000|                   SP00|      the TIO locator s', IERS 2003|
| ❌ |xy_2006a|                  XY06|      CIP, IAU 2006/2000A, from series|
| ❌ |xys_2000a|                 XYS00A|    CIP and s, IAU 2000A|
| ❌ |xys_2000b|                 XYS00B|    CIP and s, IAU 2000B|
| ❌ |xys_2006a|                 XYS06A|    CIP and s, IAU 2006/2000A|

### Fundamental arguments for nutation, etc.
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ❌ |fa_d_2003|     FAD03|     mean elongation of the Moon from the Sun|
| ❌ |fa_e_2003|     FAE03|     mean longitude of Earth|
| ❌ |fa_f_2003|     FAF03|     mean argument of the latitude of the Moon|
| ❌ |fa_ju_2003|    FAJU03|    mean longitude of Jupiter|
| ❌ |fa_l_2003|     FAL03|     mean anomaly of the Moon|
| ❌ |fa_p_2003|     FALP03|    mean anomaly of the Sun|
| ❌ |fa_ma_2003|    FAMA03|    mean longitude of Mars|
| ❌ |fa_me_2003|    FAME03|    mean longitude of Mercury|
| ❌ |fa_ne_2003|    FANE03|    mean longitude of Neptune|
| ❌ |fa_om_2003|    FAOM03|    mean longitude of the Moon's ascending node|
| ❌ |fa_pa_2003|    FAPA03|    general accumulated precession in longitude|
| ❌ |fa_sa_2003|    FASA03|    mean longitude of Saturn|
| ❌ |fa_ur_2003|    FAUR03|    mean longitude of Uranus|
| ❌ |fa_ve_2003|    FAVE03|    mean longitude of Venus|

### Star Catalog Conversions
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ❌ |fk5_to_hip|              FK52H|     transform FK5 star data into the Hipparcos system|
| ❌ |fk5_to_hipp_rot|         FK5HIP|    FK5 to Hipparcos rotation and spin|
| ❌ |fk5_to_hipp_zero|        FK5HZ|     FK5 to Hipparcos assuming zero Hipparcos proper motion|
| ❌ |hipp_to_fk5|             H2FK5|     transform Hipparcos star data into the FK5 system|
| ❌ |hipp_to_fk5_zero|        HFK5Z|     Hipparcos to FK5 assuming zero Hipparcos proper motion|
| ❌ |fk4_to_fk5|              FK425|     transform FK4 star data into FK5|
| ❌ |fk4_to_fk5_zero|         FK45Z|     FK4 to FK5 assuming zero FK5 proper motion|
| ❌ |fk5_to_fk4|              FK524|     transform FK5 star data into FK4|
| ❌ |fk5_to_fk4_zero|         FK54Z|     FK5 to FK4 assuming zero FK5 proper motion|


## Vector/Matrix Library

### Operations on Angles
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ✅ |wrap::angle_normalize_positive|   ANP|       normalize radians to range 0 to 2pi|
| ✅ |wrap::angle_normalize_pm|         ANPM|      normalize radians to range -pi to +pi|
| ✅ |to_sexagesimal::radians_to_hms|  A2TF|      decompose radians into hours, minutes, seconds|
| ✅ |to_sexagesimal::radians_to_dms|  A2AF|      decompose radians into degrees, arcminutes, arcseconds|
| ✅ |to_sexagesimal::days_to_hms|     D2TF|      decompose days into hours, minutes, seconds|
| ✅ |from_sexagesimal::dms_to_radians|    AF2A|      degrees, arcminutes, arcseconds to radians|
| ✅ |from_sexagesimal::hms_to_radians|    TF2A|      hours, minutes, seconds to radians|
| ✅ |from_sexagesimal::hms_to_days|       TF2D|      hours, minutes, seconds to days|

### Operations Involving P-Vectors and R-Matrices
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ✅ | initialize::zero_pvector|        ZP|        zero p-vector|
| ✅ | initialize::zero_rmatrix|        ZR|        initialize r-matrix to null|
| ✅ | initialize::identity_rmatrix|    IR|        initialize r-matrix to identity|
| ✅ | copy::cp|  CP|        copy p-vector|
| ✅ | copy::cr|  CR|        copy r-matrix|
| ✅ |rotations::rotate_rmatrix_about_x|    RX|        rotate r-matrix about x|
| ✅ |rotations::rotate_rmatrix_about_y|    RY|        rotate r-matrix about y|
| ✅ |rotations::rotate_rmatrix_about_z|    RZ|        rotate r-matrix about z|
| ✅ |sphere_cart_conv::spherical_to_unit_vector|  S2C|       spherical to unit vector|
| ✅ |sphere_cart_conv::unit_vector_to_spherical|  C2S|       unit vector to spherical|
| ✅ |sphere_cart_conv::spherical_to_pvector|      S2P|       spherical to p-vector|
| ✅ |sphere_cart_conv::p_vector_to_spherical|     P2S|       p-vector to spherical|
| ✅ |vec_ops::pvector_plus_pvector|          PPP|       p-vector plus p-vector|
| ✅ |vec_ops::pvector_minus_pvector|         PMP|       p-vector minus p-vector|
| ✅ |vec_ops::pvector_plus_scaled_pvector|   PPSP|      p-vector plus scaled p-vector|
| ✅ |vec_ops::pvector_dot_product|           PDP|       inner (=scalar=dot) product of two p-vectors|
| ✅ |vec_ops::pvector_cross_product|         PXP|       outer (=vector=cross) product of two p-vectors|
| ✅ |vec_ops::pvector_modulus|               PM|        modulus of p-vector|
| ✅ |vec_ops::pvector_normalize|             PN|        normalize p-vector returning modulus|
| ✅ |vec_ops::pvector_multiply_scalar|       SXP|       multiply p-vector by scalar|
| ✅ |matrix_ops::rmatrix_multiply|      RXR|       r-matrix multiply|
| ✅ |matrix_ops::transpose_rmatrix|     TR|        transpose r-matrix|
| ✅ |matrix_vec_products::rmatrix_pvector_product|           RXP|       product of r-matrix and p-vector|
| ✅ |matrix_vec_products::transpose_rmatrix_pvector_product| TRXP|      product of transpose of r-matrix and p-vector|
| ✅ |sep_position_angle::angular_separation_pvector|    SEPP|      angular separation from p-vectors|
| ✅ |sep_position_angle::angular_separation_spherical|  SEPS|      angular separation from spherical coordinates|
| ✅ |sep_position_angle::position_angle_from_pvector|   PAP|       position-angle from p-vectors|
| ✅ |sep_position_angle::position_angle_from_spherical| PAS|       position-angle from spherical coordinates|
| ✅ |rotation_vectors::rvector_to_rmatrix|    RV2M|      r-vector to r-matrix|
| ✅ |rotation_vectors::rmatrix_to_rvector|    RM2V|      r-matrix to r-vector|

### Operations Involving PV-Vectors
|  | sofar | iau | description |
| --- | --- | --- | --- |
| ✅ |initialize::zero_pvvector| ZPV|       zero pv-vector|
| ✅ |copy_extend_extract::copy_pvvector|                  CPV|       copy pv-vector|
| ✅ |copy_extend_extract::append_zvelocity_pvvector|      P2PV|      append zero velocity to p-vector|
| ✅ |copy_extend_extract::discard_velocity_pvvector|      PV2P|      discard velocity component of pv-vector|
| ✅ |sphere_cart_conv::spherical_to_pvvector|     S2PV|      spherical to pv-vector|
| ✅ |sphere_cart_conv::pvvector_to_spherical|     PV2S|      pv-vector to spherical|
| ✅ |pvector_ops::pvvector_plus_pvvector|            PVPPV|     pv-vector plus pv-vector|
| ✅ |pvector_ops::pvvector_minus_pvvector|           PVMPV|     pv-vector minus pv-vector|
| ✅ |pvector_ops::pvvector_dot_pvvector|             PVDPV|     inner (=scalar=dot) product of two pv-vectors|
| ✅ |pvector_ops::pvvector_cross_pvvector|           PVXPV|     outer (=vector=cross) product of two pv-vectors|
| ✅ |pvector_ops::pvvector_modulus|                  PVM|       modulus of pv-vector|
| ✅ |pvector_ops::pvvector_multiply_scalar|          SXPV|      multiply pv-vector by scalar|
| ✅ |pvector_ops::pvvector_multiply_two_scalar|      S2XPV|     multiply pv-vector by two scalars|
| ✅ |pvector_ops::pvvector_update|                   PVU|       update pv-vector|
| ✅ |pvector_ops::pvvector_update_discard_velocity|  PVUP|      update pv-vector discarding velocity|
| ✅ |matrix_vector_products::rmatrix_multiply_pvvector|             RXPV|      product of r-matrix and pv-vector|
| ✅ |matrix_vector_products::rmatrix_multiply_pvvector_transpose|   TRXPV|     product of transpose of r-matrix and pv-vector|
