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
