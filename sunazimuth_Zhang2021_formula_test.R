
# testing formula described in https://www.sciencedirect.com/science/article/pii/S0960148121004031
# Zhang et al., 2021
# for north-clockwise convention
# for cube "C0466670113R.cub

# ^IMAGE_HEADER = ("C0466670113R.IMG",1)

# /* Camera and Lighting Geometry                        */
# /* Resolution of HORIZONTAL_PIXEL_SCALE,               */
# /* VERTICAL_PIXEL_SCALE, and SLANT_DISTANCE            */
# /* is calculated from the light source values in       */
# /* INTERCEPT_POINT_LATITUDE, INTERCEPT_POINT_LONGITUDE,*/
# /* INTERCEPT_POINT_LINE and INTERCEPT_POINT_LINE_SAMPLE keywords  */
# /* If the target is a Ring keyword RING_RADIUS         */
# /* is substituted for INTERCEPT_POINT_LATITUDE         */
# /* If the TARGET_NAME = J RINGS, viewing geometry was */
# /* calculated using Jupiter as the target. */
# TWIST_ANGLE =   33.053
# CONE_ANGLE =  148.371
# RIGHT_ASCENSION =  334.264
# DECLINATION =   21.560
# NORTH_AZIMUTH =  173.412
# SMEAR_AZIMUTH = "UNK"
# SMEAR_MAGNITUDE = "UNK"
# HORIZONTAL_PIXEL_SCALE = 5.924380e+01
# VERTICAL_PIXEL_SCALE = 4.346910e+01
# SLANT_DISTANCE = 4.183800e+03
# SOLAR_DISTANCE = 7.424390e+08
# SUB_SOLAR_LATITUDE =    2.028
# SUB_SOLAR_LONGITUDE =  150.387
# SUB_SOLAR_AZIMUTH =  189.350
# INCIDENCE_ANGLE =   77.018
# EMISSION_ANGLE =   44.806
# PHASE_ANGLE =   32.271
# LOCAL_HOUR_ANGLE =  133.369
# INTERCEPT_POINT_LATITUDE =  -68.110
# INTERCEPT_POINT_LONGITUDE =  196.672
# INTERCEPT_POINT_LINE = 400.0
# INTERCEPT_POINT_LINE_SAMPLE = 400.0

# subsolar point
phi_s = pi/180* 2.028
lambda_s = pi/180* (360-150.387)

# observer point (on the surface)
phi_0 = pi/180* (-68.11) # INTERCEPT_POINT_LATITUDE
lambda_0 = pi/180* (360-196.672)


S_x <- cos(phi_s)*sin(lambda_s - lambda_0)
S_y <- cos(phi_0)*sin(phi_s) - sin(phi_0)*cos(phi_s)*cos(lambda_s - lambda_0)
gamma_s <- atan2(S_x, S_y)
gamma_s_deg <- gamma_s * 180/pi

gamma_arcgis <- 47.46

relative_difference <- (gamma_s_deg - gamma_arcgis)/gamma_s_deg
