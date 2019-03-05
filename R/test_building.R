# Common properties
constants <- dplyr::tibble(
  lambda_p = list(c(0.00, 0.04, 0.11, 0.25, 0.44)), # Plan density
  lambda_w = list(c(1.00, 0.90, 0.70, 0.50, 0.30)), # Sheltering factor wind
  lambda_sol=list(c(5.00, 3.00, 2.00, 1.30, 0.90)), # Shading factor solar
  U_alpha =   list(c(0.10, 0.13, 0.18, 0.25, 0.33)), # c(0.33, 0.22, 0.14, 0.10) in ASHRAE fund 2017, ch 24 tbl 1
  U_delta=   list(c(220, 257, 320, 393, 460)),      # c(460, 370, 270, 220) in ASHRAE fund 2017, ch 24 tbl 1
  # Specific heat capacity of opaque elements table B.14 ISO 52016
  # 1. Very light, 2. Light, 3. Medium, 4. Heavy, 5. Very heavy
  k_op      =list(c(50000, 75000, 110000, 175000, 250000)/3600), # Wh/(m2*K)
  p_a = 1.204,                  # density of air kg/m3 at 20c
  k_a = 1006,                   # heat capacity of air [J/(kg*K)] at 20c
  kp_a = k_a*p_a/1000,          # heat capacity of air per volume [J/(l*K)] or [Ws/(l*K)]
  Fsky_hor = 1, Fsky_ver=0.5, # Sky view factor
  abs_sol = 0.5,                # solar absorption fraction
  Rf_rf = 1.11, Rf_ew = 1.52,   # Surface Roughness Multipliers (EnergyPlus ref)
  h_re = 4.14, h_ri = 5.13,     # surface radiative heat transfer coeffs
  U10_idx_N = 10,               # number of categories to split wind speed into
  Kv_sol = 115                  # Global luminous efficacy ISO52010-1 eq 43 [lm/W]
)


# basic building and thermal properties
p <- constants %>% dplyr::mutate(
  # Basic geometrics and coeffs
  N_fl = 1,                     # number of floor
  H_fl = 2.6,                   # average internal floor height
  H_fl2 = H_fl + 0.25,          # average floor level height
  P = (20 + 10)*2,              # perimeter
  A_fl = (20 * 10)*N_fl,        # floor area
  H_b = N_fl*H_fl2,             # Building height
  class = 3,                    # Shelter class
  lat = 58.575,                 # latitude, degrees
  lng = 16.15,                  # longitude, degrees
  az = list(c(180, 90, 0, -90)),# external wall surface azimuths, # 180: north, 90: east, 0: south, -90: west
  Faz=list(c(10, 20, 10, 20)/P),# weighting factor for each az

  # weighting ratios of the bulding elements
  r_gl = 0.15,
  r_rf = (1/N_fl),
  r_gf = r_rf,
  r_ew = P*H_fl*N_fl/A_fl - r_gl, # incl. window frames
  r_im = (2 - 2/N_fl) + 1.5,
  r_si = list(c(r_rf, r_ew, r_gl, r_im, r_gf)),
  r_tot= sum(r_si[[1]]),

  # surface area fractions of bulding elements
  f_rf = r_rf/r_tot,
  f_ew = r_ew/r_tot,
  f_gl = r_gl/r_tot,
  f_im = r_im/r_tot,
  f_gf = r_gf/r_tot,

  # convective fractions. table B.11 ISO52016
  f_c_hyd = 0.4,
  f_c_int = 0.4,
  f_c_sol = 0.1,                # 0.1 in ISO52016, 0.4 in 13790

  # solar properties window and blinds
  g_gl = 0.76*0.9,              # g-value glazing incl. effect of window frames
  g_bl_max = 0.53,              # g-value blind fully drawn
  u_bl_min = 0.3,               # min control signal blinds
  u_bl_max = 0.7,               # max control signal blinds

  # thermal properties
  U_rf = 0.199,
  U_ew = 0.716,
  U_gl = 2.9,
  U_gf = 0.23,
  H_tb = 0.00,                  # thermal bridges
  H_gr_vi = 0.074,
  k_op_cl = 4,                  # areal heat capacity class of opaque element [Wh/(m2*K)]
  # internal heat capacity furniture + air volyme. ISO52016 table B.17: 10000 J/(m2*K)
  C_int = 1.87 + H_fl*1006/3600,# [Wh/(m2*K)]

  # ventilation. Constant flow
  Q_ve = 0.35,                  # air flow rate, [l/(s*m2)]
  eta_ve = 0.6*0,                 # temperature transfer efficiency

  # Infiltration
  H_mod = ifelse(N_fl > 2, N_fl/3*H_fl/2, N_fl*H_fl/2), # modified height for infiltration calc
  C_w = 0.22,                   # wind effect factor
  C_s = 0.25,                   # stack effect factor
  C_inf = 0.08,                 # infiltration or flow coefficient
  n_inf = 0.67,                 # flow exponent

  # Heating system
  T_int_set = 21,                 # internal air setpoint
  # radiator supply temperature lookup table
  xy = list(matrix(c(c(-20, -15, -10, 0, 12, 18), c(60, 55.6, 51, 43, 28, 18)), ncol=2)),
  H_hyd = 0.7,                    # radiator constant
  n_hyd = 1.28,                 # radiator exponent
  T_hyd_d = 60,                 # design supply temperature
  dT_hyd_d = 20,                # design temperature drop
  # correlation coefficients used for radiator calc
  a_hyd = (n_hyd - dT_hyd_d/200),
  b_hyd = dT_hyd_d/(T_hyd_d - 20)^a_hyd,
  T_trv_pb = 2,                 # proportional band of the TRV:s
  T_trv_set = T_int_set + T_trv_pb/2, # max temp set-point of the TRV:s

  # lighting system
  Kv_L = 25,                    # Luminous efficacy lighting [lm/W]
  Ev_min = 50,                  # Minumum used illuminance level [lx] [lm/m2]
  Ev_max = 150                  # Installed lighting illuminance [lx]
)

