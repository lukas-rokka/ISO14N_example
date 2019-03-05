# Copyright (C) 2019 Lukas Lundström


#' @description Calculate total, surface area weighte, solar irradiance on vertical surfaces. eq. 29 in [1].
#' @export
calc_I_ver_sh <- function(.df, p) {
  tmp <- .df %>% solarCalcISO52010::tidyISO52010(
    p$lat, p$lng, 0, NULL,
    surfaceAzimuths = p$az[[1]],
    surfaceTilts = rep(90, length(p$az[[1]]))) %>%
    select("I_tot_dif_s1", starts_with("I_tot_dir"))
  # for each surface: multiply I_tot_dir with shading factor F_sh_ver and add I_tot_dif,
  # row wise weighted sum with factor Faz
  .df %>% mutate(I_tot_ver_sh = (as.matrix(tmp[, 2:(ncol(tmp))] * .df$F_sh_ver + tmp[, 1]) %*% p$Faz[[1]])[ ,1])
}


#' @description Calculate effect of window blinds. Section 3.5 in [1]
#' @export
calc_blinds <- function(.df, p) {
  l = 24*1/(diff(as.numeric(.df$timestamp[1:2]))/3600)         # roll back 24 hours. 
  .df %>%  mutate(
  I_wma = RcppRoll::roll_meanr(c(I_tot_ver_sh[(l-1):1], I_tot_ver_sh), l, 1:l, fill=NULL), # weighted moving average
  u_bl = pmax(p$u_bl_min, pmin(p$u_bl_max, (I_wma - 0)/200)),  # blind position signal. eg 36 in [1]
  g_bl = 1 + (p$g_bl_max-1)*u_bl)                              # blind g-value. eq 35 in [1]
}


#' @description Calculate local wind speeds at building height and half height. Eq. 22 in [1]
#' @export
calc_U_loc <- function(.df, p) {
  delta=p$U_delta[[1]][p$class]
  alpha=p$U_alpha[[1]][p$class]
  .df %>% mutate(
  U_rf = U10 * 1.59 * (p$H_b/delta)^alpha,       # wind speed at building height
  U_ew = U10 * 1.59 * ((0.5*p$H_b)/delta)^alpha) # wind speed at building midpoint
}

#' @description Calculate weather exposed exterior surface heat transfer coefficients. Section 3.1 in [1]
#' @export
calc_h_se <- function(.df, p, dT = 5) .df %>% mutate(
  F_ww = 0.5,                  # factor wind ward facing surface. could be weighted based on wind direction
  h_n = 1.53*abs(dT)^0.36,     # natural (due to temp diff) convection due to temp diff
  # forced (due to wind) convective on vertical and horizontal surfaces
  h_cf_ver = F_ww*(3.39 - 5.03*p$lambda_p[[1]][p$class])*U_ew^0.94 + (1 - F_ww)*(1.15 + 0.82*p$lambda_p[[1]][p$class])*U_ew^0.94,
  h_cf_hor = (3.57 + 1.72*p$lambda_p[[1]][p$class])*U_rf^0.84,
  #h_cf_ver2 = F_ww*(2.38 * U_ew^0.89) + (1-F_ww)*(2.86 * U_ew^0.617), # MoWITT
  #h_cf_hor2 = (2.38 * U_rf^0.89), # MoWITT
  # calc h_se for the elements roof, external wall and glazing. In [1] radiative part, h_re, is added later. 
  h_se_rf = ((1 - p$Fsky_hor) * p$h_re + h_n + p$Rf_rf*(sqrt(h_cf_hor^2 + h_n^2) - h_n)), 
  h_se_ew = ((1 - p$Fsky_ver) * p$h_re + h_n + p$Rf_ew*(sqrt(h_cf_ver^2 + h_n^2) - h_n)),
  h_se_gl = ((1 - p$Fsky_ver) * p$h_re + sqrt(h_cf_ver^2 + h_n^2))  # Rf_gl = 1
  ) %>% select(-F_ww, -h_cf_ver, -h_cf_hor)



#' @description Calculate infiltration potential Section 3.2 in [1]
#' @export
calc_Hmod_inf <- function(.df, p) .df %>%  mutate(
  Qmod_w = p$C_w*(0.5 * p$p_a * (p$lambda_w[[1]][p$class] * U_rf)^2)^p$n_inf,                          # wind effect
  Qmod_s = p$C_s*(p$p_a * 9.806 * p$H_mod * abs(p$T_int_set - T_e) / (p$T_int_set + 273.15))^p$n_inf,  # stact effect
  Qmod_inf = (Qmod_w^(1/p$n_inf) + Qmod_s^(1/p$n_inf) - 0.33*(Qmod_w*Qmod_s)^(1/(2*p$n_inf)))^p$n_inf, # combine wind and stack effects
  Hmod_inf = p$kp_a*Qmod_inf
  ) %>% select(-Qmod_w, -Qmod_s, -Qmod_inf)

#' @description Just a simple constant flow and balanced ventilation system for now. Section 2.3 in [1]
#' @export
calc_H_ve <- function(.df, p) .df %>% mutate(
  H_ve = p$kp_a*p$Q_ve * (1 - p$eta_ve)  # P_ve=H_ve*(T_int - T_e), eq 13 in [1]
)


# @description Calculate shading reduction factor
# @export
#calc_F_sh <- function(.df, D_ovh=2, L_ovh=1, H_gl_ovh=6, H_obst=30, L_obst=50, H_gl_obst=3*2) .df %>% mutate(
#  F_sh_ovh = 1 - pmax(0, pmin(1, (D_ovh*tan(alpha_sol) - L_ovh)/H_gl_ovh)), # shading overhangs
#  F_sh_obst = (pmin(H_obst, L_obst*tan(alpha_sol)) + H_gl_obst)/(H_gl_obst + H_obst), # shading obstacles
#  F_sh = (alpha_sol > 0) * (F_sh_ovh + F_sh_obst - 1))                      # combine shading effect


# #' @description Calculate shading reduction factor obstacles. 
#calc_F_sh_obst <- function(alpha_sol, H_obst, L_obst, H1, H0=0, F_W=0.5)
#  pmin(1, pmax(0, (H1/F_W - pmax(0, H_obst - H0 - L_obst * tan(alpha_sol)))/(H1/F_W)))
#
#calc_h_obst <- function(alpha_sol, H_obst, L_obst, H1, H0=0) H_obst - H0 - L_obst * tan(alpha_sol)

#' @description Calculate shading reduction factor overhangs. Eq 32 in [1].
calc_F_sh_ovh <- function(alpha_sol, D_ovh=2.0, L_ovh=1, H_gl_ovh=6)
  1 - pmax(0, pmin(1, (D_ovh*tan(alpha_sol) - L_ovh)/H_gl_ovh))


#' @description Calculate shading effect on horizontal (roof) and vertical (walls) surfaces. Section 3.4 in [1].
#' TODO: needs a remake and validation. use viewshed approach for easier GIS integration?
#' @export
calc_F_sh <- function(.df, p) {
  H0 <- 0
  H_b <- p$H_b
  L_obst1 <- 10*p$lambda_sol[[1]][p$class]
  H_obst1 <- 15
  L_obst2 <- H_b*p$lambda_sol[[1]][p$class]
  H_obst2 <- 1.3*H_b
  message(paste0("L_obst1: ", L_obst1, " H_obst1: ", H_obst1,  " L_obst2: ", L_obst2, " H_obst2: ", H_obst2))
  .df %>% mutate(
    tan_alpha = tan(alpha_sol),
    H_sh_obst1 = pmin(H_b, pmax(0, (H_obst1 - H0 - L_obst1 * tan_alpha))), # F.13 in [2]
    H_sh_obst2 = pmin(H_b, pmax(0, (H_obst2 - H0 - L_obst2 * tan_alpha))),
    F_sh_ovh =  calc_F_sh_ovh(alpha_sol, D_ovh=2.0),          # eq 32 in [1]
    F_sh_obst = 1 - (H_sh_obst1*0.5 + H_sh_obst2*0.5) / H_b,  # eq 33 in [1], *0.5 for semitransparency. 
    F_sh_ver = (alpha_sol > 0) * (F_sh_ovh + F_sh_obst - 1),  # eq 31 in [1]
    H_rf = pmin(3, (H_b - H_sh_obst1))/3,
    F_sh_hor = H_rf*0.5 + 0.5) %>%                            # eq 34 in [1]
    select(-tan_alpha, -H_sh_obst1, -H_sh_obst2, -F_sh_obst, -F_sh_ovh, -H_rf)
}

#' @description Internal heat gain, not part of [1]
#' @export
calc_P_int <- function(.df, p) {
  l = 2*1/(diff(as.numeric(.df$timestamp[1:2]))/3600) # roll back 2 hours
  .df %>% mutate(
    q_oc = 0*S_oc,
    q_A =  1 + 1*S_A,
    Ev_sol_wma = RcppRoll::roll_meanr(Ev_sol, weights = 1:l, fill=0),
    q_L0 = pmax(p$Ev_max - Ev_sol_wma, p$Ev_min)/p$Kv_L,
    q_L =  q_L0*S_L,
    q_WA = 0,
    q_HVAC = 0,
    q_proc = 0,
    P_int = q_oc + q_A + q_L + q_WA + q_HVAC + q_proc
  ) #%>% select(-q_oc, -q_A, -q_L, -q_WA, -q_HVAC, -q_proc)
}
