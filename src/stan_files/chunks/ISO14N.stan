// Copyright (C) 2019 Lukas Lundström

//' @title Space heating use simulation
//' 
//' @description Space heating use and internal air temperature simulated using an ISO5016-1:2017 [2] based RC-network. See paper [1] for details. 
//' References:
//' 1. Lundström, L; Akander, J; Zambrano, J. Development of a Space Heating Model Suitable for the Automated Model Generation of Existing Multifamily Buildings — A Case Study in Nordic Climate; Energies 2019. https://doi.org/10.3390/en12030485
//' 2. ISO 52016-1:2017: Energy Performance of Buildings — Energy Needs for Heating and Cooling, Internal Temperatures and Sensible and Latent Heat Loads—Part1: Calculation Procedures
//'
//' @param r_el Ratios of interior surface area of each building element type and the total floor area. Vector of 5 real values, one value per building element type. The unit is m2_s/m2_fl, ie interior surface area per floor area. 
//' @param cl_el Distribution class for each building element. Array of 5 integers, one value per building element type (values for gl and im are not used, but need to be given to fill the array). See table 1 in [1].
//' @param U_el Thermal transmittance, W/(K*m2_s). Vector of 5 real values, one value per building element type.
//' @param U_gr_vi Thermal transmittance for the virtual ground layer, W/(K*m2_s). See eq. 8 [1]
//' @param H_tb Floor area normalized heat transfer coefficient for thermal bridges, W/(K*m2_fl)
//' @param dt time interval, in hours (or seconds if heat capacities are given in Joule)
//' @param k_m Areal heat capacity of the building elements, Wh/(K*m2_s) (or J/(K*m2_s) if time interval is given in seconds). Vector of 5 real values (one value per building element type)
//' @param C_int heat capacity the internal air node, Wh/(K*m2_fl)  
//' @param k_gr Areal heat capacity of the virtual ground layer, Wh/(K*m2_s)
//' @param H_hyd Floor area normalized radiator constant (aka heat transfer coef), W/(K*m2_fl). See eq. 9 [1]
//' @param n_hyd Radiator exponent, unitless. See eq. 9 [1]
//' @param T_trv_pb Proportional band of the TRV's, °C. See eq. 12 [1]
//' @param T_int_set Set-point temperature for the  TRV's (desired indoor temperature), °C. See eq. 12 [1]
//' @param a_hyd Empirical constant used to calculate the return temperature of the hydronic heating system. See eq. 11 [1]
//' @param b_hyd Empirical constant used to calculate the return temperature of the hydronic heating system. See eq. 11 [1]
//' @param T_hyd_s Supply temperatures of the hydronic heating system, °C. N length vector.
//' @param f_c_hyd Convective fraction of the hydronic heating system, unitless. See eq. 17 [1]
//' @param f_c_int Convective fraction of the internal heat gains, unitless. See eq. 17 [1]
//' @param f_c_sol Convective fraction of the solar heat gains, unitless. See eq. 17 [1]
//' @param Fsky_ver Sky view factor vertical elements (ew & gl), eq 16 in [1]
//' @param Fsky_hor Sky view factor horizontal elements (rf), eq 16 in [1]
//' @param U10_idx Discretized wind speeds. Integer array of same length as data input. 
//' @param h_se Heat transfer coefficients for weather-exposed exterior surfaces, W/(K*m2_s). Matrix of N_U10 x 3 size, where row N_U10 is the number of discretized wind speeds (eg. 10) and each column give h_se values for the building elements roof, external walls and glazing.
//' @param T_e External air temperature, °C. N length vector.
//' @param T_gr Ground temperature, °C. N length vector.
//' @param T_sky Sky temperature, °C. N length vector.
//' @param I_tot_ver_sh Orientation weighted hemispherical solar irradiance on vertical surfaces including shading effects, W/m2_s. N length vector. See eq. 28 [1]. 
//' @param I_tot_hor_sh Hemispherical solar irradiance on horizontal surfaces (the roof) including shading effects, W/m2_s. N lenght vector. See eq. 30 [1].
//' @param P_gn_sol Floor area weighted solar heat gains, including shading and window blinds, W/m2_fl. N lenght vector. See eq. 27 [1].
//' @param P_gn_int Internal heat gains (excl. solar heat gains), W/m2_fl. N length vector.
//' @param H_ve Floor area normalized heat transfer coefficient for ventilation, W/(K*m2_fl). N lenght vector.
//' @param H_inf Floor area normalized heat transfer coefficient for infiltration, W/(K*m2_fl). N length vector.
//' @param N_pl Number of planes for the opaque elements (rf, ew, gf). N_pl = 3, result in a 14 nodes system described in [1], N_pl = 5 results in a 20 nodes system.
//' @param debug Values > 0 prints additional information. 
//' @param Nout Outputs to return. 1: only P_hyd, 2: P_hyd & T_int, 3: All node temperatures and additonal variables
//' @param recalc_h_re If > 0 additional calculation to account for non-linear surface radiation effects are conducted (mainly impact during extreme cold or extreme hot periods when h_re = 4.14 is a bad approximation)

matrix ISO14N(
  vector r_el,
  int[] cl_el,
  vector U_el,
  real U_gr_vi, 
  real H_tb,
  real dt, 
  vector k_m, real C_int, real k_gr,
  real H_hyd, real n_hyd, real T_trv_pb, real T_int_set, real a_hyd, real b_hyd, vector T_hyd_s,
  real f_c_hyd, real f_c_int, real f_c_sol,
  real Fsky_ver, real Fsky_hor, 
  int[] U10_idx, matrix h_se,
  vector T_e, vector T_gr, vector T_sky,
  vector I_tot_ver_sh, vector I_tot_hor_sh, vector P_gn_sol,
  vector P_gn_int, 
  vector H_ve, vector H_inf,
  int N_pl, int debug, int Nout, int recalc_h_re) 
{
  /*****************************************************
  //  Declarations
  *****************************************************/
  
  /* Counts and indices */
  int N_U10 = rows(h_se);         // number of U10 categories
  int N = rows(T_e);              // total number of time steps
  int N_el = 5;                   // number of elements
  int pln[N_el] = {N_pl, N_pl, 2, 2, N_pl}; // number of planes/nodes per element
  int m = sum(pln) + 1;           // total number of nodes / states
  int rf = 1;                     // index element roof 
  int ew = 2;                     // index element external walls
  int gl = 3;                     // index element glazing 
  int im = 4;                     // index element internal mass
  int gf = 5;                     // index element ground floor
  int idx_se[N_el];               // indices for exterior surface nodes
  int idx_si[N_el];               // indices for interior surface nodes
  
  /* Ratios and fractions */
  //real r_im1 = 1.5;               // ratio internal walls per square meter floor, eq 3 in [1]
  //real r_im2 = r_el[im] - r_im1;   // ratio internal ceilings/roof per square meter floor
  vector[N_el] f_el = r_el / sum(r_el); // area fractions of the building elements, eq 17 in [1]
  real f_r_int  = 1.0 - f_c_int;  // Radiative fraction internal heat gains, eq 17 in [1]
  real f_r_sol  = 1.0 - f_c_sol;  // Radiative fraction solar heat gains
  real f_r_hyd  = 1.0 - f_c_hyd;  // Radiative fractions hydronic heating system
  real a_sol = 0.5;               // solar absorption coefficient
  
  /* Heat transfer coefs, resistances and heat capacities */
  row_vector[N_pl-1] H_el[N_el];     // heat transfer coefs, W/(K*m2_fl). eq 5 in [1]
  real h_ri = 5.13;               // radiative heat transfer coef at 20C, W/(K*m2). Table 25 in [2]
  real h_re = 4.14;               // radiative heat transfer coef at 0C, W/(K*m2). Table 25 in [2]
  // convective heat transfer coefs interior surfaces, eq. 4 in [1]. Surface are weighted: r*h_ri, W/(K*m2_fl)
  vector[N_el] H_ci = [           
    r_el[rf]*0.7,   // Downwards 
    r_el[ew]*2.5,   // Horizontal 
    r_el[gl]*2.5,   // Horizontal 
    r_el[im]*2.5,   // Average (could be weighted on number of floors)
    r_el[gf]*5.0]'; // Upwards
  vector[N_el] H_ri = r_el * h_ri; // radiative heat transfer coefs, W/(K*m2_fl). eq 2 in [1]
  vector[N_el] R_el;// = [1/U_rf, 1/U_ew, 1/U_gl, 1.0, 1/U_gf]' - 0.13 - 0.04; // resistances, eq 6 in [1]
  real R_gr = 0.5/2.;              // 0.5 meter soil with conductivity=2.0, eq 8 in [1]
  real H_gr_vi = r_el[gf]*U_gr_vi; // See eq. 8 in [1]
  // distribution of the heat capacities for roof and external walls. See table 1 in [1] 
  vector[N_pl] k_dist[5] = (N_pl == 3) ? {
    [0.10, 0.40, 0.50]', // Class I
    [0.50, 0.40, 0.10]', // Class E
    [0.40, 0.20, 0.40]', // Class IE
    [0.25, 0.50, 0.25]', // Class D ([0.33, 0.33, 0.33] used in paper [1])
    [0.10, 0.80, 0.10]'} : { // Class M
    [0.00, 0.00, 0.00, 0.00, 1.00]', // Class I
    [1.00, 0.00, 0.00, 0.00, 0.00]', // Class E
    [0.50, 0.00, 0.00, 0.00, 0.50]', // Class IE
    [0.125,0.25, 0.25, 0.25, 0.125]',// Class D 
    [0.00, 0.00, 1.00, 0.00, 0.00]'};// Class M
  // distribution of the heat capacities for ground floor. 
  vector[N_pl] k_dist_gf[5] = (N_pl == 3) ? {
    [0.00, 0.40, 0.60]', // Class I ([0, 0.5, 0.5] used in paper [1])
    [0.20, 0.70, 0.10]', // Class E
    [0.20, 0.40, 0.40]', // Class IE
    [0.20, 0.40, 0.40]', // Class D
    [0.20, 0.70, 0.10]'} : { // Class M
    [0.00, 0.00, 0.00, 0.00, 1.00]', // Class I
    [0.00, 0.00, 1.00, 0.00, 0.00]', // Class E
    [0.00, 0.00, 0.50, 0.00, 0.50]', // Class IE
    [0.00, 0.00, 0.25, 0.50, 0.20]',// Class D ([0.33, 0.33, 0.33] used in [1])
    [0.00, 0.00, 0.00, 1.00, 0.00]'};// Class M
  vector[m] C;                   // heat capacity of opaque element [Wh/(K*m2_fl)]
  vector[m] Cdt;                 // C/dt 

  /* Hydronic heatings system */
  vector[N] T_hyd_r;            // return temperature
  vector[N] T_hyd_lmtd;         // logarithmic mean temperature difference
  vector[N] P_hyd;              // heat output hydronic heating system
  vector[N] u_trv;              // control signal TRV's
  
  /* Equation system */
  vector[m] X;                  // Node temperatures
  vector[m] b;                  // known terms
  matrix[m, m] A = rep_matrix(0, m, m);  // system coefs
  matrix[m, m] A_inv[N_U10];             // Array of inverted A
  //matrix[m, 12] B = rep_matrix(0, m, 12);
  //matrix[m, 12] B_arr[N_U10];
  matrix[N, (Nout > 2) ? (m + 5) : Nout] out; // returned values
  
  /* temporary variables */
  vector[N_pl] r_init;// used for first guess on node states
  vector[N] P_r;      // radiative heat gains
  real T_trv_set = T_int_set + T_trv_pb/2;
  real dT;            // temperature difference;
  real tmp;           // temporary variable
  real h_se_add;      // correction term
  real h_re_new;      // correction term 
  real T_h_r_sky;     // temporary variable

  // calc index postion for interior & exterior surface nodes
  idx_si[1] = 0;
  for (el in 1:N_el) {
    idx_si[el] = idx_si[max(el-1, 1)] + pln[el];
    idx_se[el] = idx_si[el] - pln[el] + 1;
  }
  

  /*****************************************************
  //  Distribution of thermal mass and resistance
  *****************************************************/
  
  R_el = 1.0 ./ U_el - 0.13 - 0.04; // resistances, eq 6 in [1]

  if (N_pl == 3) {  // 3 planes in the opaque elements rf, ew, gf. A total of 14 nodes
    H_el[rf] = r_el[rf]*[2, 2] / R_el[rf]; 
    H_el[ew] = r_el[ew]*[2, 2] / R_el[ew];  
    H_el[gl] = r_el[gl]*[1, 0] / R_el[gl];
    H_el[im] = r_el[im]*[1, 0] / R_el[im];
    H_el[gf] = r_el[gf]*[1/(R_el[gf]/2 + R_gr), 1/(R_el[gf]/2)];
    
    C = append_row(r_el[rf]*k_m[rf]*k_dist[cl_el[rf]],  // rf
        append_row(r_el[ew]*k_m[ew]*k_dist[cl_el[ew]],  // ew
        append_row(r_el[gl]*k_m[gl]*[0.5, 0.5]',        // gl
        append_row(r_el[im]*k_m[im]*[0.9, 0.1]',        // im
        append_row(r_el[gf]*k_m[gf]*k_dist[cl_el[gf]],  // gf
        C_int)))));                                     // internal air
    C[idx_se[gf]] += r_el[gf]*k_gr;   // add virtual ground layer mass
    
    r_init =[1, 0.5, 0]';  // initial node temp weighting
    
  } else if (N_pl == 5) {  // 5 planes in the opaque elements rf, ew, gf. A total of 20 nodes
    vector[5] gf_dist;
    H_el[1] = r_el[1]*[6, 3, 3, 6] / R_el[rf]; 
    H_el[2] = r_el[2]*[6, 3, 3, 6] / R_el[ew]; 
    H_el[3] = r_el[3]*[1/R_el[gl], 0, 0, 0];
    H_el[4] = r_el[4]*[1/R_el[im], 0, 0, 0];
    H_el[5] = r_el[5]*[2/R_gr,  1/(R_el[gf]/4 + R_gr/2), 2/R_el[gf], 4/R_el[gf]]; // iso52016 eq 48
    
    C = append_row(r_el[rf]*k_m[rf]*k_dist[cl_el[rf]],  // rf
        append_row(r_el[ew]*k_m[ew]*k_dist[cl_el[ew]],  // ew
        append_row(r_el[gl]*k_m[gl]*[0.5, 0.5]',        // gl
        append_row(r_el[im]*k_m[im]*[0.9, 0.1]',        // im
        append_row(r_el[gf]*k_m[gf]*k_dist[cl_el[gf]],  // gf
        C_int)))));                                     // internal air
    C[idx_se[gf]+1] += r_el[gf]*k_gr;  // add virtual ground layer mass
        
    r_init = [1, 0.75, 0.5, 0.25, 0]'; // initial node temp weighting
  } else reject("Only 3 or 5 planes supported")
  

  
  /*****************************************************
  //  Construct system matrix A
  //  1 thermal zone with 5 elements
  //  X = A \ b can be calculated by X = inverse(A) * b
  *****************************************************/

  // build diagonal and coupling nodes
  for (el in 1:N_el) {
    // Exterior side
    A[idx_se[el], idx_se[el]]   =  H_el[el, 1]; // diagonal/capacity node
    A[idx_se[el], idx_se[el]+1] = -H_el[el, 1]; // coupling next node
    // middle/inside nodes
    if (pln[el] > 2) {
      for (pl in 2:(pln[el]-1)) {
        int idx_diag = idx_se[el]+pl-1;
        A[idx_diag, idx_diag-1] = -H_el[el, pl-1]; // coupling previous node
        A[idx_diag, idx_diag] =  + H_el[el, pl-1] + H_el[el, pl]; // diagonal/capacity node
        A[idx_diag, idx_diag+1] = -H_el[el, pl];  // coupling next node
      }
    }
    // Interior side node
    A[idx_si[el], idx_si[el]-1] = -H_el[el, pln[el]-1]; // coupling previous node
    A[idx_si[el], idx_si[el]]   =  H_el[el, pln[el]-1] + H_ci[el] + (1-f_el[el])*H_ri[el]; // diagonal/capacity node
    
    // interior surface radiative heat exchange
    for (j in 1:N_el) if (j != el) A[idx_si[el], idx_si[j]] = -H_ri[el]*f_el[j];
  }

  // Internal air thermal zone node
  A[m, m] = sum(H_ci) + H_tb;

  // interior surface convective heat exchange
  A[m, idx_si] = -H_ci';
  A[idx_si, m] = -H_ci;

  // add constant H_se to exterior node ground floor (gf)
  A[idx_se[5], idx_se[5]] += H_gr_vi;
  
  // precalculate inverse A matrises for every U10 category, using
  // variable H_se for roof, external walls and glazing
  Cdt = C/dt;
  A += diag_matrix(Cdt);
  for (i in 1:N_U10) {
    A[idx_se[rf], idx_se[rf]] = Cdt[idx_se[rf]] + H_el[rf, 1] + r_el[rf]*(h_se[i, rf] + Fsky_hor*h_re);
    A[idx_se[ew], idx_se[ew]] = Cdt[idx_se[ew]] + H_el[ew, 1] + r_el[ew]*(h_se[i, ew] + Fsky_ver*h_re);
    A[idx_se[gl], idx_se[gl]] = Cdt[idx_se[gl]] + H_el[gl, 1] + r_el[gl]*(h_se[i, gl] + Fsky_ver*h_re);
    A_inv[i] = inverse_spd(A);
  }
  
  /*****************************************************
  // Initiation
  *****************************************************/


  if (debug > 3) {
    for (i in 1:m) print(round(A[i, 1:m]*100)/100);
    print("m: ", m, ", N: ", N);
    print("idx_si:", idx_si);
    print("idx_se:", idx_se);
    print("r_el: ", r_el);
    print("f_el: ", f_el);
    print("f_r_sol: ", f_r_sol, " f_r_int: ", f_r_int, " f_r_hyd: ", f_r_hyd);
    print("H_ci: ", H_ci);
    print("H_el: ", H_el);
    print("C: ", C);
  }
  
  // initial node states
  tmp = mean(T_e[1:24]);
  X[idx_se[1]:idx_si[1]] = tmp*r_init + (1-r_init)*T_int_set;
  X[idx_se[2]:idx_si[2]] = tmp*r_init + (1-r_init)*T_int_set;
  X[idx_se[3]:idx_si[3]] = [tmp, T_int_set]';
  X[idx_se[4]:idx_si[4]] = [T_int_set, T_int_set]';
  X[idx_se[5]:idx_si[5]] = T_gr[1]*r_init + (1-r_init)*T_int_set;
  X[m] = T_int_set;
  u_trv[1] = 0.8;


  /*****************************************************
  // start the simulation
  *****************************************************/

  for (t in 1:N) {

    //if (debug > 10 && t>1 && t < debug && fmod(t, 1) == 0) { print("t: ", t, " X[1]: ", X[1]);  }
    
    
    /*****************************************************
    // Hydronic heating system
    //
    *****************************************************/


    if (T_hyd_s[t] <= X[m]) { // no heating need
      T_hyd_r[t] = T_hyd_s[t];
      T_hyd_lmtd[t]= 0;
      P_hyd[t] = 0;
      u_trv[t] = 0;
    } else { // heating needed
      /* calc radiators "free" heat output, without impact from TRV */
      T_hyd_r[t] = T_hyd_s[t] - b_hyd*fmax(0, T_hyd_s[t] - X[m])^a_hyd;
      T_hyd_lmtd[t]= (T_hyd_s[t] - T_hyd_r[t]) / log((T_hyd_s[t] - X[m]) / (T_hyd_r[t] - X[m]));
      P_hyd[t] = H_hyd*T_hyd_lmtd[t]^n_hyd;
      
      /* TRV (thermostatic radiator valves) calculation */
      //u_trv[t] = fmax(0, fmin(1, (T_trv_set - X[m])/T_trv_pb)); 
      u_trv[t] = 0.5 + sin(3.1416*fmax(-0.5, fmin(0.5, (T_int_set - X[m])/T_trv_pb)))*0.5;  // IDA ICE's smooth version
      
      // update radiator heat output with TRV impact
      P_hyd[t] = u_trv[t]*P_hyd[t];
    }
  
    
    /*****************************************************
    // Construct b vector
    *****************************************************/

    // precalculations 
    P_r[t] = f_r_sol*P_gn_sol[t] + f_r_int*P_gn_int[t] + f_r_hyd*P_hyd[t]; // radiative heat towards interior surfaces
    T_h_r_sky = T_sky[t]*h_re; // sky_temp*h_re
  
    //the dependency on previous time step
    b = Cdt .* X; // C * X_k-1 / dt
    
    // add input to extorior and interior nodes
    
    // Roof
    b[idx_se[rf]] += r_el[rf]*(h_se[U10_idx[t], rf]*T_e[t] + I_tot_hor_sh[t]*a_sol + Fsky_hor*T_h_r_sky);
    b[idx_si[rf]] += f_el[rf]*P_r[t];

    // Opaque external walls
    b[idx_se[ew]] += r_el[ew]*(h_se[U10_idx[t], ew]*T_e[t] + I_tot_ver_sh[t]*a_sol + Fsky_ver*T_h_r_sky);
    b[idx_si[ew]] += f_el[ew]*P_r[t];

    // Window glazing
    b[idx_se[gl]] += r_el[gl]*(h_se[U10_idx[t], gl]*T_e[t] + Fsky_ver*T_h_r_sky);
    b[idx_si[gl]] += f_el[gl]*P_r[t];

    // Internal mass
    b[idx_si[im]] += f_el[im]*P_r[t];
    
    // Ground floor
    b[idx_se[gf]] += H_gr_vi*T_gr[t];
    b[idx_si[gf]] += f_el[gf]*P_r[t];
    
    // Internal air node
    b[m] += H_tb*T_e[t] - (H_ve[t] + H_inf[t])*(X[m] - T_e[t]) +
     f_c_int*P_gn_int[t] + f_c_sol*P_gn_sol[t] + f_c_hyd*P_hyd[t];
  
    // Compensate for non-linear h_re values
    if (recalc_h_re>0) {
      // compensate for the constant h_re = 4.14 in A matrix. add -0.06 per degree 
      // deviation from 0C, multiplied 0.5 to average two temperatures, 0.06*0.5=0.03.
      b[idx_se[rf]] -= r_el[rf]*Fsky_hor*0.030*(X[idx_se[rf]] + T_sky[t])*(X[idx_se[rf]] - T_sky[t]);
      b[idx_se[ew]] -= r_el[ew]*Fsky_ver*0.030*(X[idx_se[ew]] + T_sky[t])*(X[idx_se[ew]] - T_sky[t]);

      // additional natural surface heat transfer due to temperature differences
      // larger than what is accounted for in the precalculate h_se values
      // h_se_add = sqrt((h_se - (1-Fsky)*h_re_old)^2 - h_n_old^2 + h_n^2) - h_se + (1-Fsky)*h_re_new
      // for dT_old=5: h_n_old^2 = (1.53*5^0.36)^2 = 2.34*5^72 = 7.46
  
      dT = (T_e[t] - X[idx_se[1]]);  // temp diff external air and surface node of roof
      if (fabs(dT) > 100) {
        // h_re given at 0C, add 0.06 per degree deviation
        h_re_new = 4.14 + 0.02*(X[idx_se[rf]] + X[idx_se[ew]] + T_e[t]); 
        
        h_se_add = sqrt((h_se[U10_idx[t], 1] - (1-Fsky_hor)*4.14)^2 - 7.46 + 2.34*fabs(dT)^0.72) - h_se[U10_idx[t], 1] + (1-Fsky_hor)*h_re_new;
        b[1] += r_el[1]*h_se_add*dT;
        // do same for external walls. Roof almost always higer dT, so the check for
        // external walls can be within the conditional check of the roof element
        dT = (T_e[t] - X[idx_se[2]]);
        if (fabs(dT) > 10) {
          h_se_add = sqrt((h_se[U10_idx[t], 2] - (1-Fsky_ver)*4.14)^2 - 7.46 + 2.34*fabs(dT)^0.72) - h_se[U10_idx[t], 2] + (1-Fsky_ver)*h_re_new;
          b[idx_se[2]] += r_el[2]*h_se_add*dT;
        }
      }
    }
    
    
    /*****************************************************
    // Time update (the effect of system dynamics) 
    *****************************************************/

    // solve the equation system. A_inv[U10_idx[t]] are pre-calculated
    X = (A_inv[U10_idx[t]] * b);
    // X = mdivide_left_spd(A, b);   // higher precision than inverse(A)*b but slower
    // X = mdivide_left_ldlt(??);    // mdivide_left_ldlt can be accessed on C++ level only
    
    
    /*****************************************************
    // Measurements update (the effect of measurements)
    *****************************************************/
    
    
    /*****************************************************
    // Output
    *****************************************************/
    
    if (Nout==1) {
      out[t, 1] = P_hyd[t];
    } else if (Nout==2) {
      out[t, 1] = P_hyd[t];
      out[t, 2] = X[m];
    } else {
      out[t, 1:m] = X';
      out[t, m+1] = P_hyd[t];                     // heat load radiators
      out[t, m+2] = u_trv[t];                     // control signal radiators
      out[t, m+3] = dot_product(X[idx_si], f_el); // mean radiant temp, T_mrt
      out[t, m+4] = (out[t, m+3] + out[t, m])/2;  // operative temp, T_op
      //out[t, m+5] = X[idx_se[2]];
    }

  }
  return out;
}
