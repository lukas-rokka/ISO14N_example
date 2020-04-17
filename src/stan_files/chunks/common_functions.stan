vector calc_T_hyd_s(vector Te, vector x, vector y) { //}  NumericVector tbl_Te, NumericVector tbl_Trad_s) {
  // linear interpoplation of the supply temperature from a xy-table,
  // x in xy-table required to be in increasing order
  int N = dims(Te)[1];
  vector[N] T_hyd_s;
  int N_xy = dims(x)[1]; // need to be at least 3
  
  int i;

  for (n in 1:N) {

    if (Te[n] >= x[N_xy]) {
      T_hyd_s[n] = y[N_xy]; // use last y-value
    } else if (Te[n] <= x[1]) {
      T_hyd_s[n] = y[1]; // use first y-value
    } else {
      i = N_xy; 
      while (x[i] > Te[n]) i += -1;
      
      T_hyd_s[n] = ((Te[n] - x[i])*y[i+1] + (x[i+1] - Te[n])*y[i]) / (x[i+1] - x[i]);
    } 
  }
    
  return T_hyd_s;
}

