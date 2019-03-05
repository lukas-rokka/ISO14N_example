#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector calc_Trad_s(NumericVector Te, NumericMatrix xy) { //}  NumericVector tbl_Te, NumericVector tbl_Trad_s) {
  // linear interpoplation of the supply temperature from a xy-table,
  // x in xy-table required to be in increasing order
  int N = Te.size();
  NumericVector Trad_s(N);
  int Nxy = xy.nrow() - 1;

  //Rcout << Nxy << std::endl;
  //Rcout << xy << std::endl;

  for (int n=0; n < N; n++) {
    int i = Nxy; // start search from end, because more likely

    if (Te(n) >= xy(Nxy, 0)) {
      Trad_s(n) = xy(Nxy, 1); // use last y-value
    } else if ((Te[n] <= xy(0, 0))) {
      Trad_s(n) = xy(0, 1); // use first y-value
    } else {
      while (xy(i, 0) > Te[n]) i += -1;
      //Rcout << "x: " << Te(n) << std::endl; // x
      //Rcout << "xL: " << xy(i, 0) << std::endl; // xL
      //Rcout << "yR: " << xy(i+1, 1) << std::endl; // yR
      //Rcout << "xR: " << xy(i+1, 0) << std::endl; // xR
      //Rcout << "yL: " << xy(i, 1) << std::endl; // yL
      // linear interpoplations: ((x - xL)*yR + (xR - x)*yL) / (xR - xL)
      Trad_s[n] = ((Te(n) - xy(i, 0))*xy(i+1, 1) + (xy(i+1, 0) - Te(n))*xy(i, 1)) / (xy(i+1, 0) - xy(i, 0));;
    }
  }
  return Trad_s;
}
