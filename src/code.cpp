#include <Rcpp.h>

using namespace Rcpp;

int gcd(int a, int b) {
    if (b == 0) {
      return a; 
    }
    return gcd(b, a % b);
}

// [[Rcpp::export]]
int vector_gcd(NumericVector vec) {
    int current_gcd = vec[0];
    for (int i = 1; i < vec.size(); ++i) {
        current_gcd = gcd(vec[i], current_gcd);

        if (current_gcd == 1) {
            return (1);
        }
    }
    return current_gcd;
}

bool all_sug(LogicalVector x) {
  return is_true(all(x == TRUE));
}

int check_equal(NumericVector x, NumericVector y) {
  if (all_sug(x == y)) {
    return 1;
  } else {
    return 0;
  }
}

// [[Rcpp::export]]
int count_occurrences(DataFrame df_treatments, NumericVector curr_treatment){
  if (df_treatments.length() == 1) { // If the number of columns is 1
    return 0;
  }
  
  int count = 0;
  for (int i = 1; i < df_treatments.length(); ++i) {
    if (check_equal(df_treatments[i], curr_treatment)){
      count++;
    }
  }

  return count;
}
