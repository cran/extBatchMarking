/*** R

library(Rcpp)
library(RcppArmadillo)
library(RcppDist)

*/

//' Transition State Probability
//' 'gamma_gt' computes the transition probability matrix
//'
//' @param R integer number of marked individuals released per occasion
//' @param phi double number. Survival probability of individuals
//' @param cores The number of cores on your machine.
//' @return pR Returns the transition matrix
//' @name gamma_gt



#include <RcppArmadillo.h>
//#include <omp.h>
//#include <math.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat gamma_gt(int R, double phi){

  // Parallel
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams
  //omp_set_num_threads(cores);

  // declare our variables

  arma::mat PR = zeros(R+1, R+1);

  double f;
  double b;
  double y;

  for(int i = 0; i <= R; i++){
    for(int j = 0; j <= i; j++){
      f = lgamma(i+1) - lgamma(i-j+1) - lgamma(j+1);
      b = j*log(phi) + (i-j)*log1p(-phi);
      y = exp(f + b);
      PR(i,j) = y;
    }

  }

  return PR;

}
