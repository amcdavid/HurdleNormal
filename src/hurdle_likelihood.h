#include <RcppArmadillo.h>
//#define DEBUG 0
#include "assert.h"
#define LARGE 25
// class HurdleLikelihood
// data x, y
// Vector th
// keep 

class HurdleLikelihood {
  public:
  arma::vec y, yI; 			// response
  arma::mat x, xI;			// covariates
  arma::ivec grp;			// parameter grouping
  double g0, h0, k0;
  arma::vec gba, hba,  hab, kba;			// current parameter values. gba, hba, hab, kba should be length max(grp)

  //parts of likelihood 
  //crossproducts with covariates
  arma::mat xI_gba, x_hba, xI_hab, x_kba;
  //sums used repeatedly
  arma::vec gpart, hpart, cpart, gplusc, cumulant, cumulant2;
  //linear sufficient statistics (that allow us to sum over the data)
  double Sy2, Sy, SyI;
  int Sn, k;
  arma::vec SxI, Sx, SyIxI, SyIx, SyxI, Syx;
  //group sums for penalties and scaling
  arma::vec pengrp, lambda;
  bool populatePar(int grp, arma::vec th);
  void updateCrossproducts(int grp);
  void updateGroupSums(bool updateSums, bool updateGrad);

  //public:
  HurdleLikelihood (arma::vec y_, arma::mat x_, arma::ivec grp_, arma::vec th, arma::vec lambda_, double tol);
  double LL(arma::vec th, int grp);
  double LL(arma::vec th);
  arma::vec grad(arma::vec th, int grp, bool penalize);
};


using namespace Rcpp;
 /// create an external pointer to a HurdleLikelihood
RcppExport SEXP HurdleLikelihood__new(SEXP y, SEXP x, SEXP grp, SEXP th, SEXP lambda, SEXP tol) {
// convert inputs to appropriate C++ types
  DPRINT("Starting wrapper\n");
  arma::vec y_ = as<arma::vec>(y);
  arma::vec th_ = as<arma::vec>(th);
  arma::mat x_ = as<arma::mat>(x);
  arma::ivec grp_ = as<arma::ivec>(grp);
  arma::vec lambda_ = as<arma::vec>(lambda);
  double tol_ = as<double>(tol);
  DPRINT("Finished wrapper\n");
  Rcpp::XPtr<HurdleLikelihood> ptr( new HurdleLikelihood( y_, x_, grp_, th_, lambda_, tol_), true );

// return the external pointer to the R side
return ptr;

}

/// invoke the loglikelihood method
RcppExport SEXP HurdleLikelihood__LL(SEXP xp, SEXP thR_, SEXP grpR_) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  // convert the parameter to int
  arma::vec th = as<arma::vec>(thR_);
  int grp = as<int>(grpR_);
  double res = ptr->LL( th , grp);
  // return the result to R
  return wrap(res);
}

// invoke the loglikelihood (don't curry)
RcppExport SEXP HurdleLikelihood__LLall(SEXP xp, SEXP thR_) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  arma::vec th = as<arma::vec>(thR_);
  double res = ptr->LL( th );
  return wrap(res);
}


/// invoke the gradient method (piecewise)
RcppExport SEXP HurdleLikelihood__grad(SEXP xp, SEXP th_, SEXP grp_, SEXP penalize_) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  arma::vec th = as<arma::vec>(th_);
  int grp = as<int>(grp_);
  bool penalize = as<bool>(penalize_);
  arma::vec res = ptr->grad(th , grp, penalize);
  return wrap(res);
}

//debugging code follows
RcppExport SEXP HurdleLikelihood__gpart(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->gpart);
}

RcppExport SEXP HurdleLikelihood__xI_gba(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->xI_gba);
}

RcppExport SEXP HurdleLikelihood__cumulant(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->cumulant);
}

RcppExport SEXP HurdleLikelihood__gplusc(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->gplusc);
}
