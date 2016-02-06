/* Computation for hurdle model:
 * likelihood function, gradient and hessian
 * 
 */
//#define DEBUG 0
//#include "grplasso.h"
#include <RcppArmadillo.h>

class HurdleLikelihood {
 private:
  static const double large=30;
 public:
  arma::vec y, yI; 			// response
  arma::mat xd;//, xc;			// design matrix, assumed to be in group order
  arma::ivec grp;			// parameter grouping
  int ngrp; //number of groups
  arma::vec th_d, th_c;			// parameter values for discrete and continuous
  double kbb; //normal precision

  //parts of likelihood 
  //sums used repeatedly
  arma::vec gpart, hpart, gplusc, cumulant, cumulant2;
  // working variables used in inner loops
  arma::uvec small;
  arma::vec dc, pen;
  //linear sufficient statistics (that allow us to sum over the data)
  const double Sy2, Sy;
  double SyI;
  const int Sn, k;
  //mapping from variables to indices in th
  const arma::uvec pm;
  arma::vec Syxd, SyIxd;
  //group sums for penalties and scaling
  arma::vec pengrp;
  arma::vec lambda;
  void populatePar(const arma::vec& th);
  void updateGroupSums();

  //public:
  HurdleLikelihood (const arma::vec& y_, const arma::mat& xd_, const arma::ivec& grp_, const arma::vec& th, const arma::vec& lambda_, double tol);
  double LL(const arma::vec& th, bool penalize);
  double LL(bool penalize);
  arma::vec grad(const arma::vec& th, bool penalize);
  arma::vec grad(bool penalize);
  //arma::mat hessian(const arma::vec& th, int grp);
  void setLambda(const arma::vec& lambda_);
  static arma::uvec::fixed<7> parmap(int k_);

};
