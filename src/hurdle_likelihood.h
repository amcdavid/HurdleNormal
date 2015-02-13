/* Computation for hurdle model:
 * likelihood function, gradient and hessian
 * 
 */
//#define DEBUG 0
#include "grplasso.h"

class HurdleLikelihood : public PenalizedLikelihood {
 private:
  static const int large=30;
  static const int  G0=0, H0=1, K0=2;
  static const int GBA = 0, HBA = 1, HAB =2, KBA=3;
  public:
  arma::vec y, yI; 			// response
  arma::mat x, xI;			// covariates
  arma::ivec grp;			// parameter grouping
  double g0, h0, k0;
  arma::vec gba, hba,  hab, kba;			// current parameter values. gba, hba, hab, kba should be length max(grp)

  //parts of likelihood 
  //sums used repeatedly
  arma::vec gpart, hpart, gplusc, cumulant, cumulant2;
  // working variables used in inner loops
  arma::uvec small;
  arma::vec dc;
  arma::vec3 Sdc3;
  const arma::vec3 pen3;
  arma::vec4 Sdc4, pen4;
  //linear sufficient statistics (that allow us to sum over the data)
  const double Sy2, Sy;
  double SyI;
  const int Sn, k;
  arma::vec SxI, Sx, SyIxI, SyIx, SyxI, Syx;
  //group sums for penalties and scaling
  arma::vec pengrp;
  arma::vec lambda;
  bool populatePar(int grp, const arma::vec& th);
  void updateCrossproducts(int grp);
  void updateGroupSums(bool updateSums);

  //public:
  HurdleLikelihood (const arma::vec& y_, const arma::mat& x_, const arma::ivec& grp_, const arma::vec& th, const arma::vec& lambda_, double tol);
  double LL(const arma::vec& th, int grp);
  double LL(const arma::vec& th);
  arma::vec grad(const arma::vec& th, int grp, bool penalize);
  arma::mat hessian(const arma::vec& th, int grp);
  void setLambda(const arma::vec& lambda_);
};
