//#define DEBUG 0
#include "hurdle_likelihood.h"
#include "assert.h"
using namespace arma;

const double HurdleLikelihood::large = 30;

// [[Rcpp::depends(RcppArmadillo)]]
HurdleLikelihood::HurdleLikelihood (const arma::vec& y_, const arma::mat& xd_, // const arma::mat& xc_, 
				    const arma::ivec& grp_, const arma::vec& th, const arma::vec& lambda_, double tol)
  : y(y_), yI(y_.n_rows, fill::ones),
    xd(xd_),  grp(grp_), 
    th_d(xd_.n_cols), th_c(xd_.n_cols), kbb(999),
    gpart(xd_.n_rows, fill::zeros), hpart(xd_.n_rows, fill::zeros), 
    gplusc(xd_.n_rows, fill::zeros), cumulant(xd_.n_rows, fill::zeros), cumulant2(xd_.n_rows, fill::zeros),
    small(xd_.n_rows, fill::zeros), dc(xd_.n_rows, fill::zeros),
    pen(grp_.n_rows, fill::zeros),
    Sy2(sum(square(y_))), Sy(sum(y_)), SyI(999),
  Sn(xd_.n_rows),k(xd_.n_cols), pm(parmap(k)),
  Syxd(xd_.n_cols), SyIxd(xd_.n_cols),
    pengrp(xd_.n_cols, fill::zeros), lambda(lambda_)
{
  // set up indicator matrices
   yI.elem(find(abs(y)<tol)).zeros();
   SyI = sum(yI);
   SyIxd = xd.t() * yI;
   Syxd = xd.t() * y;
  // populate params
   populatePar(th);

   DPRINT("Finished construction. dim(gpart)=" << gpart.n_rows << "gpart(0)=" << gpart(0)<<"\n");
}


void HurdleLikelihood::populatePar(const vec& th){
  th_d = th.subvec(0, k-1);
  th_c = th.subvec(k, 2*k-1);
  kbb = th(2*k);
  gpart = xd*th_d; //N x M * M * 1
  hpart = xd*th_c;
  gplusc = gpart + -.5*log(kbb/(2*datum::pi))+pow(hpart, 2)/(2*kbb);
  small = find(gplusc<large);
  cumulant = gplusc;
  cumulant.elem(small) = log(1+exp(cumulant.elem(small)));
  cumulant2.ones();
  cumulant2.elem(small) = 1/(1+exp(-gplusc.elem(small)));
}



// Replace all coordinates
double HurdleLikelihood::LL(const vec& th, bool penalize){
  populatePar(th);
  return(LLpen(penalize));
}



//th parameters, grp groups
//all others held fixed
double HurdleLikelihood::LLpen(bool penalize){
  double pen = 0;
  if(kbb<=0) return(99999999.0);
  double nloglik = sum(yI % gpart + y % hpart - cumulant) -.5*Sy2*kbb;

  // for(g=0; g<ngrp, g++){
  //   pengrp(g) = 
  // }
  // pengrp = sqrt(square(gba)+square(hba)+square(hab)+square(kba));
  
  if(penalize) pen = 1;
  double negll = -(nloglik/Sn -pen);
  //DPRINT("negll=" << negll<<"\n");
  return(negll);
}


/* Only the part in front of cumulant2 (e^gplusc/(1+e^gplusc) 
   needs to be calculate per observation.  This includes
   dgba/dtheta and dCba/dtheta, so we include those in `dc`.
   The rest collapses into a sum, so we use the linear sufficient stats,
   which are added to `Sdc` after we sum over `dc`.
 */
vec HurdleLikelihood::gradAll(const vec& th, bool penalize){
  populatePar(th);
  vec g = gradFixed(penalize);
  return(g);
}

vec HurdleLikelihood::gradFixed(bool penalize){
  arma::vec gvec(2*k+1, fill::zeros);

  gvec.subvec(0, k-1) = SyIxd- trans(trans(cumulant2) * xd);
  gvec.subvec(k, 2*k-1) = Syxd - trans(trans(cumulant2 % hpart / kbb) * xd);
  gvec(2*k) = as_scalar(trans(cumulant2) * (square(hpart)/ kbb+1)) / (2*kbb) - .5*Sy2;
  // if(false){ //penalize
  //     ASSERT_TRUE(any(abs(th)>0) || lambda(grp)<datum::eps)
  //     pen4 = th / pengrp(grp) * lambda(grp);
  //     return(-Sdc4/Sn + pen4);
  //   }
    return(-gvec/Sn);
  }

void HurdleLikelihood::setLambda(const vec& lambda_){
  lambda = lambda_;
}

// arma::mat HurdleLikelihood::hessian(const arma::ivec& idx){
//   arma::vec tmp = exp(gplusc);
//   cumulant3= 1/tmp;
//   cumulant3.elem(small) = tmp.elem(small)/square(1+tmp.elem(small));
//   x * w

// }

// map from parameters into starting indices of theta
arma::uvec::fixed<7> HurdleLikelihood::parmap(int k_) {
  arma::uvec::fixed<7> pm;
  // pm[tG0] = 0;
  // pm[tGBA] = pm[tG0]+1;
  // pm[tH0] = pm[tGBA]+k_;
  // pm[tHBA] = pm[tH0]+1;
  // pm[tHAB] = pm[tHBA]+k_;
  // pm[tKBA] = pm[tHAB]+k_;
  // pm[tK0] = pm[tKBA]+k_;
  return(pm);
}


/* 
 * R INTERFACE CODE
 */
using namespace Rcpp;

RCPP_MODULE(HurdleLikelihoodMod){
  class_<HurdleLikelihood>("HurdleLikelihoodMod")

    .constructor<arma::vec, arma::mat, arma::ivec, arma::vec, arma::vec, double>()
    .method("LLall", &HurdleLikelihood::LL, "Eval log-likelihood at a new parameter theta.")
    .method("gradAll", &HurdleLikelihood::gradAll, "Eval grad at new theta.")
    .method("gradFixed", &HurdleLikelihood::gradFixed, "Eval grad at (implicit) current theta.")
    .field("lambda", &HurdleLikelihood::lambda, "mutate penalty parameter/vector.")
    .field_readonly("cumulant2", &HurdleLikelihood::cumulant2, "Access squared cumulant.")
    .field_readonly("gpart", &HurdleLikelihood::gpart, "G crossproduct.")
    .field_readonly("hpart", &HurdleLikelihood::hpart, "H crossproduct.")
    .field_readonly("gplusc", &HurdleLikelihood::gplusc, "G + C.")
    ;
}
