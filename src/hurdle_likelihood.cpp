//#define DEBUG 0
#include "hurdle_likelihood.h"
#include "assert.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
HurdleLikelihood::HurdleLikelihood (const arma::vec& y_, const arma::mat& x_, const arma::ivec& grp_, const arma::vec& th, const arma::vec& lambda_, double tol)
  : y(y_), yI(y_.n_rows, fill::ones),
    x(x_), xI(x_.n_rows, x_.n_cols,fill::ones),
    grp(grp_), g0(0), h0(0), k0(1),
    gba(x_.n_cols, fill::zeros), hba(x_.n_cols, fill::zeros), hab(x_.n_cols, fill::zeros), kba(x_.n_cols, fill::zeros),
    gpart(x_.n_rows, fill::zeros), hpart(x_.n_rows, fill::zeros),
    gplusc(x_.n_rows, fill::zeros), cumulant(x_.n_rows, fill::zeros), cumulant2(x_.n_rows, fill::zeros),
    small(x_.n_rows, fill::zeros), dc(x_.n_rows, fill::zeros),
    Sdc3(fill::zeros), pen3(fill::zeros), Sdc4(fill::zeros), pen4(fill::zeros),
    Sy2(sum(square(y_))), Sy(sum(y_)), SyI(999),
  Sn(x_.n_rows),k(x_.n_cols), pm(parmap(k)),
  grpm3(fill::zeros), grpm4(fill::zeros),
  SxI(x_.n_cols), Sx(x_.n_cols), SyIxI(x_.n_cols), SyIx(x_.n_cols),
  SyxI(x_.n_cols), Syx(x_.n_cols),
    pengrp(x_.n_cols, fill::zeros), lambda(lambda_)
{
  // set up indicator matrices
  yI.elem(find(abs(y)<tol)).zeros();
  xI.elem(find(abs(x)<tol)).zeros();
  //for(uword i=0;i<x.n_elem;i++) if(std::abs(x(i))<tol) xI(i)=0; 
  // set up linear sufficient statistics
   SxI = sum(xI, 0).t();
   Sx = sum(x, 0).t();
   SyI = sum(yI);
   SyIxI = xI.t() * yI;
   SyIx = x.t() * yI;
   SyxI = xI.t() * y;
   Syx = x.t() *y;
   //set up coordinate maps
   grpm3(G0) = pm(tG0);grpm3(H0) = pm(tH0);grpm3(K0) = pm(tK0);
   grpm4(GBA)=pm(tGBA); grpm4(HBA)=pm(tHBA); grpm4(HAB)=pm(tHAB); grpm4(KBA)=pm(tKBA);
  // populate params
   populatePar(th);

    DPRINT("Finished construction. dim(gpart)=" << gpart.n_rows << "gpart(0)=" << gpart(0)<<"\n");
}

  // pack parameters group `grp` with th
//do online updates of crossproducts, gpart and hpart 
  bool HurdleLikelihood::populatePar(int grp, const vec& th){
    //    DPRINT("grp=" << grp << " th=" << th << "\n");
    bool changed=true;
    //update param and cross products
    if(grp == -1){ //intercept
      changed = std::abs(g0-th(G0))> datum::eps || 
      	std::abs(h0-th(H0))>datum::eps || std::abs(k0-th(K0))>datum::eps;
      if(changed){
	//	DPRINT("gpart(0)=" << gpart(0));
	gpart += (th(G0)-g0);
	hpart += (th(H0)-h0);
	g0 = th(G0), h0=th(H0), k0=th(K0);
	//	DPRINT("gpart(0)=" << gpart(0));
      }
    } else{ //covariates
      changed = std::abs(gba(grp)-th(GBA))> 0 || 
      std::abs(hba(grp)-th(HBA))>0 || 
      std::abs(hab(grp)-th(HAB))>0 ||
      std::abs(kba(grp)-th(KBA))>0;
	if(changed){
	  DPRINT("gpart(0)=" << gpart(0));
	  gpart += (2*xI.col(grp)*(th(GBA)-gba(grp)) + 
		    x.col(grp)*(th(HBA)-hba(grp)));
	  hpart +=  (xI.col(grp)*(th(HAB)-hab(grp)) - 
		     x.col(grp)*(th(KBA)-kba(grp)));
	  // //	  DPRINT("gpart(0)=" << gpart(0));
	  gba(grp)=th(GBA), hba(grp)=th(HBA), hab(grp)=th(HAB), kba(grp)=th(KBA);
	  pengrp(grp) = sqrt(sum(square(th)));
	  DPRINT("Updated\n");
	} else{
	  DPRINT("Unchanged\n");
	}
    } //end covariates
    DEXEC(double gbug = g0+sum(2*gba.t() % xI.row(0) + hba.t() % x.row(0)));
    DEXEC(double hbug = h0+sum(hab.t() % xI.row(0) - kba.t() % x.row(0)));
    //    DPRINT("gpart(0)=" << gpart(0)<< " gbug="<<gbug<<"\n");
    //DPRINT("hpart(0)=" << hpart(0)<< " hbug="<<hbug<<"\n");
    ASSERT_TRUE(std::abs(gpart(0)- gbug)<.0001)
    ASSERT_TRUE(std::abs(hpart(0)-hbug)<.0001)

    return(changed);
  }

void HurdleLikelihood::populatePar(const vec& th){
  g0 = th(pm[tG0]);
  gba = th.subvec(pm[tGBA],pm[tGBA]+k-1);
  h0 = th(pm[tH0]);
  hba = th.subvec(pm[tHBA], pm[tHBA]+k-1);
  hab = th.subvec(pm[tHAB], pm[tHAB]+k-1);
  kba = th.subvec(pm[tKBA], pm[tKBA]+k-1);
  k0 =th(pm[tK0]);
  pengrp = sqrt(square(gba)+square(hba)+square(hab)+square(kba));
  gpart = g0+2*xI*gba + x*hba;
  hpart = h0+xI*hab -x*kba;
}



void HurdleLikelihood::updateGroupSums(){
    gplusc = gpart + -.5*log(k0/(2*datum::pi))+pow(hpart, 2)/(2*k0);
    small = find(gplusc<large);
    cumulant = gplusc;
    cumulant.elem(small) = log(1+exp(cumulant.elem(small)));
    cumulant2.ones();
    cumulant2.elem(small) = 1/(1+exp(-gplusc.elem(small)));
}

// Replace all coordinates
double HurdleLikelihood::LL(const vec& th, bool penalize){
  populatePar(th);
  return LL(th, -2, penalize);
}

//th parameters, grp groups
//all others held fixed
double HurdleLikelihood::LL(const vec& th, int grp, bool penalize){
  bool needUpdate = true;
  double pen = 0;
  if(grp > -2) needUpdate = populatePar(grp, th);
  if(k0<=0) return(99999999.0);
  if(needUpdate)   updateGroupSums();
  double nloglik = sum(yI % gpart + y % hpart - cumulant) -.5*Sy2*k0;
  if(penalize) pen = sum(lambda % pengrp);
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
vec HurdleLikelihood::grad(const vec& th, int grp, bool penalize, bool updatePar){
  bool needUpdate = false;
  if(updatePar){
    needUpdate = populatePar(grp, th);
  }
  if(needUpdate){
    updateGroupSums();
  }
  if(grp==-1){//intercepts
    Sdc3(G0) = SyI - sum(cumulant2);//gbb
    dc = (hpart/k0) % cumulant2; //hbb
    Sdc3(H0) = Sy - sum(dc); //yb * dh/dhbb + dc
    dc = -(square(hpart)+k0)/(std::pow(k0,2)*2) % cumulant2;
    Sdc3(K0) = -.5*Sy2 - sum(dc);
    return(-Sdc3/Sn);
  } else{
    dc = ( 2*xI.col(grp) + 0 ) % cumulant2; //gba
    Sdc4(GBA) = 2*SyIxI(grp)-sum(dc);
    dc = (x.col(grp) + 0) % cumulant2;  //hba
    Sdc4(HBA) = SyIx(grp) - sum(dc);
    dc = (0 + xI.col(grp)/k0 % hpart) % cumulant2; //hab
    Sdc4(HAB) = SyxI(grp) - sum(dc);
    dc = -(x.col(grp)/k0 % hpart) % cumulant2; //kba
    Sdc4(KBA) = -Syx(grp) - sum(dc);
    if(penalize){
      ASSERT_TRUE(any(abs(th)>0) || lambda(grp)<datum::eps)
      pen4 = th / pengrp(grp) * lambda(grp);
      return(-Sdc4/Sn + pen4);
    }
    return(-Sdc4/Sn);
  }
}

vec HurdleLikelihood::grad(const vec& th, bool penalize){
  populatePar(th);
  updateGroupSums();
  arma::vec gvec(th.n_elem, fill::zeros);
  arma::uvec4 thIndex(grpm4); //current index into theta
  // we should be taking subvectors of th, but it doesn't matter here
  gvec.elem(grpm3) = grad(th, -1, penalize, false);
  for(uint u=0;u<k; u++){
    gvec.elem(thIndex)=grad(th.elem(thIndex), u, penalize, false);
    thIndex++;
  }
  return gvec;
}

void HurdleLikelihood::setLambda(const vec& lambda_){
  lambda = lambda_;
}

arma::mat HurdleLikelihood::hessian(const arma::vec& th, int grp){
  arma::mat22 a(fill::zeros);
  return(a);
}

// map from parameters into starting indices of theta
arma::uvec::fixed<7> HurdleLikelihood::parmap(int k_) {
  arma::uvec::fixed<7> pm;
  pm[tG0] = 0;
  pm[tGBA] = pm[tG0]+1;
  pm[tH0] = pm[tGBA]+k_;
  pm[tHBA] = pm[tH0]+1;
  pm[tHAB] = pm[tHBA]+k_;
  pm[tKBA] = pm[tHAB]+k_;
  pm[tK0] = pm[tKBA]+k_;
  return(pm);
}


/* 
 * R INTERFACE CODE
 */
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
RcppExport SEXP HurdleLikelihood__LL(SEXP xp, SEXP thR_, SEXP grpR_, SEXP penalize_) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  // convert the parameter to int
  arma::vec th = as<arma::vec>(thR_);
  int grp = as<int>(grpR_);
  bool penalize = as<bool>(penalize_);
  double res = ptr->LL( th , grp, penalize);
  // return the result to R
  return wrap(res);
}

// invoke the loglikelihood (don't curry)
RcppExport SEXP HurdleLikelihood__LLall(SEXP xp, SEXP thR_, SEXP penalize_) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  arma::vec th = as<arma::vec>(thR_);
  bool penalize = as<bool>(penalize_);
  double res = ptr->LL( th , penalize);
  return wrap(res);
}


/// invoke the gradient method (piecewise)
RcppExport SEXP HurdleLikelihood__grad(SEXP xp, SEXP th_, SEXP grp_, SEXP penalize_) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  arma::vec th = as<arma::vec>(th_);
  int grp = as<int>(grp_);
  bool penalize = as<bool>(penalize_);
  arma::vec res = ptr->grad(th, grp, penalize, true);
  return wrap(res);
}

//grad (all coordinates)
RcppExport SEXP HurdleLikelihood__gradAll(SEXP xp, SEXP th_, SEXP penalize_){
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  arma::vec th = as<arma::vec>(th_);
  bool penalize = as<bool>(penalize_);
  arma::vec res = ptr->grad(th, penalize);
  return wrap(res);
}

RcppExport void HurdleLikelihood__setLambda(SEXP xp, SEXP lambda_){
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  arma::vec lambda = as<arma::vec>(lambda_);
  ptr->setLambda(lambda);
}

//debugging code follows
RcppExport SEXP HurdleLikelihood__gpart(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->gpart);
}

RcppExport SEXP HurdleLikelihood__cumulant(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->cumulant);
}

RcppExport SEXP HurdleLikelihood__gplusc(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->gplusc);
}

RcppExport SEXP HurdleLikelihood__xI(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->xI);
}

RcppExport SEXP HurdleLikelihood__x(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->x);
}

RcppExport SEXP HurdleLikelihood__gba(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->gba);
}

RcppExport SEXP HurdleLikelihood__hba(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->hba);
}

RcppExport SEXP HurdleLikelihood__hab(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->hba);
}

RcppExport SEXP HurdleLikelihood__kba(SEXP xp) {
  Rcpp::XPtr<HurdleLikelihood> ptr(xp);
  return wrap(ptr->kba);
}

