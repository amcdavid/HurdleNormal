//#define DEBUG 0
#include "hurdle_likelihood.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
HurdleLikelihood::HurdleLikelihood (const arma::vec& y_, const arma::mat& x_, const arma::ivec& grp_, const arma::vec& th, const arma::vec& lambda_, double tol)
  : y(y_), yI(y_.n_rows, fill::ones),
    x(x_), xI(x_.n_rows, x_.n_cols,fill::ones),
    grp(grp_), g0(0), h0(0), k0(1),
    gba(x_.n_cols, fill::zeros), hba(x_.n_cols, fill::zeros), hab(x_.n_cols, fill::zeros), kba(x_.n_cols, fill::zeros),
    gpart(x_.n_rows, fill::zeros), hpart(x_.n_rows, fill::zeros),
    gplusc(x_.n_rows, fill::zeros), cumulant(x_.n_rows, fill::zeros), cumulant2(x_.n_rows, fill::zeros),
    Sy2(sum(square(y_))), Sy(sum(y_)), SyI(999),
    Sn(x_.n_rows),k(x_.n_cols),
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
  // populate params
  for(int i=-1; i<k; i++){
    populatePar(i, th.elem(find(grp==i)));
  }
    gpart = g0+2*xI*gba + x*hba;
    hpart = h0+xI*hab -x*kba;
    DPRINT("Finished construction. dim(gpart)=" << gpart.n_rows << "gpart(0)=" << gpart(0)<<"\n");
}

  // pack parameters group `grp` with th
//do online updates of crossproducts, gpart and hpart 
  bool HurdleLikelihood::populatePar(int grp, const vec& th){
    //    DPRINT("grp=" << grp << " th=" << th << "\n");
    bool changed=true;
    //update param and cross products
    if(grp == -1){ //intercept
      changed = std::abs(g0-th(0))> datum::eps || 
      	std::abs(h0-th(1))>datum::eps || std::abs(k0-th(2))>datum::eps;
      if(changed){
	//	DPRINT("gpart(0)=" << gpart(0));
	gpart += (th(0)-g0);
	hpart += (th(1)-h0);
	g0 = th(0), h0=th(1), k0=th(2);
	//	DPRINT("gpart(0)=" << gpart(0));
      }
    } else{ //covariates
      changed = std::abs(gba(grp)-th(0))> 0 || 
      std::abs(hba(grp)-th(1))>0 || 
      std::abs(hab(grp)-th(2))>0 ||
      std::abs(kba(grp)-th(3))>0;
	if(changed){
	  DPRINT("gpart(0)=" << gpart(0));
	  gpart += (2*xI.col(grp)*(th(0)-gba(grp)) + 
		    x.col(grp)*(th(1)-hba(grp)));
	  hpart +=  (xI.col(grp)*(th(2)-hab(grp)) - 
		     x.col(grp)*(th(3)-kba(grp)));
	  // //	  DPRINT("gpart(0)=" << gpart(0));
	  gba(grp)=th(0), hba(grp)=th(1), hab(grp)=th(2), kba(grp)=th(3);
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



// void HurdleLikelihood::updateCrossproducts(int i){
//   xI_gba.col(i) = xI.col(i)*gba(i);
//   x_hba.col(i) = x.col(i)*hba(i);
//   xI_hab.col(i) = xI.col(i)*hab(i);
//   x_kba.col(i) = x.col(i)*kba(i);
//   ASSERT(x_kba.n_cols, x.n_cols)
// }


void HurdleLikelihood::updateGroupSums(bool updateSums, bool updateGrad){
  if(updateSums){
    gplusc = gpart + -.5*log(k0/(2*datum::pi))+pow(hpart, 2)/(2*k0);
    // log(1 + e^x) = x unless x < large 
    // for(int i=0; i<Sn; i++){
    //   if(gplusc[i]<large){
    // 	cumulant[i] = std::log(1+std::exp(gplusc[i]));
    // 	cumulant2[i] = 1/(1+std::exp(-gplusc[i]));
    // 	  } else{
    // 	cumulant[i] = gplusc[i];
    // 	cumulant2[i]=1;
    //   }
    // }
    uvec small = find(gplusc<large);
    cumulant = gplusc;
    cumulant.elem(small) = log(1+exp(cumulant.elem(small)));
    cumulant2.ones();
    cumulant2.elem(small) = 1/(1+exp(-gplusc.elem(small)));
  }
}

// Replace all coordinates
double HurdleLikelihood::LL(const vec& th){
  //DPRINT("replace all. x.n_cols=" << x.n_cols <<"\n");
  for(int i=-1; i<k; i++){
    //    DPRINT("i=" << i << "\n");
    vec subth = th.elem(find(grp==i));
    populatePar(i, subth);
  }
  return(LL(th, -2));
}

//th parameters, grp groups
//all others held fixed
double HurdleLikelihood::LL(const vec& th, int grp){
  bool needUpdate = true;
  if(grp > -2) needUpdate = populatePar(grp, th);
  if(k0<=0) return(99999999);
  updateGroupSums(needUpdate, false);
  double nloglik = sum(yI % gpart + y % hpart - cumulant) -.5*Sy2*k0;
  double pen = sum(lambda % pengrp);
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
vec HurdleLikelihood::grad(const vec& th, int grp, bool penalize){
  vec dc(Sn, fill::zeros); //temp storage for derivatives in front of cumulant
  vec Sdc(th.n_elem, fill::zeros);
  vec pen(th.n_elem, fill::zeros);

  bool needUpdate = populatePar(grp, th);
  updateGroupSums(needUpdate, true);
  if(grp==-1){//intercepts
    //dc.col(0)=1 - cumulant2;
    Sdc(0) = SyI - sum(cumulant2);//gbb
    dc = (hpart/k0) % cumulant2; //hbb
    Sdc(1) = Sy - sum(dc); //yb * dh/dhbb + dc
    dc = -(square(hpart)+k0)/(std::pow(k0,2)*2) % cumulant2;
    Sdc(2) = -.5*Sy2 - sum(dc);
  } else{
    dc = ( 2*xI.col(grp) + 0 ) % cumulant2; //gba
    Sdc(0) = 2*SyIxI(grp)-sum(dc);
    dc = (x.col(grp) + 0) % cumulant2;  //hba
    Sdc(1) = SyIx(grp) - sum(dc);
    dc = (0 + xI.col(grp)/k0 % hpart) % cumulant2; //hab
    Sdc(2) = SyxI(grp) - sum(dc);
    dc = -(x.col(grp)/k0 % hpart) % cumulant2; //kba
    Sdc(3) = -Syx(grp) - sum(dc);
    if(penalize){
      ASSERT_TRUE(any(abs(th)>0) || lambda(grp)<datum::eps)
      pen = th / pengrp(grp) * lambda(grp);
    }
  }
  //handle penalty
  return(-Sdc/Sn + pen);
}
