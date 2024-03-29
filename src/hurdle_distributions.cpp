#include <Rcpp.h>
#include "assert.h"
#include "hurdle.h"
//#define DEBUG 0
using namespace Rcpp;

double expit(double x){
  return exp(x)/(1+exp(x));
}

// [[Rcpp::export]]
double cpp_rCondHurdle(const NumericVector& x, int j, const NumericMatrix& G,  const NumericMatrix& H,  const NumericMatrix& K, double tol){
  int p=G.nrow();
  NumericVector xI = ifelse(abs(x)>tol, 1.0, 0.0);
  double Gba = G(j,j), Hba=H(j,j), Kba=K(j,j);
  int ii;
  for(int i=0;i<x.length(); i++){ //index into x
    //which entry in the parameters is x?
    ii = (j+i+1) % p; 
    DEXEC(char c[100]; snprintf(c, 100, "K(%d,%d)=%f, x(%d)=%f \n", j, ii, K(ii,j), i, x[i]));
    DPRINT(c);
    Gba = Gba + (2*G(j,ii)*xI[i]+ x[i]*H(j,ii));
    Hba = Hba + (H(ii,j)*xI[i] - K(ii,j)*x[i]);
    DEXEC(snprintf(c, 100, "Hba=%f\n", Hba));
    DPRINT(c);
    //std::cout << "j=" << j << " i=" << i << ", Gba=" << Gba;
  }
  ASSERT((j-ii+p) % p, 1)
  double logitP = Gba - .5*log(Kba/(2*M_PI))+pow(Hba,2)/(2*Kba);
  //std::cout << ". logitP=" << logitP << "\n";
  double mu = Hba/Kba;
  DEXEC(Rprintf("x=(%f,%f), j=%d, mu=%f\n", x[0], x[1], j, mu));
  double y = 0;
  double yI = Rf_runif(0.0, 1.0);//runif(1)[0];
  DEXEC(Rcout << ". expit(logitP)=" << expit(logitP)  << "\n");
  if(yI<expit(logitP)) y = Rf_rnorm(mu, pow(Kba, -.5));//rnorm(1, mu, pow(Kba, -.5))[0];
  return y;	
}

//' @export
// [[Rcpp::export]]
NumericMatrix cpp_rGibbsHurdle(const NumericMatrix& G, const NumericMatrix& H, const NumericMatrix& K, int Nt, double tol){
  int p = G.nrow(), j=0;
  //although samp is a vector, really think of it as a p * Nt matrix
  NumericVector samp(p*Nt);
  std::fill(samp.begin(), samp.end(), 0);
  //we initialized with zeros, now condition on first row
  for(int i=p;i<samp.length(); i++){
    j=i % p;
    NumericVector tmp = samp[seq(i-p+1, i-1)];
    // for( int k=0; k<tmp.length(); k++){
    //   Rprintf("%f, ", tmp[k]);
    // }
    // Rprintf("\n");
    samp[i] = cpp_rCondHurdle(tmp, j, G, H, K, tol);
  }
  return NumericMatrix(p, Nt, samp.begin());
}

