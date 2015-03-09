/* Declarations for Meier and Buhlmann's algorithm for quadratic
 * approximate coordinate descent for group-L1 penalized likelihoods
 *
 *
 */

#include <RcppArmadillo.h>

class PenalizedLikelihood
{
public:
  virtual double LL(const arma::vec& th, int grp, bool penalize) = 0;
  virtual arma::vec grad(const arma::vec& th, int grp, bool penalize, bool updatePar) = 0;
  virtual arma::mat hessian(const arma::vec& th, int grp) = 0;
  virtual void setLambda(const arma::vec& lambda_) = 0;
};

class GrpLassoPath
{
private:
double slambda;
double lineSearch(arma::vec direction, int grp);
arma::vec descentDirection(arma::vec point, int grp);
std::list<int> activeSet; //invariant: union(activeSet, inactiveSet) = set
std::list<int> inactiveSet;
unsigned maxit; //max touches of groups
double ctol;    //kkt convergence tolerance
//armijo parameters?

public:
PenalizedLikelihood* pl;  

};
