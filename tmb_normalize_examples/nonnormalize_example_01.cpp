// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;


// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
  SparseMatrix<Type> Q;
  Type kappa2 = exp(2. * logkappa);
  Type kappa4 = kappa2*kappa2;
  Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  // INPUTS
  DATA_VECTOR(y);          // obs successes per binomial experiment
  DATA_VECTOR(n);          // trials per experiment

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M1);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M2);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(Aproj); // used to project from spatial mesh to data locations

  // Prior means
  DATA_SCALAR(alphamean); // mean of prior for alpha
  DATA_SCALAR(logtaumean); // mean of prior for logtau
  DATA_SCALAR(logkappamean); // mean of prior for logkappa

  // Parameters
  PARAMETER(alpha);   // intercept
  PARAMETER(logtau);           // log of INLA tau param (precision of spatial covariance mat)
  PARAMETER(logkappa);         // log of INLA kappa - related to spatial correlation and range

  PARAMETER_VECTOR(S);    // spatial random effects


  // Make transformations of some of our parameters
  int num = y.size();
  vector<Type> projS(num);                // value of gmrf at data points

  SparseMatrix<Type> Q   = spde_Q(logkappa, logtau, M0, M1, M2);

  Type nll = 0;

  // Hyperpriors
  nll -= dnorm(logtau, logtaumean , Type(1.0), true);
  nll -= dnorm(logkappa, logkappamean , Type(1.0), true);

  
  // Spatial random effects
  nll += GMRF(Q)(S);

  // Fixed effects
  nll -= dnorm(alpha, alphamean, Type(1.0), true);

  // Likelihood
  projS = Aproj * S; // Project S at mesh pts to data pts

  for (int i = 0; i < num; i++) {
    nll -= dbinom_robust( y(i), n(i), alpha + projS(i), true );
  }
  
  return nll;
}