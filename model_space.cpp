// /////////////////////////////////////////////////////////////////////////////
// RB AOZ
// 2018
// Template file for fitting space only models
// /////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// 3. ref https://github.com/nmmarquez/re_simulations/blob/master/inla/sta.cpp
//        https://github.com/nmmarquez/re_simulations/blob/master/inla/SPDEAR1AR1.R
// /////////////////////////////////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
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

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}



template<class Type>
Type objective_function<Type>::operator() ()
{

  // ~~~~~~~~~------------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------------~~
  // FIRST, we define params/values/data that will be passed in from R
  // ~~~~~~~~~~~------------------------------------------------------
  // ~~~~~~~~~------------------------------------------------------~~

  // normalization flag
  DATA_INTEGER( flag ); // flag == 1 => no data contribution added to jnll
  
  // Indices
  DATA_INTEGER( num_i );   // Number of data points in space
  DATA_INTEGER( num_s );   // Number of mesh points in space mesh

  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_i );   // Num occurrences (deaths) per binomial experiment at point i (cluster)
  DATA_VECTOR( n_i );   // Trials per cluster
  DATA_MATRIX( X_alpha );  // Covariate 'design matrix' for just intercept (i.e. 1col matrix of all 1s)
  DATA_MATRIX( X_betas );   // Covariate design matrix excluding intercept columm

  // SPDE objects
  DATA_SPARSE_MATRIX( M0 );
  DATA_SPARSE_MATRIX( M1 );
  DATA_SPARSE_MATRIX( M2 );
  DATA_SPARSE_MATRIX( Aproj );   // Used to project spatial mesh to data locations

  // Options
  DATA_VECTOR( options );
  // options[0] == 1 : turn on adreport
  // options[1] == 1 : use priors
  // options[2] == 1 : fit with intercept
  // options[3] == 1 : fit with cov effects
  // options[4] == 1 : fit with nugget
  // options[5] == 0 : use normal data lik
  // options[5] == 1 : use binom  data lik 

  // Fixed effects
  PARAMETER( alpha );            // Intercept
  PARAMETER_VECTOR( betas );     // Covariate coefficients
  PARAMETER( log_obs_sigma);     // if using normal likelihood, sd of single normal obs
  PARAMETER( log_tau );          // Log of INLA tau param (precision of space covariance matrix)
  PARAMETER( log_kappa );        // Log of INLA kappa (related to spatial correlation and range)
  PARAMETER( log_nugget_sigma ); // Log of SD for irreducible nugget variance
  PARAMETER_VECTOR( nug_i );     // random effect nugget for each observation
  
  // Random effects
  PARAMETER_ARRAY( Epsilon_s );  // Random effect for each spatial mesh location. Currently a 1d array of num_s

  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~
  
  // objective function -- joint negative log-likelihood
  Type jnll = 0;
  //parallel_accumulator<Type> jnll(this); // parallelize jnll NOTE: seems to break with omp>1 and mkl >=1

  // print parallel info
  // max_parallel_regions = omp_get_max_threads();
  max_parallel_regions = 1;
  // printf("This is thread %d\n", max_parallel_regions);

  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(log_kappa, log_tau, M0, M1, M2);

  // Transform some of our parameters
  Type range = sqrt(8.0) / exp(log_kappa);
  Type sigma = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * log_tau) * exp(2.0 * log_kappa));
  Type nugget_sigma  = exp(log_nugget_sigma);

  // Define objects for derived values
  vector<Type> fe_i(num_i);              // main effect alpha + X_betas %*% t(betas)
  vector<Type> latent_field_i(num_i);      // Logit estimated prob for each point i
  vector<Type> epsilon_s(num_s);         // Epsilon_s (array) unlisted into a vector for easier matrix multiplication
  vector<Type> projepsilon_i(num_i);     // value of gmrf at data points
  
  // evaluate fixed effects for intercept and covs of applicable
  fe_i = X_alpha * Type(0.0); // initialize
  if(options[2] == 1){
    fe_i = X_alpha * alpha; // add on intercept if using
  }
  if(options[3] == 1){
    fe_i = X_betas * betas.matrix(); // add on covariate effects if using
  }
  
  // Transform GMRFs and make vector form
  for(int s = 0; s < num_s; s++){
    epsilon_s[s] = Epsilon_s(s); // This is probably unnecssary in space-only...
  }

  // Project GP approx from mesh points to data points 
  projepsilon_i = Aproj * epsilon_s.matrix();

  // ~~~~~~~~~------------------------------------------------~~-
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) priors
  // 2) GP field
  // 3) data
  // ~~~~~~~~~------------------------------------------------~~-
  // ~~~~~~~~~------------------------------------------------~~-

  ///////// 
  // (1) // 
  /////////
  // Prior contributions to joint likelihood (if options[0]==1)
  if(options[1] == 1) {

    // add in priors for spde gp
    jnll -= dnorm(log_tau,   Type(0.0), Type(1.0), true); // N(0,1) prior for logtau
    jnll -= dnorm(log_kappa, Type(0.0), Type(1.0), true); // N(0,1) prior for logkappa

    // prior for intercept
    if(options[2] == 1){
      jnll -= dnorm(alpha, Type(0.0), Type(3), true); // N(0, 3)
    }

    // prior for covariate coefs
    if(options[3] == 1){
      for( int j = 0; j < betas.size(); j++){
	jnll -= dnorm(betas(j), Type(0.0), Type(3), true); // N(0, 3)
      }
    }

    // prior for log(nugget sd)
    if(options[4] == 1){
      jnll -= dnorm(log_nugget_sigma, Type(-4.0), Type(2.0), true); // N(-4, 2)
    }

    // prior for log(obs sd) of using normal data lik
    if(options[5] == 0){
      jnll -= dnorm(log_obs_sigma, Type(-4.0), Type(2.0), true); // N(-4, 2)
    }
    
  } 

  /////////
  // (2) //
  /////////
  // 'GP' field contribution (i.e. log-lik of Gaussian-Markov random fields, GMRFs)
  // NOTE: likelihoods from namespace 'density' already return NEGATIVE log-liks so we add
  //       other likelihoods return positibe log-liks
  if(flag == 1){
    // then we are not calculating the normalizing constant in the inner opt.
    // that norm constant means taking an expensive determinant of Q_ss
    jnll += GMRF(Q_ss, false)(epsilon_s);
  }else{
    jnll += GMRF(Q_ss)(epsilon_s);
  }

  // nugget contribution to the likelihood
  if(options[4] == 1 ){
    for (int i = 0; i < num_i; i++){
      jnll -= dnorm(nug_i(i), Type(0.0), nugget_sigma, true);
    }
  }

  /////////
  // (3) //
  /////////
  // Likelihood contribution from each datapoint i

  if (flag == 1) return jnll; // return without data ll contrib to avoid unneccesary log(det(Q)) calcs

  for (int i = 0; i < num_i; i++){

    // latent field estimate at each obs
    latent_field_i(i) = fe_i(i) + projepsilon_i(i);
    // add on nugget if using
    if(options[4] == 1){ 
      latent_field_i(i) = latent_field_i(i) + nug_i(i);
    }

    // and add data contribution to jnll
    if(!isNA(y_i(i))){

      // if normal
      if(options[5] ==  0){
	// Uses the dbinom_robust function, which takes the logit probability
	jnll -= dnorm( y_i(i), latent_field_i(i), exp(log_obs_sigma)/n_i(i), true );
      }

      // if binom
      if(options[5] ==  1){
	// Uses the dbinom_robust function, which takes the logit probability
	jnll -= dbinom_robust( y_i(i), n_i(i), latent_field_i(i), true );
      }
      
    } // !isNA

  } // for( i )
  

  // ~~~~~~~~~~~
  // ADREPORT
  // ~~~~~~~~~~~
  if(options[0] == 1){
    ADREPORT(betas);
    ADREPORT(Epsilon_s);
  }

  return jnll;
}




