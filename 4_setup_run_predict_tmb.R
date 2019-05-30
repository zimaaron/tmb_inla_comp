## ######
## ######
## TMB ##
## ######
## ######

## ########
## SETUP ##
## ########

## build design mats for int and covs
## this is done seperately to allow indep. turning each on/off
X_alpha <- matrix(rep(1, nrow(dt), ncol = 1))
X_betas <- as.matrix(dt[, covs[, name], with=FALSE])

templ <- "model_space"
setwd("/homes/azimmer/tmb_inla_comp/")
TMB::compile(paste0('./', templ,".cpp"))
dyn.load( dynlib(templ) )


## function to convert from data lik string to integer
## allows easily adding more options even though overkill for just 2
tmb.lik.dict <- function(x){
  dict <- list(normal = 0,
               binom = 1)
  dict[[x]]
}

## setup data to feed into the model
data_full <- list(num_i = nrow(dt),  # Total number of observations
                  num_s = mesh_s$n,  # Number of vertices in SPDE mesh
                  y_i   = dt[,Y],    # Number of observed deaths in the cluster
                  n_i   = dt[,N],    # Number of exposures in the cluster
                  X_alpha  = X_alpha,# Covariate design matrix
                  X_betas  = X_betas,# Covariate design matrix
                  M0    = spde$param.inla$M0, # SPDE sparse matrix
                  M1    = spde$param.inla$M1, # SPDE sparse matrix
                  M2    = spde$param.inla$M2, # SPDE sparse matrix
                  Aproj = A.proj,             # Projection matrix
                  options = c(1, ## if 1, run adreport 
                              1, ## if 1, use priors
                              ifelse(is.null(alpha), 0, 1), ## if 1, run with intercept
                              ifelse(is.null(betas), 0, 1), ## if 1, run with covs
                              ifelse(is.null(nug.var), 0, 1), ## if 1, run with nugget
                              tmb.lik.dict(data.lik), ## if 0, normal data. if 1, binom data lik
                              1 ## use normalization trick?
                              ),
                  flag = 1, # normalization flag. when 0, prior is returned. when 1 data is included in jnll
                  norm_prec_pri = norm.prec.pri, ## gamma on log(prec)
                  nug_prec_pri = nug.prec.pri, ## gamma on log(prec)
                  alphaj_pri = alphaj.pri, ## normal
                  logtau_pri = spde.theta1.pri, ## normal logtau
                  logkappa_pri = spde.theta1.pri ## normal logkappa
                  )

## Specify starting values for TMB parameters for GP
tmb_params <- list(alpha = 0.0, # intercept
                   betas = rep(0, ncol(X_betas)), # cov effects
                   log_obs_sigma = 0.0, # log(data sd) if using normal dist
                   log_tau   = 1.0, # Log inverse of tau (Epsilon)
                   log_kappa = 1.0, # Matern range parameter
                   log_nugget_sigma = -1.0, # log of nugget sd
                   nug_i = rep(0, nrow(dt)), # vector of nugget random effects
                   Epsilon_s = matrix(1, nrow=nodes, ncol=1) # GP value at obs locs
                   )

## make a list of things that are random effects
rand_effs <- c('Epsilon_s')

## NULL out things that aren't in this run and identify extra rand effs if using them
ADmap <- list()
if(is.null(alpha)){
  ADmap[['alpha']] <- factor(NA)
}
if(is.null(betas)){
  ADmap[['betas']] <- rep(factor(NA), ncol(X_betas))
}
if(is.null(nug.var)){
  ADmap[['log_nugget_sigma']] <- factor(NA)
  ADmap[['nug_i']] <- rep(factor(NA), nrow(dt))
}else{
  rand_effs <- c(rand_effs, 'nug_i')
}
if(data.lik == 'binom'){
  ADmap[['log_obs_sigma']] <- factor(NA)
}

## make the autodiff generated liklihood func & gradient
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 map = ADmap, 
                 hessian=TRUE,
                 DLL=templ)

## should we use the normalization flag?
if(data_full$options[7] == 1){
  obj <- normalize(obj, flag="flag", value = 0) ## value: Value of 'flag' that signifies to not include the data term.
}

## Run optimizer
ptm <- proc.time()[3]
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr,
                                        #                                lower = rep(-10, ncol(X_xp) + 3), ## TODO
                                        #                                upper = rep( 10, ncol(X_xp) + 3), ## TODO
                              control     =    list(trace=1)))
fit_time_tmb <- proc.time()[3] - ptm

## Get standard errors
SD0 = TMB::sdreport(obj, getJointPrecision=TRUE,
                    bias.correct = bias.correct, 
                    bias.correct.control = list(sd = sd.correct)) 
tmb_total_fit_time <- proc.time()[3] - ptm 
tmb_sdreport_time <-  tmb_total_fit_time - fit_time_tmb


## make a flag to see if tmb converged with a PD cov structure
## we check in a minute to see if it actually is true
tmb.pd.converge <- TRUE

## ##########
## PREDICT ##
## ##########

message('making predictions')

## now we can take draws and project to space-time raster locs
mu <- c(SD0$par.fixed,SD0$par.random)

## simulate draws
ptm2 <- proc.time()[3]
rmvnorm_prec <- function(mu, chol_prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- chol_prec #Cholesky(prec, super=TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}

L <- try(suppressWarnings(Cholesky(SD0$jointPrecision, super = T)), silent = TRUE)
if(class(L) == "try-error"){
  tmb.pd.converge <- FALSE ## the jointPrec was not PD
  message('TMB PRECISION IS NOT! PD - mapping to nearest PD precision ')
  message('TMB PRECISION IS NOT! PD - mapping to nearest PD precision ')
  message('TMB PRECISION IS NOT! PD - mapping to nearest PD precision ')
  SD0$jointPrecision <- Matrix(nearPD(SD0$jointPrecision)$mat, sparse = T)
  L <- Cholesky(SD0$jointPrecision, super = T)
}
tmb_draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = ndraws)
tmb_get_draws_time <- proc.time()[3] - ptm2

## separate out the tmb_draws
parnames <- c(names(SD0$par.fixed), names(SD0$par.random))
epsilon_tmb_draws  <- tmb_draws[parnames == 'Epsilon_s',]
alpha_tmb_draws    <- matrix(tmb_draws[parnames == 'alpha',], nrow = 1)
betas_tmb_draws    <- tmb_draws[parnames == 'betas',]
if(!is.matrix(betas_tmb_draws)) betas_tmb_draws <- matrix(betas_tmb_draws, nrow = 1)
log_kappa_tmb_draws <- tmb_draws[parnames == 'log_kappa',]
log_tau_tmb_draws  <- tmb_draws[parnames == 'log_tau',]
log_nugget_sigma_draws <- tmb_draws[parnames == 'log_nugget_sigma', ]
log_gauss_sigma_draws <- tmb_draws[parnames == 'log_obs_sigma', ]

## values of S at each cell (long by nperiods)
## rows: pixels, cols: posterior draws
pred_tmb <- as.matrix(A.pred %*% epsilon_tmb_draws)

if(!is.null(alpha)){
  ## add on intercept, one alpha draw per row
  pred_tmb <- sweep(pred_tmb, 2, alpha_tmb_draws, '+')
}

if(!is.null(betas)){
  ## add on covariate values by draw
  tmb_vals <- list()
  for(p in 1:nperiods) tmb_vals[[p]] <- cov_vals[[p]] %*% betas_tmb_draws

  cell_b <- do.call(rbind, tmb_vals)

  ## add together linear and st components
  pred_tmb <- cell_b + pred_tmb
}

if(!is.null(nug.var)){
  ## simulate nugget noise
  cell_nug <- do.call(cbind, lapply(1:ndraws,
                                    FUN = function(x){rnorm(n = nrow(pred_tmb),
                                                            mean = 0,
                                                            sd = exp(log_nugget_sigma_draws)[x])})
                      )
  ## add it on
  pred_tmb <- cell_nug + pred_tmb
}

## save prediction timing
totalpredict_time_tmb <- proc.time()[3] - ptm2

## #######
## SAVE ##
## #######

## summarize the latent field
summ_tmb <- cbind(median = (apply(pred_tmb, 1, median)),
                  sd     = (apply(pred_tmb, 1, sd)))

ras_med_tmb <- insertRaster(simple_raster, matrix(summ_tmb[, 1], ncol = nperiods))
ras_sdv_tmb <- insertRaster(simple_raster, matrix(summ_tmb[, 2], ncol = nperiods))

saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%04d_iter%04d_tmb_preds_median_raster.rds', out.dir, par.iter, iii), object = ras_med_tmb)
saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%04d_iter%04d_tmb_preds_stdev_raster.rds', out.dir, par.iter, iii), object = ras_sdv_tmb)

if(data.lik == 'binom'){
  ## convert to prevalence space and summarize, rasterize, and save again
  pred_tmb_p <- plogis(pred_tmb)
  
  summ_tmb_p <- cbind(median = (apply(pred_tmb_p, 1, median)),
                      sd     = (apply(pred_tmb_p, 1, sd)))

  ras_med_tmb_p <- insertRaster(simple_raster, matrix(summ_tmb_p[, 1], ncol = nperiods))
  ras_sdv_tmb_p <- insertRaster(simple_raster, matrix(summ_tmb_p[, 2], ncol = nperiods))

  saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%04d_iter%04d_tmb_preds_median_raster_PREV.rds', out.dir, par.iter, iii),
          object = ras_med_tmb_p)
  saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%04d_iter%04d_tmb_preds_stdev_raster_PREV.rds', out.dir, par.iter, iii),
          object = ras_sdv_tmb_p)
}