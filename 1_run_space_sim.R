## this script simulates some realistic datasets for comparison between INLA and TMB
## it leverages existing architecture that the LBD team at IHME has already created
## written by AOZ
## last editted Oct 3, 2018

## options(error = recover)

## ## for working on local laptop
## load('/homes/azimmer/scratch/tmb_space_debug.RData')
## load('~/Desktop/tmb_inla_comp/scratch/tmb_space_debug.RData')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par.iter <- 1  ## as.numeric(  commandArgs()[4]) ## all we need is to grab the (parallel) iteration of this run
run_date <- "2018_06_12_15_54_38" ## as.character(commandArgs()[5]) ## and the run_date so we know where to load from
run_date <- "2018_10_03_12_51_58"
run_date <- "2018_12_01_16_28_20"
run_date <- "2019_02_10_12_59_31"

#############################################
## setup the environment for singularity R ##
#############################################

## Set core_repo location and tmb_repo loc
user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)
pull_tmb_git <- FALSE

## grab libraries and functions from MBG code
setwd(core_repo)
commondir    <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, 'mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(scales)
library(RandomFields)

## Now we can switch to the TMB repo
setwd(tmb_repo)
if(pull_tmb_git) system(sprintf('cd %s\ngit pull %s %s', core_repo, remote, branch))
source('./realistic_sim_utils.R')

## setup is now done. setup some parameters for this simulation

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################################################################
## load in the loopvars from launch and setup all the params for this job ##
############################################################################
out.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', run_date)
loopvars <- readRDS(file = paste0(out.dir, '/loopvars.rds'))

reg             <- as.character(loopvars[par.iter, 1])
year_list       <- eval(parse(text = loopvars[par.iter, 2]))
cov_names       <- eval(parse(text = as.character(loopvars[par.iter, 3])))
cov_measures    <- eval(parse(text = as.character(loopvars[par.iter, 4])))
betas           <- eval(parse(text = as.character(loopvars[par.iter, 5])));if(is.na(betas)) betas <- NULL

alpha           <- as.numeric(loopvars[par.iter, 6]);if(is.na(alpha)) alpha <- NULL
sp.range        <- as.numeric(loopvars[par.iter, 7])
sp.var          <- as.numeric(loopvars[par.iter, 8])
sp.alpha        <- as.numeric(loopvars[par.iter, 9])
nug.var         <- as.numeric(loopvars[par.iter, 10]);if(is.na(nug.var)) nug.var <- NULL

t.rho           <- as.numeric(loopvars[par.iter, 11])
mesh_s_max_edge <- as.character(loopvars[par.iter, 12])
n.clust         <- as.numeric(loopvars[par.iter, 13])
m.clust         <- as.numeric(loopvars[par.iter, 14])
sample.strat    <- eval(parse(text = as.character(loopvars[par.iter, 15])))
obs.loc.strat   <- sample.strat[['obs.loc.strat']]
urban.pop.pct   <- sample.strat[['urban.pop.pct']]
urban.strat.pct <- sample.strat[['urban.strat.pct']]

cores           <- as.numeric(loopvars[par.iter, 16]) 
ndraws          <- as.numeric(loopvars[par.iter, 17])
alphaj.pri      <- eval(parse(text = as.character(loopvars[par.iter, 18]))) ## normal mean and sd ## TODO pass this to INLA and TMB
nug.pri         <- eval(parse(text = as.character(loopvars[par.iter, 19])))  ## gamma for nug preciion with shape and inv-scale ## TODO pass this to INLA and TMB
inla.int.strat  <- as.character(loopvars[par.iter, 20]) ## can be one of: 'eb', 'ccd', 'grid'

inla.approx     <- as.character(loopvars[par.iter, 21]) ## can be one of: 'eb', 'ccd', 'grid'
l.tau.pri       <- NULL  ## taken from INLA spde mesh obj
l.kap.pri       <- NULL  ## taken from INLA spde mesh obj
Nsim <-  as.numeric(loopvars[par.iter, 22]) ## number of times to repeat simulation
data.lik <- as.character(loopvars[par.iter, 23])
sd.norm <-  as.numeric(loopvars[par.iter, 24])

## TODO? add in some validation options? or maybe just always do them all


## from these imputs, make a table of covariate names and measures
covs <- data.table(name = cov_names, meas = cov_measures)

## I hardcode a few other options that are useful sometimes when running interactively
## these can probably be deleted...

## to make it easier to run real data from this code, usually have SIM=TRUE, REAL=FALSE
use_real_data <- FALSE
use_sim_data  <- TRUE

## should we save this data to mbg input_data? useful if we want to run sim data through lbd_core
save.as.input.data <- FALSE
data.tag <- '_allyrs_nug'

## end of user inputs
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## transform some inputs into other useful quantities

## convert sp.range and sp.var into sp.kappa for rspde function
sp.kappa   <- sqrt(8) / sp.range
logkappa   <- log(sp.kappa)
sp.tau     <- sqrt(1 / (4 * pi * sp.kappa ^ 2 * sp.var)) ## sp.var = 1/(4*pi*kappa^2*tau^2)
logtau     <- log(sp.tau)
trho_trans <- log((-1 - t.rho) / (t.rho - 1))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################
## load in region/counry shapes and covs ##
###########################################
dir.create(sprintf('%s/simulated_obj', out.dir), recursive = TRUE)

## load in the region shapefile and prep the boundary
gaul_list           <- get_adm0_codes(reg)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']]

###################
## simulate data ##
###################

##TODO fix up truths now that things can be shut off... (betas, alpha, covs, ...)

## make an object with true param values
true.param.names <- true.param.vals <- c()
if(!is.null(alpha){
  
}

true.params <- data.table(param = c('int',
                                    cov_names,
                                    'nom. range',
                                    'nom. var',
                                    'time rho'
                                    ),
                          truth = c(alpha,
                                    betas,
                                    sp.range,
                                    sp.var,
                                    t.rho))

## true.par.vec <- c(alpha, betas, logtau, logkappa, trho_trans)
## names(true.par.vec) <- c(rep('alpha_j', length(betas) + 1), 'logtau', 'logkappa', 'trho_trans')
## if(length(year_list) == 1) true.par.vec <- true.par.vec[-length(true.par.vec)]

saveRDS(file = sprintf('%s/simulated_obj/true_param_table.rds', out.dir),
        object = true.params)


###########
###########
###########

for(iii in 1:Nsim){ ## repeat Nsim times

  ## TODO get this logic working to speed things up...
  ## if(iii == 1){ ## first time, must load covs, after that, we can reuse them
  
  sim.obj <- sim.realistic.data(reg = reg,
                                  year_list = year_list,
                                  data.lik = data.lik,
                                  sd.norm = sd.norm, 
                                  betas = betas,
                                  sp.kappa = sp.kappa,
                                  sp.alpha = sp.alpha,
                                  t.rho = t.rho,
                                  nug.var = nug.var, 
                                  n.clust = n.clust,
                                  m.clust = m.clust,
                                  covs = covs,
                                  cov_layers = NULL, 
                                  simple_raster = simple_raster,
                                  simple_polygon = simple_polygon,
                                  pop_raster = pop_raster, 
                                  obs.loc.strat = obs.loc.strat,
                                  urban.pop.pct = urban.pop.pct,
                                  urban.strat.pct = urban.strat.pct, 
                                  out.dir = paste(out.dir, iii, sep = '/'),
                                  sp.field.sim.strat = 'SPDE', 
                                  seed = NULL)

  ## }else{
  ##    sim.obj <- sim.realistic.data(reg = reg,
  ##                                 year_list = year_list,
  ##                                 betas = betas,
  ##                                 sp.kappa = sp.kappa,
  ##                                 sp.alpha = sp.alpha,
  ##                                 t.rho = t.rho,
  ##                                 nug.var = nug.var, 
  ##                                 n.clust = n.clust,
  ##                                 m.clust = m.clust,
  ##                                 covs = covs,
  ##                                 cov_layers = cov_list, ## which is created from each sim.obj after stripping GP from the list. ~line228
  ##                                 simple_raster = simple_raster,
  ##                                 simple_polygon = simple_polygon,
  ##                                 pop_raster = pop_raster, 
  ##                                 obs.loc.strat = obs.loc.strat,
  ##                                 urban.pop.pct = urban.pop.pct,
  ##                                 urban.strat.pct = urban.strat.pct, 
  ##                                 out.dir = paste(out.dir, iii, sep = '/'),
  ##                                 sp.field.sim.strat = 'SPDE', 
  ##                                 seed = NULL)
  ## }

  saveRDS(file = sprintf('%s/simulated_obj/sim_obj_%i.rds', out.dir, iii),
          object = sim.obj)

  dt <- sim.obj$sim.dat ## simulated data, lat-long, year, covs, true surface
  covs.gp <- sim.obj$cov.gp.rasters   ## rasters of covs and true simulated gp field
  true.gp <- covs.gp[['gp']]
  cov_list <- covs.gp[!grepl('gp',names(covs.gp))]
  true.rast <- sim.obj$true.rast

  ## save (if desired) this simulated dataset to .../mbg/input_data for mbg pipeline
  if(save.as.input.data){
    df <- data.table(longitude = dt$long,
                     latitude = dt$lat,
                     year = dt$year,
                     country = reg,
                     N = dt$N,
                     simulation = dt$Y, 
                     weight = 1)
    write.csv(file = sprintf('/share/geospatial/mbg/input_data/simulation%s.csv', data.tag),
              x = df)
  }


  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## ##################################
  ## ##################################
  ## setup for tmb and INLA modeling ##
  ## ##################################
  ## ##################################
  dir.create(sprintf('%s/modeling/inputs', out.dir), recursive = TRUE, showWarnings = F)
  dir.create(sprintf('%s/modeling/tmb/outputs', out.dir), recursive = TRUE, showWarnings = F)
  dir.create(sprintf('%s/modeling/inla/outputs', out.dir), recursive = TRUE, showWarnings = F)

  ## ######################################
  ## load in some real data (if desired) ##
  ## ######################################

  if(use_real_data){

    reg <- 'sssa'
    indicator = 'hiv_test'
    indicator_group = 'hiv'
    age <- holdout <- test <- 0
    yearload <- 'annual'
    withtag <- TRUE
    datatag <- '_survey'
    use_share <- 'FALSE'
    year_list <- 2000:2016
    
    pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

    
    ## load in the region shapefile and prep the boundary
    gaul_list           <- get_gaul_codes(reg)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
    subset_shape        <- simple_polygon_list[[1]]
    simple_polygon      <- simple_polygon_list[[2]]

    ## Load list of raster inputs (pop and simple)
    raster_list        <- build_simple_raster_pop(subset_shape)
    simple_raster      <- raster_list[['simple_raster']]
    pop_raster         <- raster_list[['pop_raster']]

    run_date <- make_time_stamp(TRUE)
    dt <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                          simple      = simple_polygon,
                          agebin      = age,
                          removeyemen = TRUE,
                          pathaddin   = pathaddin,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list)

    ## just to make sure everything goes smoothly, add in single datapoints to missing years
    missing.yrs <- setdiff(year_list, unique(dt[, year]))

    if(length(missing.yrs) > 0){
      for(yy in missing.yrs){
        new.row <- dt[1, ]
        new.row[, weight := 0]
        new.row[, year := yy]
        dt <- rbind(dt, new.row)
      }
    }
    

    dt[, Y:= hiv_test]

  }


  ## ############################
  ## SETUP SOME SHARED OBJECTS ##
  ## ############################

  ## ~~~
  ## required for model fit
  ## ~~~

  dt[, id := 1:.N]
  dt[, period_id := as.numeric(as.factor(dt[, year]))]
  dt.coords <- cbind(dt$long, dt$lat)
  dt.pers   <- dt[, period_id]
  nperiods  <- length(year_list)

  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d           = dt,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = "c(1, 5)")

  pdf(sprintf('%s/modeling/inputs/mesh_%i.pdf', out.dir, iii))
  plot(mesh_s)
  plot(simple_raster, add = TRUE) ## just to show loc of simple_raster under mesh for scale
  plot(mesh_s, add = TRUE)
  dev.off()

  nodes <- mesh_s$n ## get number of mesh nodes
  spde <- inla.spde2.matern( mesh_s,alpha=2 )
  ## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
  ## ^ this gives us a linear reduction of \Sigma^{-1} as:
  ## \Sigma^{-1} = \kappa^4 M_0 + 2\kappa^2M_1 + M_2
  ## M_2 = M_1M_0^{-1}M_1
  ## Where the Ms are all sparse matrices stored as "dgTMatrix"
  ## names(spde$param.inla)

  ## use inla helper functions to project the spatial effect from mesh points to data points
  A.proj <- inla.spde.make.A(mesh  = mesh_s,
                             loc   = dt.coords,
                             group = dt.pers)

  ## save relevant objects
  saveRDS(file = sprintf('%s/modeling/inputs/mesh_%i.rds', out.dir, iii), mesh_s)
  saveRDS(file = sprintf('%s/modeling/inputs/spde_%i.rds', out.dir, iii), spde)

  ## ~~~
  ## required for predict
  ## ~~~

  ## get space-time-locs grid to predict onto
  f_orig <- data.table(cbind(coordinates(simple_raster), t=1))
                                        # add time periods
  fullsamplespace <- copy(f_orig)
  if(nperiods > 1){
    for(p in 2:nperiods){
      tmp <- f_orig
      tmp[,t := p]
      fullsamplespace <- rbind(fullsamplespace,tmp)
    }
  }

  ## get surface to project on to
  pcoords <- cbind(x=fullsamplespace$x, y=fullsamplespace$y)
  groups_periods <- fullsamplespace$t

  ## use inla helper functions to project the spatial effect.
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = pcoords,
    group = groups_periods)

  ## extract cell values  from covariates, deal with timevarying covariates here

  
  ## pull out covariates in format we expect them
  ## a list of length periods with a brick of named covariates inside
  new_cl <- list()
  for(p in 1:nperiods){
    new_cl[[p]] <- list()
    for(n in names(cov_list)){
      if(dim(cov_list[[n]])[3] == 1){
        new_cl[[p]][[n]] <- cov_list[[n]]
      }else{
        new_cl[[p]][[n]] <- cov_list[[n]][[p]]
      }
    }
    new_cl[[p]] <- brick(new_cl[[p]])
  }

  cov_vals <- list()
  for(p in 1:nperiods){
    cov_vals[[p]] <- raster::extract(new_cl[[p]], pcoords[1:(nrow(fullsamplespace)/nperiods),])
    cov_vals[[p]] <- (cbind(int = 1, cov_vals[[p]]))
  }
  

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
                    options = c(0, ## if 1, run adreport 
                                1, ## if 1, use priors
                                ifelse(is.null(alpha), 0, 1), # if 1, run with intercept
                                ifelse(is.null(betas), 0, 1), # if 1, run with covs
                                ifelse(is.null(nug.var), 0, 1), # if 1, run with nugget
                                tmb.lik.dict(data.lik)  # if 0, normal data. if 1, binom data lik
                                ),
                    flag = 0 # normalization flag. if 1, use norm trick
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
  
  ## make the autodiff generated liklihood func & gradient
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=rand_effs,
                   map = ADmap, 
                   hessian=TRUE,
                   DLL=templ)

  ## should we use the normalization flag?
  if(data_full$flag == 1){
    obj <- normalize(obj, flag="flag") 
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
  SD0 = TMB::sdreport(obj, getJointPrecision=TRUE)##,
##                      bias.correct = TRUE,
##                      bias.correct.control = list(sd = TRUE))
  tmb_total_fit_time <- proc.time()[3] - ptm 
  tmb_sdreport_time <-  tmb_total_fit_time - fit_time_tmb

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
    message('\nTMB PRECISION IS NOT! PD - mapping to nearest PD precision ')
    message('\nTMB PRECISION IS NOT! PD - mapping to nearest PD precision ')
    message('\nTMB PRECISION IS NOT! PD - mapping to nearest PD precision ')
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
  log_nugget_sigma <- tmb_draws[parnames == 'log_nugget_sigma', ]

  ## values of S at each cell (long by nperiods)
  ## rows: pixels, cols: posterior draws
  pred_tmb <- as.matrix(A.pred %*% epsilon_tmb_draws)

  if(!is.null(alpha)){
    ## add on intercept, one alpha draw per row
    pred_tmb <- sweep(pred_tmb, 2, alpha_tmb_draws, '+')
  }

  if(!is.null(betas)){
    ## add on ovariate values by draw
    tmb_vals <- list()
    for(p in 1:nperiods) tmb_vals[[p]] <- cov_vals[[p]] %*% betas_tmb_draws

    cell_b <- do.call(rbind, tmb_vals)

    ## add together linear and st components
    pred_tmb <- cell_b + pred_tmb
  }

  ## save prediction timing
  totalpredict_time_tmb <- proc.time()[3] - ptm2

  ## #######
  ## SAVE ##
  ## #######

  ## make some useful files first
  summ_tmb <- cbind(median = (apply(pred_tmb,1,median)),
                    sd = (apply(pred_tmb,1,sd)))

  ras_med_tmb  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,1]),
                                            t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")
  ras_sdv_tmb  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,2]),
                                            t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")

  saveRDS(file = sprintf('%s/modeling/tmb/outputs/tmb_preds_median_raster_%i.rds', out.dir, iii), object = ras_med_tmb)
  saveRDS(file = sprintf('%s/modeling/tmb/outputs/tmb_preds_stdev_raster_%i.rds', out.dir, iii), object = ras_sdv_tmb)


  ## #######
  ## #######
  ## INLA ##
  ## #######
  ## #######


  ## ######
  ## FIT ##
  ## ######
  A <- inla.spde.make.A(
    mesh = mesh_s,
    loc = dt.coords,
    group = dt.pers
  )
  space   <- inla.spde.make.index("space",
                                  n.spde = spde$n.spde,
                                  n.group = nperiods)

  inla.covs <- covs$name
  design_matrix <- data.frame(int = rep(1, nrow(dt)), dt[, inla.covs, with=F])
  stack.obs <- inla.stack(tag='est',
                          data=list(Y=dt$Y), ## response
                          A=list(A,1), ## proj matrix, not sure what the 1 is for
                          effects=list(
                            space,
                            design_matrix)
                          )

  formula <- formula(paste('Y ~ -1',
                           ifelse(is.null(alpha), '', 'int'), 
                           ifelse(is.null(betas), '', (paste(inla.covs, collapse = ' + '))),
                           'f(space, model = spde, group = space.group, control.group = list(model = \'ar1\'))',
                           sep = ' + '))


  
  ## function to convert from data lik string to integer
  ## allows easily adding more options even though overkill for just 2
  inla.lik.dict <- function(x){
    dict <- list(normal = 'normal',
                 binom = 'binomial')
    dict[[x]]
  }
  inla.setOption("enable.inla.argument.weights", TRUE)

  ptm <- proc.time()[3] 
  res_fit <- inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           ## link = 1, ## removed after looking at NMM
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla'),
                  control.inla = list(strategy = inla.approx,
                                      int.strategy = inla.int.strat ##,
                                      ## h = 1e-3, ## removed after looking at NMM
                                      ## tolerance = 1e-6 ## removed after looking at NMM
                                      ),
                  control.compute=list(config = TRUE),
                  family = inla.lik.dict(data.lik),
                  num.threads = cores, #
                  Ntrials = dt$N,
                  ## weights = rep(1, nrow(dt)),
                  verbose = TRUE,
                  keep = FALSE)
  fit_time_inla <- proc.time()[3] - ptm



  ## ##########
  ## PREDICT ##
  ## ##########

  ptm <- proc.time()[3]
  inla_draws <- inla.posterior.sample(ndraws, res_fit)
  inla_get_draws_time <- proc.time()[3] - ptm

  ## get parameter names
  par_names <- rownames(inla_draws[[1]]$latent)

  ## index to spatial field and linear coefficient samples
  s_idx <- grep('^space.*', par_names)
  l_idx <- which(!c(1:length(par_names)) %in% grep('^space.*|Predictor', par_names))

  ## get spatial draws as matrices and project to deaws at locations 
  pred_s <- sapply(inla_draws, function (x) x$latent[s_idx])
  pred_inla <- as.matrix(A.pred %*% pred_s)

  ## get intercept and coef draws and convert to covariate effects
  if(!is.null(alpha) | !is.null(betas)){
    pred_l <- sapply(inla_draws, function (x) x$latent[l_idx])
    if(!is.matrix(pred_l)){
      pred_l <- matrix(pred_l, ncol = length(pred_l))
    }
    rownames(pred_l) <- res_fit$names.fixed

    ## extract cell values  from covariates, deal with timevarying covariates here
    non_space_names <- par_names[l_idx]
    cov_effs <- list()
    for(p in 1:nperiods)  cov_effs[[p]] <- cov_vals[[p]][, non_space_names] %*% pred_l
    
    cov_effs <- do.call(rbind, cov_effs)
    
    pred_inla <- pred_inla + cov_effs
  }

  ## make them into time bins
  len = nrow(pred_inla)/nperiods
  totalpredict_time_inla <- proc.time()[3] - ptm

  ## #######
  ## SAVE ##
  ## #######

  ## make some useful files first
  summ_inla <- cbind(median = (apply(pred_inla,1,median)),
                     sd = (apply(pred_inla,1,sd)))

  ras_med_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,1]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
  ras_sdv_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,2]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")

  saveRDS(file = sprintf('%s/modeling/inla/outputs/inla_preds_median_raster_%i.rds', out.dir, iii), object = ras_med_inla)
  saveRDS(file = sprintf('%s/modeling/inla/outputs/inla_preds_stdev_raster_%i.rds', out.dir, iii), object = ras_sdv_inla)



  ## #############
  ## #############
  ## VALIDATION ##
  ## #############
  ## #############

  ## 1) summarize fitted params
  ## 2) big plots showing difference in fits
  ## 3) calcualte and summarize predictive metrics

  source('./5_run_validation.R')


} ## end iii loop repeating iterations over 1:Nsim

write.csv(complete.surface.metrics, sprintf('%s/validation/surface_metrics_complete.csv',out.dir))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## (i) compare data to fitted surface
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#############
## SCRATCH ##
#############

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## !! moved into realistic_sim_utils.R/sim.realistic.data() on 0ct 3, 2018 !!
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## top_pop_urban <- 0.01
## urban_strat <- 0.4
## urban_thresh <- quantile(probs = (1 - top_pop_urban), na.omit(values(pop_raster)))
## u_r_raster <- pop_raster[[1]] ## urban is 1, rural is 0
## u_r_raster[pop_raster[[1]] < urban_thresh] <- 0
## u_r_raster[pop_raster[[1]] >= urban_thresh] <- 1

## pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)
## u_r.pts <- rasterToPoints(u_r_raster, spatial = TRUE)

## ## reproject sp obj
## geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
## pix.pts <- spTransform(pix.pts, CRS(geo.prj))
## u_r.pts <- spTransform(u_r.pts, CRS(geo.prj))
## proj4string(pix.pts)
## proj4string(u_r.pts)

## ## get coords
## pix.pts@data <- data.frame(pix.pts@data, long=coordinates(pix.pts)[,1],
##                            lat=coordinates(pix.pts)[,2])
## pix.pts.numeric <- as.data.frame(pix.pts@data)

## u_r.pts@data <- data.frame(u_r.pts@data, long=coordinates(u_r.pts)[,1],
##                            lat=coordinates(u_r.pts)[,2])
## u_r.pts.numeric <- as.data.frame(u_r.pts@data)

## sim.rows <- sample(x = 1:nrow(pix.pts.numeric), size = n.clust * length(year_list),
##                    replace = TRUE)
## sim.dat <- as.data.table(pix.pts.numeric[, -1])
## sim.dat <- sim.dat[sim.rows, ]


## u.rows <- sample(x = which(u_r.pts.numeric[, 1] == 1), size = round(n.clust * urban_strat),
##                  replace = TRUE)
## r.rows <- sample(x = which(u_r.pts.numeric[, 1] == 0), size = round(n.clust * (1 - urban_strat)),
##                  replace = TRUE)

## sim.dat <- as.data.table(pix.pts.numeric[, -1])
## sim.dat <- sim.dat[c(u.rows, r.rows), ]


## pdf(sprintf('%s/validation/pop_strat.pdf',out.dir), height=10,width=20)
## par(mfrow = c(1, 2))
## raster::plot(pop_raster[[1]])
## points(sim.dat)
## raster::plot(u_r_raster)
## dev.off()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



