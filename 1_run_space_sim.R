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
betas           <- eval(parse(text = as.character(loopvars[par.iter, 5])))

alpha           <- as.numeric(loopvars[par.iter, 6])
sp.range        <- as.numeric(loopvars[par.iter, 7])
sp.var          <- as.numeric(loopvars[par.iter, 8])
sp.alpha        <- as.numeric(loopvars[par.iter, 9])
nug.var         <- as.numeric(loopvars[par.iter, 10])

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

## make an object with true param values
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
true.par.vec <- c(alpha, betas, logtau, logkappa, trho_trans)
names(true.par.vec) <- c(rep('alpha_j', length(betas) + 1), 'logtau', 'logkappa', 'trho_trans')
if(length(year_list) == 1) true.par.vec <- true.par.vec[-length(true.par.vec)]

saveRDS(file = sprintf('%s/simulated_obj/true_param_table.rds', out.dir),
        object = true.params)


###########
###########
###########

for(iii in 1:Nsim){ ## repeat Nsim times

  ## if(iii == 1){ ## first time, must load covs, after that, we can reuse them
    sim.obj <- sim.realistic.data(reg = reg,
                                  year_list = year_list,
                                  betas = betas,
                                  sp.kappa = sp.kappa,
                                  sp.alpha = sp.alpha,
                                  t.rho = t.rho,
                                  nug.var = nug.var, 
                                  n.clust = n.clust,
                                  m.clust = m.clust,
                                  covs = covs,
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
  cov_list <- covs.gp$gp <- NULL
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
  dir.create(sprintf('%s/modeling/inputs', out.dir), recursive = TRUE)
  dir.create(sprintf('%s/modeling/tmb/outputs', out.dir), recursive = TRUE)
  dir.create(sprintf('%s/modeling/inla/outputs', out.dir), recursive = TRUE)

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

  ## pull out covariates in format we expect them
  ## a list of length periods with a brick of named covariates inside
  cov_list <- covs.gp
  cov_list$gp <- NULL
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
  X_xp = as.matrix(cbind(1, dt[,covs[, name], with=FALSE])) ## design mat

  templ <- "model_space"
  setwd("/homes/azimmer/tmb_transition/realistic_sims")
  TMB::compile(paste0('./', templ,".cpp"))
  dyn.load( dynlib(templ) )

  ## setup data to feed into the model
  data_full <- list(num_i = nrow(dt),  # Total number of observations
                    num_s = mesh_s$n,  # Number of vertices in SPDE mesh
                    y_i   = dt[,Y],    # Number of observed deaths in the cluster
                    n_i = dt[,N],    # Number of exposures in the cluster
                    X_ij  = X_xp,               # Covariate design matrix
                    M0    = spde$param.inla$M0, # SPDE sparse matrix
                    M1    = spde$param.inla$M1, # SPDE sparse matrix
                    M2    = spde$param.inla$M2, # SPDE sparse matrix
                    Aproj = A.proj,             # Projection matrix
                    options = c(1, # use priors 
                                0, # adreport
                                1  # fit with nugget 
                                ),
                    flag = 1 # normalization flag
                    )

  ## Specify starting values for TMB parameters
  tmb_params <- list(alpha_j   = rep(0,ncol(X_xp)), # Alphas for FE parameters
                     log_tau   = 1.0, # Log inverse of tau (Epsilon)
                     log_kappa = 1.0, # Matern range parameter
                     Epsilon_s = matrix(1,nrow=nodes,ncol=1), # GP locations
                     log_nugget_sigma = -1.0, # log of nugget sd
                     nug_i = rep(0, nrow(dt)) # vector of nugget random effects
                     )


  ## make the autodiff generated liklihood func & gradient
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=c('Epsilon_s', 'nug_i'),
                   hessian=TRUE,
                   DLL=templ)

  ## should we use the normalization flag?
  obj <- normalize(obj, flag="flag")


  ## Run optimizer
  ptm <- proc.time()[3]
  opt0 <- do.call("nlminb",list(start       =    obj$par,
                                objective   =    obj$fn,
                                gradient    =    obj$gr,
                                lower = rep(-10, ncol(X_xp) + 3),
                                upper = rep(10, ncol(X_xp) + 3),
                                control     =    list(trace=1)))
  fit_time_tmb <- proc.time()[3] - ptm

  ## Get standard errors
  SD0 = TMB::sdreport(obj, getJointPrecision=TRUE,
                      bias.correct = TRUE,
                      bias.correct.control = list(sd = TRUE))
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
  alpha_tmb_draws    <- tmb_draws[parnames == 'alpha_j',]
  log_kappa_tmb_draws <- tmb_draws[parnames == 'log_kappa',]
  log_tau_tmb_draws  <- tmb_draws[parnames == 'log_tau',]
  log_nugget_sigma <- tmb_draws[parnames == 'log_nugget_sigma', ]


  ## values of S at each cell (long by nperiods)
  cell_s <- as.matrix(A.pred %*% epsilon_tmb_draws)

  ## covariate values by alpha draws
  tmb_vals <- list()
  for(p in 1:nperiods)  tmb_vals[[p]] <- cov_vals[[p]] %*% alpha_tmb_draws

  cell_l <- do.call(rbind,tmb_vals)

  ## add together linear and st components
  pred_tmb <- cell_l + cell_s

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

  formula <- formula(paste0('Y ~ -1+int+',
  (paste(inla.covs, collapse = ' + ')),
  ' + f(space, model = spde, group = space.group, control.group = list(model = \'ar1\'))'))

  ptm <- proc.time()[3]

  inla.setOption("enable.inla.argument.weights", TRUE)
  res_fit <- inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           ## link = 1, ## removed after looking at NMM
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla'),
                  control.inla = list(strategy = inla.approx,
                                      int.strategy = inla.int.strat,
                                      ## h = 1e-3, ## removed after looking at NMMn
                                      ## tolerance = 1e-6 ## removed after looking at NMM
                                      ),
                  control.compute=list(config = TRUE),
                  family = 'binomial',
                  num.threads = cores, #
                  Ntrials = dt$N,
                  weights = rep(1, nrow(dt)),
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

  ## get samples as matrices
  pred_s <- sapply(inla_draws, function (x) x$latent[s_idx])
  pred_l <- sapply(inla_draws, function (x) x$latent[l_idx])
  rownames(pred_l) <- res_fit$names.fixed

  ## get samples of s for all coo locations
  s <- as.matrix(A.pred %*% pred_s)

  ## extract cell values  from covariates, deal with timevarying covariates here
  inla_vals <- list()
  for(p in 1:nperiods)  inla_vals[[p]] <- cov_vals[[p]] %*% pred_l

  l <- do.call(rbind,inla_vals)

  pred_inla <- s+l

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

  dir.create(sprintf('%s/validation', out.dir))

  ## ###################################
  ## 1) summarize fitted param values ##
  ## ###################################

  res <- data.table(st_mesh_nodes = rep(nrow(epsilon_tmb_draws),2))
  res[,cores           := rep(cores,2)]
  res[,s_mesh_max_edge := rep(eval(parse(text = mesh_s_max_edge))[1],2)]
  res[,periods         := c(nperiods,nperiods)]
  res[,draws           := c(ndraws,ndraws)]

  ## time variables
  res[,fit_time  := c(fit_time_inla,fit_time_tmb)]
  res[,pred_time := c(totalpredict_time_inla,totalpredict_time_tmb)]
  res[,pt_tmb_sdreport_time := c(NA,tmb_sdreport_time)]
  res[,pt_get_draws_time := c(inla_get_draws_time,tmb_get_draws_time)]

  ## fe coefficients
  for(i in 1:length(res_fit$names.fixed)){
    fn <- res_fit$names.fixed[i]
    res[[paste0('fixedeff_',fn,'_mean')]] <- c(res_fit$summary.fixed$mean[i], SD0$par.fixed[i])
    res[[paste0('fixedeff_',fn,'_sd')]]   <- c(res_fit$summary.fixed$sd[i], sqrt(SD0$cov.fixed[i,i]))
  }

  ## hyperparameters
  res[,hyperpar_logtau_mean := c(res_fit$summary.hyperpar[1,1], SD0$par.fixed['log_tau']) ]
  res[,hyperpar_logtau_sd := c(res_fit$summary.hyperpar[1,2], sqrt(SD0$cov.fixed['log_tau','log_tau'])) ]

  res[,hyperpar_logkappa_mean := c(res_fit$summary.hyperpar[2,1], SD0$par.fixed['log_kappa']) ]
  res[,hyperpar_logkappa_sd := c(res_fit$summary.hyperpar[2,2], sqrt(SD0$cov.fixed['log_kappa','log_kappa'])) ]

  ## res[,hyperpar_rho_mean := c(res_fit$summary.hyperpar[3,1], SD0$par.fixed['trho']) ]
  ## res[,hyperpar_rho_sd := c(res_fit$summary.hyperpar[3,2], sqrt(SD0$cov.fixed['trho','trho'])) ]

  ## truth
  res.true.params <- c(rep(NA, 9),
                       alpha, NA, ## intercept
                       c(rbind(betas, rep(NA, length(betas)))), ## fixed effect coefs
                       -0.5 * log(4 * pi * sp.var * sp.kappa^2), NA, ## log tau
                       log(sp.kappa), NA) ## log kappa
  if(nperiods > 1){
    res.true.params <- c(res.true.params, c(t.rho, NA))
  }

  rr <- data.table(item=colnames(res))
  rr <- cbind(rr,res.true.params, t(res))
  names(rr) <- c('quantity','TRUE', 'R-INLA','TMB')
  rr$diff <- rr[,3]-rr[,4]

  ## we can now plot this table with: grid.table(rr)

  ## ####################
  ## 2) setup big plot ##
  ## ####################
  pdf(sprintf('%s/validation/inla_tmb_summary_comparison_plot_%i.pdf',out.dir, iii), height=15,width=30)


  ## ~~~
  ## plot the table of results
  ## ~~~
  grid.table(rr)

  ## ~~~
  ## plot priors and posteriors
  ## ~~~

  ## assume:
  ## 1) alphas are normal
  ## 2) logtau and logkappa are normal
  ## 3) log nugget precision is gamma

  ## NOTE: also assume that intercept is the first 'beta' and that betas are listed first!
  params <- c(rep('beta', length(betas) + 1), 'logkappa', 'logtau')#, 'nug.var')
  num.dists <- length(params)

  par(mfrow = rep(ceiling(sqrt(num.dists)), 2))

  for(ii in 1:num.dists){

    param <- params[ii]
    
    ## get prior curves and posterior draws
    if(param == 'beta'){
      prior.mean <- alphaj.pri[1]
      prior.sd   <- sqrt(alphaj.pri[2])
      xlim       <- c(-3, 3) ## TODO for all plots ## prior.mean + c(-4, 4) * prior.sd
      x.prior    <- seq(xlim[1], xlim[2], by = 0.01)
      y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

      tmb.post.draws  <- alpha_tmb_draws[ii, ]
      inla.post.draws <- pred_l[ii, ]
      tmb.post.median <- median(tmb.post.draws)
      inla.post.median <- median(inla.post.draws)
      param.name <- names(res_fit$marginals.fixed)[ii]

      if(param.name == 'int'){
        true.val <- alpha
      }else{
        true.val <- betas[ii - 1]
      }

    }
    if(param == 'logkappa'){
      mesh.pri   <- param2.matern.orig(mesh_s) ## theta1 = logtau, theta2 = logkappa 
      prior.mean <- mesh.pri$theta.prior.mean[2]
      prior.prec <- mesh.pri$theta.prior.prec[2, 2]
      prior.sd   <- 1 / sqrt(prior.prec)
      x.prior    <- seq(xlim[1], xlim[2], by = 0.01)
      y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

      tmb.post.draws  <- log_kappa_tmb_draws
      inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[['Theta2 for space']][, 1],
                                      size = ndraws,
                                      replace = TRUE, 
                                      res_fit$marginals.hyperpar[['Theta2 for space']][, 2])
      tmb.post.median  <- median(tmb.post.draws)
      inla.post.median <- median(inla.post.draws)
      param.name <- "log kappa"
      true.val <- logkappa
    }
    if(param == 'logtau'){
      mesh.pri   <- param2.matern.orig(mesh_s) ## theta1 = logtau, theta2 = logkappa 
      prior.mean <- mesh.pri$theta.prior.mean[1]
      prior.prec <- mesh.pri$theta.prior.prec[1, 1]
      prior.sd   <- 1 / sqrt(prior.prec)
      x.prior    <- seq(xlim[1], xlim[2], by = 0.01)
      y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

      tmb.post.draws <- log_tau_tmb_draws
      inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[['Theta1 for space']][, 1],
                                      size = ndraws,
                                      replace = TRUE, 
                                      res_fit$marginals.hyperpar[['Theta1 for space']][, 2])
      tmb.post.median <- median(tmb.post.draws)
      inla.post.median <- median(inla.post.draws)
      param.name <- "log tau"
      true.val <- logtau
    }
    if(param == 'nug.var'){
      prior.shape  <- nug.pri[1]
      prior.iscale <- nug.pri[2]
      x.prior <- seq(0, 5, by = 0.01) ## TODO automate this?
      y.prior <- dgamma(x.prior, shape = prior.shape, scale = 1 / prior.iscale)

    }

    ## get posterior samples (we'll use density curves)
    tmb.dens <- density(tmb.post.draws)
    inla.dens <- density(inla.post.draws)
    xrange <- range(c(x.prior, tmb.dens$x, inla.dens$x))
    yrange <- range(c(y.prior, tmb.dens$y, inla.dens$y))

    prior.col <- "black"
    tmb.col <- "red"
    inla.col <- "blue"

    ## setup plot and plot prior
    plot(x.prior, y.prior, pch = ".", col = prior.col, main = param.name,
         xlim = xrange, ylim = yrange)
    lines(x.prior, y.prior, col = prior.col)

    ## plot tmb post and data
    lines(tmb.dens$x,tmb.dens$y, col = tmb.col)
    points(tmb.post.draws, rep(0, ndraws), col = alpha(tmb.col, 0.25), cex = 2, pch = '|')
    abline(v = tmb.post.median, col = tmb.col, lwd = 2)
    
    ## plot inla post and data
    lines(inla.dens$x,inla.dens$y, col = inla.col)
    points(inla.post.draws, rep(0, ndraws), col = alpha(inla.col, 0.25), cex = 2, pch = '|')
    abline(v = inla.post.median, col = inla.col, lwd = 2)

    ## plot truth
    abline(v = true.val, col = prior.col, lwd = 2)
    
    ## add legend
    legend("topright", legend = c("prior", "tmb", "inla"), col = c(prior.col, tmb.col, inla.col), lwd = rep(2, 3))
  }




  ## plot results in logit space or in prevalence space?
  plot.in.logit.space <- FALSE

  ## setup layout of main plots
  layout(matrix(1:(nperiods * 2 * 5), (nperiods * 2), 5, byrow = TRUE))

  ## randomly select pixels for plotting tmb v inla scatter
  samp <- sample(cellIdx(ras_med_inla[[1]]),1e4) 

  for(i in 1:nperiods){ ## for s-t

    ## TODO: allow plotting in logit space
    if(plot.in.logit.space){
      
    }else{

      for(thing in c('median','stdev')){ 
        
        if(thing=='median'){
          rinla <- ras_med_inla[[i]]
          rtmb  <- ras_med_tmb[[i]]
          true  <- invlogit(true.rast[[i]])
        }
        if(thing=='stdev'){
          rinla <- ras_sdv_inla[[i]]
          rtmb  <- ras_sdv_tmb[[i]]
        }
        
        tmp <- subset(dt, period_id==i) ## for s-t
        
        ## rasters: true, tmb, inla, inla - tmb with data locs
        par(mar = c(0, 0, 1.4, 2),bty='n')
        maxes <- max(c(as.vector(rtmb), as.vector(rinla), as.vector(true)), na.rm=TRUE)
        mins  <- min(c(as.vector(rtmb), as.vector(rinla), as.vector(true)), na.rm=TRUE)
        zrange <- c(mins, maxes)
        if(thing == 'median'){
          plot(true, maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend=T, legend.width = 3, main=paste0('TRUE ', thing) ,zlim=zrange)
          plot(simple_polygon, add = T)
        }else{
          plot.new();abline(a = 0, b = 1, lwd = 2);abline(a = 1, b = -1, lwd = 2) ## no true sd surface to plot....
        }
        plot(rtmb,  maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend.width = 3,
             legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main=paste0('TMB: ', thing), zlim=zrange)
        plot(simple_polygon, add = T)
        plot(rinla, maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend=FALSE, main=paste0('R-INLA: ', thing), zlim=zrange)
        plot(simple_polygon, add = T)
        plot(rinla-rtmb, maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend=T,  legend.width = 3, main=paste0('DIFFERENCE: ',thing))
        plot(simple_polygon, add = T)
        points( x=tmp$long,y=tmp$lat, pch=19, cex=(tmp$N / max(tmp$N)) )
        
        ## pixel scatter
        par(mar = c(4, 4, 2, 2),bty='n')
        plot(x=as.vector(rinla)[samp],y=as.vector(rtmb)[samp],xlab='R-INLA',ylab='TMB',cex=.01,pch=19,main=paste0('(sub)SCATTER COMPARE: ', thing))
        lines(x=zrange,y=zrange,col='red')

        ## residual
        ##  tmp$inla<-extract(ras_med_inla[[i]],cbind(tmp$longitude,y=tmp$latitude))
        ##  tmp$tmb<-extract(ras_med_tmb[[i]],cbind(tmp$longitude,y=tmp$latitude))
        ##  tmp$dat <- tmp$died/tmp$N
        ##  tmp$resid_inla <- tmp$dat-tmp$inla
        ##  tmp$resid_tmb  <- tmp$dat-tmp$tmb
        ##  tmp<-subset(tmp,dat<quantile(tmp$dat,.9))
        ##  plot(x=tmp$dat,y=tmp$resid_inla, pch=19,col='red',cex=.1,main='RESIDUALS')
        ##  points(x=tmp$dat,y=tmp$resid_tmb, pch=19,col='blue',cex=.1)
      }
    }
  }


  ## now make caterpillar plots

  layout(matrix(1, 1, 1, byrow = TRUE))

  ## Compare mean and distribution of random effects
  summ_gp_tmb  <- t(cbind((apply(epsilon_tmb_draws,1,quantile,probs=c(.1,.5,.9)))))
  summ_gp_inla <- t(cbind((apply(pred_s,1,quantile,probs=c(.1,.5,.9)))))
  ## all time-space random effects

  plot_d <- data.table(tmb_median = summ_gp_tmb[,2],inla_median = summ_gp_inla[,2],
                       tmb_low    = summ_gp_tmb[,1],inla_low    = summ_gp_inla[,1],
                       tmb_up     = summ_gp_tmb[,3],inla_up     = summ_gp_inla[,3])

  plot_d$period <- factor(rep(1:nperiods, each=nrow(plot_d)/nperiods))
  plot_d$loc    <- rep(1:(nrow(plot_d)/nperiods), rep=nperiods)

  if(nrow(plot_d)>2500)
    plot_d <- plot_d[sample(nrow(plot_d),2500,replace=F),]


  ggplot(plot_d, aes(x=tmb_median,y=inla_median,col=period)) + theme_bw() +
  geom_point() + geom_line(aes(group=loc)) + geom_abline(intercept=0,slope=1,col='red') +
  ggtitle('Posterior Medians of Random Effects at Mesh Nodes, TMB v R-INLA. Connected dots same location different periods. ')

  ## plot locations where they are different, are they near or far from data?
  plot_d[, absdiff := abs(tmb_median-inla_median)]
  nodelocs <- do.call("rbind", replicate(4, mesh_s$loc, simplify = FALSE))
  biggdiff <- unique(nodelocs[which(plot_d$absdiff>quantile(plot_d$absdiff,prob=0.80)),])

  nodelocs <- cbind(nodelocs,plot_d)
  if(nrow(nodelocs)>2500)
    nodelocs <- nodelocs[sample(nrow(nodelocs),2500,replace=FALSE),]

  par(mfrow=rep(ceiling(sqrt(nperiods)),2))
  for(i in 1:nperiods){
    plot(simple_polygon, main='Mesh nodes sized by abs difference TMB and R-INLA')
    points(x=dt$longitude[dt$period_id==i],y=dt$latitude[dt$period_id==i], pch=19, cex=0.1)
    points(x=nodelocs$V1[nodelocs$period==i],y=nodelocs$V2[nodelocs$period==i], pch=1, cex=nodelocs$absdiff[nodelocs$period==i]*5, col='red')
    ## add data locations
    points( x=tmp$long,y=tmp$lat, cex=(tmp$N / max(tmp$N)), pch = 16)

  }

  ## catterpillar plot
  plot_d <- plot_d[order(period,tmb_median)]
  plot_d[,i := seq(1,.N), by = period]
  ggplot(plot_d, aes(i, tmb_median, col=i)) + theme_bw() + # [seq(1, nrow(plot_d), 5)]
  geom_linerange(aes(ymin = tmb_low, ymax = tmb_up), col='blue', size=.8, alpha=.3) +
  geom_linerange(aes(x=i,ymin = inla_low, ymax = inla_up), col='red', size=.8, alpha=.3) +
  facet_wrap(~period) +
  ggtitle('Comparison of random effects (10% to 90% quantiles) ... RED == R-INLA ... BLUE == TMB')


  dev.off()

  ## ###############################################
  ## 3) generate and summarize predictive metrics ##  
  ## ###############################################

  ## there are two types of predictive metrics we might look at:
  ##  (i) metrics comparing true surfaces
  ## (ii) metrics comparing data to estimated surfaces

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## (i) compare true surface to fitted surface
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  all.preds <- data.table(rbind(pred_tmb, pred_inla))

  ## make a data.table with prediction draws, model type, and truth
  ## NOTE! the truth and the median.fit are in logit-space!
  d <- data.table(truth        = rep(values(true.rast), 2),
                  truth.p      = ilogit(rep(values(true.rast), 2)),
                  model        = c(rep('tmb', nrow(pred_tmb)),
                                   rep('inla', nrow(pred_inla))), 
                  median.fit   = apply(all.preds, 1, median),
                  median.p.fit = apply(ilogit(all.preds), 1, median))


  ## get some coverage probs
  coverage_probs <- c(50, 80, 95)
  for(c in coverage_probs){
    message(paste0('For ',c,'% coverage.'))
    coverage <- c/100
    li       <- apply(all.preds, 1, quantile, p=(1-coverage)/2, na.rm=T)
    ui       <- apply(all.preds, 1, quantile, p=coverage+(1-coverage)/2, na.rm=T)
    d[,paste0('pixel_covered_',c)] = d[['truth']]>=li & d[['truth']] <= ui
  }

  ## get error and pred var
  d[, error := truth - median.fit]
  d[, var := apply(all.preds, 1, var)]
  d[, error.p := truth.p - median.p.fit]
  d[, var.p := apply(ilogit(all.preds), 1, var)]


  ## summarize metrics across surfaces and models

  ## simple correlation function
  my.cor <- function(x, y){
    x <- na.omit(x)
    y <- na.omit(y)
    s.x <- sum(x)
    s.x2 <- sum(x ^ 2)
    s.y <- sum(y)
    s.y2 <- sum(y ^ 2)
    s.xy <- sum(x * y)
    n <- length(x)
    return((n * s.xy - s.x * s.y) / (sqrt(n * s.x2 - s.x ^ 2) * sqrt(n * s.y2 - s.y ^ 2)))
  }

  ## continuous-rank probability score
  ## NOTE! since this is a normal distribution, I calcualte it on the logit scale which should be close to gaussian in our predictions
  crpsNormal <- function(truth, my.est, my.var){
    
    sig = sqrt(my.var)
    x0 <- (truth - my.est) / sig
    res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
    
    ## sign as in Held (2008)
    res <- -res
    
    return(res)
  }

  surface.metrics <- data.table(cbind(mean_l     = d[, .(truth = mean(truth, na.rm = T)), by = c('model')], 
                                      mean_lhat  = d[, .(est = mean(median.fit, na.rm = T)), by = c('model')]$est, 
                                      bias       = d[, .(bias = mean(error, na.rm = T)), by = c('model')]$bias,
                                      rmse       = d[, .(rmse = sqrt(mean(error ^ 2, na.rm = T))), by = c('model')]$rmse,
                                      cor        = d[, .(cor = my.cor(truth, median.fit)), by = c('model')]$cor,
                                      CoV        = d[, .(cov = mean(error / var, na.rm = T)), by = c('model')]$cov,
                                      crps       = d[, .(crps = mean(crpsNormal(truth, median.fit, var), na.rm = T)), by = c('model')]$crps, 
                                      mean_p     = d[, .(truth.p = mean(truth.p, na.rm = T)), by = c('model')]$truth.p, 
                                      mean_phat  = d[, .(est.p = mean(median.p.fit, na.rm = T)), by = c('model')]$est.p, 
                                      bias_p     = d[, .(bias.p = mean(error.p, na.rm = T)), by = c('model')]$bias.p,
                                      rmse_p     = d[, .(rmse.p = sqrt(mean(error.p ^ 2, na.rm = T))), by = c('model')]$rmse.p,
                                      cor_p      = d[, .(cor.p = my.cor(truth.p, median.p.fit)), by = c('model')]$cor.p,
                                      CoV_p      = d[, .(cov.p = mean(error.p / var.p, na.rm = T)), by = c('model')]$cov.p,
                                      cov50      = d[, .(cov50 = mean(pixel_covered_50, na.rm = T)), by = c('model')]$cov50,
                                      cov80      = d[, .(cov80 = mean(pixel_covered_80, na.rm = T)), by = c('model')]$cov80,
                                      cov95      = d[, .(cov95 = mean(pixel_covered_95, na.rm = T)), by = c('model')]$cov95
                                      ))


  surface.metrics[, iter := iii]

  write.csv(surface.metrics, sprintf('%s/validation/surface_metrics_%i.csv',out.dir, iii))
  
  if(iii == 1){
    complete.surface.metrics <- surface.metrics
  }else{
    complete.surface.metrics <- rbind(complete.surface.metrics, surface.metrics)
  }

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



