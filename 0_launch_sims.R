## this script can be used to launch 1_run_simulation.R in parallel on the IHME cluster
## written by aoz
## dec 1 2018

## DO THIS!
################################################################################
## ADD A NOTE! to help identify what you were doing with this run
logging_note <- 'no nug and no covs'
################################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#############################################
## setup the environment for singularity R ##
#############################################

## Set core_repo location and tmb_repo loc
user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_transition', user)
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

## now we can setup our main directory to save these results and log our note
run_date <- make_time_stamp(TRUE)
out.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', run_date)
dir.create(out.dir)
fileConn <- file(sprintf("%s/run_notes.txt", out.dir))
writeLines(logging_note, fileConn)
close(fileConn)

## Now we can switch to the TMB repo
setwd(tmb_repo)
if(pull_tmb_git) system(sprintf('cd %s\ngit pull %s %s', core_repo, remote, branch))
source('./realistic_sims/realistic_sim_utils.R')

###############################
## setup things to loop over ##
###############################

## NOTES
## this list should align with the args in 1_run_simulation.R
## NULLs can't be passed this way, so NAs are stand-ins for NULL and this gets fixed in 1_run_space_sim.R


## loopvars 1
reg <- 'nga'

## loopvars 
year_list <- 2000 

## loopvars 3
cov_names <- "c('access2', 'distrivers', 'evi'   , 'mapincidence')" 

## loopvars 4
cov_measures <- "c('mean'   , 'mean'      , 'median', 'mean')"

## loopvars 5
betas <- NA ## "c(.5, -1, 1, -.5)"

## loopvars 6
alpha <- 0

## loopvars 7
sp.range <-  sqrt(8)     ## kappa=sqrt(8)/sp.range, so sp.range=sqrt(8) -> kappa=1 -> log(kappa)=0 (for R^2 domain)

## loopvars 8
sp.var <- 0.5            ## sp.var = 1/(4*pi*kappa^2*tau^2) (for R^2 domain)

## loopvars 9
sp.alpha <- 2.0          ## matern smoothness = sp.alpha - 1 (for R^2 domain)

## loopvars 10
nug.var <- NA ## .5 ^ 2        ## nugget variance

## loopvars 11
t.rho <-  0.8            ## annual temporal auto-corr

## loopvars 12
mesh_s_max_edge <- c("c(0.2,5)",
                     "c(0.3,5)")

## loopvars 13
n.clust <-  50                   ## clusters PER TIME slice

## loopvars 14
m.clust <- 35                    ## mean number of obs per cluster (poisson)

## loopvars 15
## each entry must be a character string with the syntax for a 3 element R list containing:
## 1) obs.loc.strat: (either 'rand' or 'pop.strat')
## 2) urban.pop.pct:   a number between 0 and 100. the % of population that belongs to urban pixels
## 3) urban.strat.pct: a number between 0 and 100. the % of observations that should come fom urban pixels
sample.strat <- "list(obs.loc.strat='rand',
                      urban.pop.pct=5,
                      urban.strat.pct=40)"  ## random or by population for now. ## TODO cluser design

## loopvars 16
cores <- 5

## loopvars 17
ndraws <- 250

## loopvars 18
alphaj.pri <- "c(0, 3)" ## normal mean and sd

## loopvars 19
nug.pri <- "c(1, 1e-5)" ## gamma for nug preciion with shape and inv-scale

## loopvars 20
inla.int.strat <- 'eb' ## can be 'eb', 'ccd', or 'grid'

## loopvars 21
inla.approx <- 'simplified.laplace' ## can be 'gaussian', 'simplified.laplace' (default) or 'laplace'

## loopvars 22
Nsim <- 5 ## number of times to repeat simulation

## loopvars 23
data.lik <- 'normal' ## either 'binom' or 'normal'

## loopvars 24
sd.norm <- 0.1

## loopvars 25

## TODO always add all vars to exand.grid() 
loopvars <- expand.grid(reg, ## 1
                        year_list,
                        cov_names,
                        cov_measures,
                        betas, ## 5
                        alpha,
                        sp.range,
                        sp.var,
                        sp.alpha,
                        nug.var, ## 10
                        t.rho,
                        mesh_s_max_edge,
                        n.clust,
                        m.clust,
                        sample.strat, ## 15
                        cores,
                        ndraws,
                        alphaj.pri,
                        nug.pri,
                        inla.int.strat, ## 20
                        inla.approx, 
                        Nsim,
                        data.lik,
                        sd.norm)

loopvars$run_date <- NA ## keep track of run_dates to later compare runs

## prepare a set of run_dates so we can write the complete loopvars to each run_date dir
## TODO for long loopvars this is stupidly slow... manually create the rundate list
for(ii in 1:nrow(loopvars)){
  loopvars$run_date[ii] <- make_time_stamp(TRUE)
  Sys.sleep(1.1)
}

for(ii in 1:nrow(loopvars)){

  ## make a run_date and setup output directory
  run_date <- loopvars$run_date[ii]

  ## now we can setup our main directory to save these results and log our note and stdouts
  out.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', run_date)
  dir.create(out.dir)
  dir.create(paste0(out.dir, '/errors'))
  dir.create(paste0(out.dir, '/output'))
  
  ## write the log note
  fileConn <- file(sprintf("%s/run_notes.txt", out.dir))
  writeLines(logging_note, fileConn)
  close(fileConn)

  ## save loopvars to this dir to reload into the parallel env
  saveRDS(file = paste0(out.dir, '/loopvars.rds'), obj = loopvars)

  ## ## save and reload loopvars in parallel env. that way, we only need to pass in iter/row #
  ## qsub.string <- qsub_sim(iter = ii, ## sets which loopvar to use in parallel
  ##                         run_date = run_date,
  ##                         slots = 4, 
  ##                         codepath = '/homes/azimmer/tmb_transition/realistic_sims/1_run_simulation.R', 
  ##                         singularity = 'default',
  ##                         singularity_opts = NULL,
  ##                         logloc = NULL ## defaults to input/output dir in sim run_date dir
  ##                         )

  ## ## launch the job
  
  ## system(qsub.string)
}
