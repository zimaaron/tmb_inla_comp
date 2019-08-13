## this script can be used to launch 1_run_simulation.R in parallel on the IHME cluster
## written by aoz
## dec 1 2018

## DO THIS!
################################################################################
## ADD A NOTE! to help identify what you were doing with this run
logging_note <- 'testing new parallelization on new cluster'

## make a master run_date to store all these runs in a single location
main.dir.name  <- NULL ## if NULL, run_date is made, OW uses name given
extra.job.name <- 'normal'
################################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## specify queue, project, and job requirements
q.q <- 'long.q'
q.m <- '20G' ## 5 gigs
q.t <- '00:10:00:00' ## DD:HH:MM:SS
q.p <- 0 ## priority: -1023 (low) - 0 (high)
cores <- 1 ## used for OMP/MKL in qsub_sim call TODO - parallel version?

#############################################
## setup the environment for singularity R ##
#############################################

## Set core_repo location and tmb_repo loc
user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)

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

## Now we can switch to the TMB repo
setwd(tmb_repo)
source('./realistic_sim_utils.R')
source('./qsub_utils.R')

## name overall directory with run_date if it is not named
if(is.null(main.dir.name)) main.dir.name <- make_time_stamp(TRUE)
print(main.dir.name)

## setup the main dir to store ALL experiments - i.e. each row in loopvar is an experiment
## all will be stored in main.dir, indexed by cooresponding row of loopvars 
main.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', main.dir.name)
dir.create(main.dir)

## write the log note
fileConn <- file(sprintf("%s/run_notes.txt", main.dir))
writeLines(logging_note, fileConn)
close(fileConn)

###############################
## setup things to loop over ##
###############################

## NOTES
## this list should align with the args in 1_run_simulation.R
## NULLs can't be passed this way, so NAs are stand-ins for NULL and this gets fixed in 1_run_space_sim.R
## NAs passed in can be used to turn things off (e.g. nug.var == NULL) removes the nugget
## 



## loopvars 1
reg <- 'nga'

## loopvars 
year_list <- 2000 

## loopvars 3
cov_names <- "c('access2', 'distrivers', 'evi'   , 'mapincidence')" 

## loopvars 4
cov_measures <- "c('mean'   , 'mean'      , 'median', 'mean')"

## loopvars 5
betas <- NA ## "c(.5, 1), 1, -.5)" ## IF THIS IS NA, NO COVS

## loopvars 6
alpha <- -1

## loopvars 7
sp.range <-  sqrt(8)     ## kappa=sqrt(8)/sp.range, so sp.range=sqrt(8) -> kappa=1 -> log(kappa)=0 (for R^2 domain)

## loopvars 8
sp.var <- 1.0 ^ 2        ## sp.var = 1/(4*pi*kappa^2*tau^2) (for R^2 domain)

## loopvars 9
sp.alpha <- 2.0          ## matern smoothness = sp.alpha - 1 (for R^2 domain)

## loopvars 10
nug.var <- NA# c(NA, .1 ^ 2, .1, .5 ^ 2, .5, 1) ##0.1 ^ 2       ## nugget variance

## loopvars 11
t.rho <-  0.8            ## annual temporal auto-corr

## loopvars 12
mesh_s_params <- c("c(0.1, 1, 5)") ## cutoff, largest allowed triangle edge length inner, and outer

## loopvars 13
n.clust <-  c(1000)         ## clusters PER TIME slice

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
cores <- 1

## loopvars 17
ndraws <- 500

## loopvars 18
alphaj.pri <- "c(0, 3)" ## normal mean and sd

## loopvars 19
nug.prec.pri <- "c(1, 1e-5)" ## gamma for nug precision with shape and inv-scale

## loopvars 20
inla.int.strat <- c('eb') ## can be 'eb', 'ccd', or 'grid'

## loopvars 21
inla.approx <- 'simplified.laplace' ## can be 'gaussian', 'simplified.laplace' (default) or 'laplace'

## loopvars 22
n.sim <- 3 ## number of times to repeat simulation

## loopvars 23
data.lik <- c('normal', 'binom') ## either 'binom' or 'normal'

## loopvars 24
norm.var <- c(.1^2)  ## sd of observations if normal. norm.var >= 0 and non-NA

## loopvars 25
norm.prec.pri <- "c(1, 1e-5)" ## gamma for normal obs  precision with shape and inv-scale

## loopvars 26
bias.correct <- c(TRUE) ## applies to both INLA and TMB!

## loopvars 27
sd.correct <- c(TRUE) ## only applies to TMB

## TODO always add all vars to exand.grid()
## NOTE: I use a named list here to ensure the columns in loopvars are named
loopvars <- expand.grid(list(reg = reg, ## 1
                             year_list = year_list,
                             cov_names = cov_names,
                             cov_measures = cov_measures,
                             betas = betas, ## 5
                             alpha = alpha,
                             sp.range = sp.range,
                             sp.var = sp.var,
                             sp.alpha = sp.alpha,
                             nug.var = nug.var, ## 10
                             t.rho = t.rho,
                             mesh_s_params = mesh_s_params,
                             n.clust = n.clust,
                             m.clust = m.clust,
                             sample.strat = sample.strat, ## 15
                             cores = cores,
                             ndraws = ndraws,
                             alphaj.pri = alphaj.pri,
                             nug.prec.pri = nug.prec.pri,
                             inla.int.strat = inla.int.strat, ## 20
                             inla.approx = inla.approx, 
                             n.sim = n.sim,
                             data.lik = data.lik,
                             norm.var = norm.var,
                             norm.prec.pri = norm.prec.pri, ## 25
                             bias.correct = bias.correct,
                             sd.correct = sd.correct))

## add a unique hash (for ease in tracking all iterations w/in an experiment) to each loopvar row
loopvars$qsub.hash <- stringi::stri_rand_strings(n=nrow(loopvars), length=5)

## save loopvars to this dir to reload into the parallel env
write.csv(file = paste0(main.dir, '/loopvars.csv'), x = loopvars,
          row.names = FALSE)

## make a data.table to save the job.ids
jid.dt <- data.table('exp'  = character(), 
                     'iter' = character(), 
                     'jid'  = character())

for(ll in 1:nrow(loopvars)){
  for(ii in 1:n.sim){
    
    message(sprintf('ON EXPERIMENT LOOP.ID %04d', ll))
    message(sprintf('----- on iter %04d', ii))
    
    if(ii == 1){
      ## now we can setup our main directory where all outputs will be saved
      dir.create(main.dir, showWarnings = F, recursive = TRUE)
      hold <- 0 ## to avoid ifelse() error
      hold.jid <- NULL
    }
    
    ## save and reload loopvars in parallel env. that way, we only need to pass in iter/row #
    qsub.string <- qsub_sim(exp.lvid = ll, ## sets which loopvar to use in parallel
                            exp.iter = ii,
                            exp.hash = loopvars$qsub.hash[ll],
                            main.dir = main.dir.name,
                            codepath = '/homes/azimmer/tmb_inla_comp/1_run_space_sim.R', 
                            singularity = 'default',
                            singularity_opts = NULL,
                            extra_name = extra.job.name,
                            mem = q.m,
                            time = q.t,
                            queue = q.q,
                            priority = q.p,
                            hold.jid = switch(hold + 1, NULL, hold.jid), ## NULL if hold==0, hold.jid if hold==1
                            logloc = NULL) ## defaults to input/output dir in sim run_date dir
    
    ## launch the job and catch the message
    sub.msg    <- system(qsub.string, intern=TRUE)
    print(sub.msg)
    
    ## save the job ids
    subbed.jid <- strsplit(sub.msg, split = ' ')[[1]][[3]]
    jid.dt <- rbind(jid.dt, list(sprintf('%04d', ll),
                                 sprintf('%04d', ii),
                                 subbed.jid))

    if(ii == 1){
      ## to use as holds for all others
      hold.jid <- subbed.jid
      
      ## all iters (except the 1st) hold on the 1st
      hold <- 1 
    }
  
  } ## iteration
}   ## experiment/loopvar row

## track job progress
in_q <- 1
while(in_q > 0) {
  jt <- track.exp.iter (jid.dt, main.dir)
  print(jt[['summ.tracker']])
  in_q <- sum(jt[['summ.tracker']]$in_q)
  Sys.sleep(60)
}

