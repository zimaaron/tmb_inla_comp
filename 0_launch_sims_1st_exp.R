## this script can be used to launch 1_run_simulation.R in parallel on the IHME cluster
## written by aoz
## 2020JAN03
## source('/homes/azimmer/tmb_inla_comp/0_launch_sims_1st_exp.R')

## DO THIS!
################################################################################
## ADD A NOTE! to help identify what you were doing with this run
logging_note <- 
'STUDY 01: vary number of clusters, cluster effect, and normal data variance. 
TRIAL 22: testing new pc.prior for matern'

## make a master run_date to store all these runs in a single location
main.dir.name  <- NULL ## IF NULL, run_date is made, OW uses name given
extra.job.name <- 'study01trial22'
################################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## specify queue, project, and job requirements
q.q   <- 'geospatial.q' ## all.q ## long.q
q.m   <- '25G' ## e.g. 10G
q.t   <- '00:3:30:00' ## DD:HH:MM:SS
q.p   <- -100 ## priority: -1023 (low) - 0 (high)
cores <- 1 ## used for OMP/MKL in qsub_sim call TODO - parallel version? NOTE!! this is overwritten in arg 16

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
main.dir  <- sprintf('/ihme/scratch/users/azimmer/tmb_inla_sim/%s', main.dir.name)
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
## NAs passed in can be used to turn things off (e.g. clust.var == NULL) removes the cluster RE


## loopvars 1: region to model over
reg <- 'nga'

## loopvars: year to use (mostly for covariates)
year_list <- 2000 

## loopvars 3: vector of covariates to use in model
cov_names <- "c('access2')" 

## loopvars 4: vector of covariate measures to use in conjunction with cov_names to load covs
cov_measures <- "c('mean')"

## loopvars 5: cov effects. either NA (NO Covs), or a vector of cov effects to use wtih cov_names 
betas <- NA

## loopvars 6 ## global intercept
alpha <- -1.0

## loopvars 7: spatial range as defined by INLA folks
## units are in degrees lat-long!
## Nigeria is approx 12 degrees wide and 10 degrees tall
## kappa=sqrt(8)/sp.range, so sp.range=sqrt(8) -> kappa=1 -> log(kappa)=0 (for R^2 domain)
## so, with kappa=1, 90% of the correlation drops by 2.8 degrees, or about 1/4 of the heigth/width 
sp.range <- sqrt(8)     

## loopvars 8: spatial nominal field variance as defiend by INLA folks
## sp.var = 1/(4*pi*kappa^2*tau^2) (for R^2 domain)
sp.var <- 0.5 ^ 2        

## loopvars 9: matern smoothness = sp.alpha - 1 -> sp.alpha = matern smooth + 1 (for R^2 domain)
sp.alpha <- 2.0          

## loopvars 10: cluster RE variance. NA means no effect
clust.var <-  c(NA, (c(1, 2, 4) / 10) ^ 2) 

## loopvars 11: temporal auto-correlation (NOT USED IN SPACE-ONLY MODEL)
t.rho <-  0.8           

## loopvars 12: R2 mesh args: largest allowed triangle edge length inner, and outer
mesh_s_params <- c("c(0.4, 5)") 

## loopvars 13: number of clusters to simulate per year
n.clust <- c(250, 500, 750, 1000, 2500, 5000)

## loopvars 14: mean number of individuals sim'ed per cluster using poisson(m.clust)
m.clust <- 35                   

## loopvars 15
## each entry must be a character string with the syntax for a 3 element R list containing:
## 1) obs.loc.strat: (either 'rand' or 'pop.strat')
## 2) urban.pop.pct:   a number between 0 and 100. the % of population that belongs to urban pixels
## 3) urban.strat.pct: a number between 0 and 100. the % of observations that should come fom urban pixels
sample.strat <- "list(obs.loc.strat='rand',
                      urban.pop.pct=5,
                      urban.strat.pct=40)"  ## random or by population for now. ## TODO cluster design

## loopvars 16: cores to use in laun
cores <- 1

## loopvars 17: number of fitted model draws to take
ndraws <- 500

## loopvars 18: mean and sd for normal prior on fixed effects (alpha and betas)
alphaj.pri <- "c(0, 3)" ## N(mean, sd)

## loopvars 19: pc.prior on clust RE precision
## (u, a) s.t. P(1/sqrt(prec) > u) = a, i.e. P(SD > u) = a
clust.prec.pri <- "c(.5, .05)" 

## loopvars 20: INLA hyperparam integration strategy. can be 'eb', 'ccd', or 'grid'
inla.int.strat <- c('eb')

## loopvars 21: INLA marginal posterior approx strategy: can be 'gaussian', 'simplified.laplace' (default) or 'laplace'
inla.approx <- 'simplified.laplace' 

## loopvars 22: number of times to repeat an experiment (monte carlo simulations)
n.sim <- 100

## loopvars 23: data distribution: either 'binom' or 'normal'
data.lik <- c('normal', 'binom') 

## loopvars 24: ONLY FOR data.lik=='normal'. variance of INDIVIDUAL normal data obs.
norm.var <- (c(1, 2, 4, 5) / 10) ^ 2

## loopvars 25: pc.prior on normal individual level precision
## (u, a) s.t. P(1/sqrt(prec) > u) = a, i.e. P(SD > u) = a
norm.prec.pri <- "c(1, .01)"

## loopvars 26: bias correct the mean estimates. NOTE: applies to both INLA and TMB!!
bias.correct <- c(TRUE) 

## loopvars 27: perform sd correction. NOTE!! for TMB only
sd.correct <- c(TRUE)

## loopvars 28: pc.prior on spde parameters
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
matern.pri <- "c(10, .95, 1., .05)" ## a, b, c, d

## TODO always add all vars to exand.grid()
## NOTE: I use a named list here to ensure the columns in loopvars are named
loopvars <- data.table(expand.grid(list(reg = reg, ## 1
                             year_list = year_list,
                             cov_names = cov_names,
                             cov_measures = cov_measures,
                             betas = betas, ## 5
                             alpha = alpha,
                             sp.range = sp.range,
                             sp.var = sp.var,
                             sp.alpha = sp.alpha,
                             clust.var = clust.var, ## 10
                             t.rho = t.rho,
                             mesh_s_params = mesh_s_params,
                             n.clust = n.clust,
                             m.clust = m.clust,
                             sample.strat = sample.strat, ## 15
                             cores = cores,
                             ndraws = ndraws,
                             alphaj.pri = alphaj.pri,
                             clust.prec.pri = clust.prec.pri,
                             inla.int.strat = inla.int.strat, ## 20
                             inla.approx = inla.approx, 
                             n.sim = n.sim,
                             data.lik = data.lik,
                             norm.var = norm.var,
                             norm.prec.pri = norm.prec.pri, ## 25
                             bias.correct = bias.correct,
                             sd.correct = sd.correct,
                             matern.pri)))

## drop wasteful combinations that don't need to be run

## drop varying norm.var with binomial
loopvars <- loopvars[!(data.lik=='binom' & norm.var > min(norm.var)),]

message(sprintf('YOU ARE ABOUT TO LAUNCH %i EXPERIMENTS', nrow(loopvars)))
message(sprintf('-- EACH WITH %i ITERATIONS', n.sim))
message(sprintf("---- THAT'S %i JOBS!", n.sim*nrow(loopvars)))

## add a unique hash (for ease in tracking all iterations w/in an experiment) to each loopvar row
loopvars$qsub.hash <- stringi::stri_rand_strings(n=nrow(loopvars), length=5)

## save loopvars to this dir to reload into the parallel env
write.table(file = paste0(main.dir, '/loopvars.csv'), x = loopvars,
            row.names = FALSE, sep=',')

## make a data.table to save the job.ids
jid.dt <- data.table('exp'  = character(), 
                     'iter' = character(), 
                     'jid'  = character(),
                     'qsub' = character())

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
                                 subbed.jid,
                                 qsub.string))

    if(ii == 1){
      ## to use as holds for all others
      hold.jid <- subbed.jid
      
      ## all iters (except the 1st) hold on the 1st
      hold <- 1 
    }
  
  } ## iteration
}   ## experiment/loopvar row

time.start.running <- proc.time()

## track job progress
in_q <- 1
while(in_q > 0) {
  message(paste0('\n\n', Sys.time()))
  jt1 <- track.exp.iter(jid.dt, main.dir)
  js1 <- jt1[['summ.tracker']]
  jf1 <- jt1[['full.tracker']]
  print(js1, nrow(js1))
  print(jf1[, c(on_track_per=mean(on_track)*100,
          in_q=sum(in_q),
          running=sum(running),
          errored=sum(errored),
          completed=sum(completed),
          completed_per=mean(completed)*100)],
          digits=3)
  in_q <- sum(js1$in_q)
  Sys.sleep(60)
}

message(sprintf('%.2f%% of your experiments completed successfully', 
                mean(js1[, completed]==100)*100))
message(sprintf('%.2f%% of your iterations across all experiments completed successfully', 
                mean(jt1[['full.tracker']][, errored]==0)*100))
message('These experiments had some iterations fail:')
print(js1[errored > 0,])
message('These are the failed experiment iterations:')
print(jt1[['full.tracker']][errored==1, .(exp, iter, jid)])

## save the environment - mostly the job ideas and qsub commands
save(list=ls(), file = sprintf('%s/completed_env.rdata', main.dir))


## relaunch the jobs that failed to see if that helps
if(mean(jt1[['full.tracker']][, errored]==0) != 1){
  
  ## relaunch failed jobs - IF we think this will help
  ## load(sprintf('%s/completed_env.rdata', main.dir))
  failed.j <- jt1[['full.tracker']][errored==1, .(exp, iter)]
  for(ff in 1:nrow(failed.j)){
    ll <- as.numeric(failed.j[ff, exp])
    ii <- as.numeric(failed.j[ff, iter])
    message(sprintf('ON EXPERIMENT LOOP.ID %04d', ll))
    message(sprintf('----- on iter %04d', ii))
    
    if(ii == 1){
      hold <- 0
      hold.jid <- NULL
    }else{
      hold <- 1
      hold.jid <- jid.dt[exp==failed.j[ff, exp] & iter=='0001', jid]
    }
    
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
    
    ## replace the job id and the new qsub string
    subbed.jid <- strsplit(sub.msg, split = ' ')[[1]][[3]]
    jid.dt[exp==failed.j[ff, exp] & iter==failed.j[ff, iter], 
           `:=`(jid=subbed.jid, qsub=qsub.string)]
  }
  
  ## track job progress
  in_q <- 1
  while(in_q > 0) {
    message(paste0('\n\n', Sys.time()))
    jt2 <- track.exp.iter(jid.dt, main.dir)
    js2 <- jt2[['summ.tracker']]
    jf2 <- jt2[['full.tracker']]
    print(js2, nrow(js2))
    print(jf2[, c(on_track_per=mean(on_track)*100,
                  in_q=sum(in_q),
                  running=sum(running),
                  errored=sum(errored),
                  completed=sum(completed),
                  completed_per=mean(completed)*100)],
          digits=3)
    in_q <- sum(js2$in_q)
    Sys.sleep(60)
  }
  
  
  message(sprintf('%.2f%% of your experiments completed successfully', 
                  mean(js2[, completed]==100)*100))
  message(sprintf('%.2f%% of your iterations across all experiments completed successfully', 
                  mean(jt2[['full.tracker']][, errored]==0)*100))
  # message('These experiments had some iterations fail:')
  # print(js2[errored > 0,])
  message('These are the failed experiment iterations:')
  print(jt2[['full.tracker']][errored==1, .(exp, iter, jid)])
  
  ## save the environment - mostly the job ideas and qsub commands
  save(list=ls(), file = sprintf('%s/completed_env_2.rdata', main.dir))
}

time.done.running <- proc.time()

message('running time:')
print(time.done.running - time.start.running)

# ## print columns of loopvar that vary
# ## so we can easily see what's going on in the experiments that fail...
# loopvars[jt2[['full.tracker']][errored==1, as.numeric(unique(exp))],
#          !apply(loopvars, MARGIN = 2,
#                  FUN=function(x){col.var <- sort(x, decreasing=F)[1] == sort(x, decreasing=T)[1]
#                  if(is.na(col.var)){
#                    return(TRUE)
#                  }else{
#                    return(col.var)}
#                  }), with=F]
