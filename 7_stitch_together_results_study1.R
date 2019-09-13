## this script may be run in a clean env:
# ## Set core_repo location and tmb_repo loc
# user      <- Sys.info()['user']
# core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
# tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)
# 
# ## grab libraries and functions from MBG code
# setwd(core_repo)
# commondir    <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
# package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))
# 
# ## Load MBG packages and functions
# message('Loading in required R packages and MBG functions')
# source(paste0(core_repo, 'mbg_central/setup.R'))
# mbg_setup(package_list = package_list, repos = core_repo)
# 
# library(TMB)
# library(gridExtra)
# library(grid)
# library(RColorBrewer)
# library(viridis)

## this script pulls together some overall results from an experiment run
main.dir.names <- c('2019_09_12_08_06_46') ## study 1

## utility function from the interwebs
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
    ## New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm),
          l.ci = stats::quantile(xx[[col]], probs = (1-conf.interval)/2, na.rm=na.rm), 
          u.ci = stats::quantile(xx[[col]], probs = 1-(1-conf.interval)/2, na.rm=na.rm) 
        )
      },
      measurevar
    )

    ## Rename the "mean" column    
    datac <- plyr::rename(datac, c("mean" = measurevar,
                                   "l.ci.2.5%" = 'l.ci',
                                   "u.ci.97.5%" = 'u.ci'))
                                

    ## get with width of the ci
    datac$w.ci <- datac$u.ci - datac$l.ci

    return(datac)
}

for(main.dir.name in main.dir.names){

  message(sprintf('ON MAIN DIR %i of %i: %s', which(main.dir.names %in% main.dir.name), length(main.dir.names), main.dir.name))
  
  main.dir <- sprintf('/ihme/scratch/users/azimmer/tmb_inla_sim/%s', main.dir.name)
  
  compar.dir <- sprintf('%s/comparisons/', main.dir)
  dir.create(compar.dir, recursive = T, showWarnings = F)

  ## read in all experiment parameters
  loopvars <- fread(file = paste0(main.dir, '/loopvars.csv'), stringsAsFactors = F)
  
  # ## print columns of loopvar that vary
  # ## so we can easily see what's going on in the experiments that fail...
  loopvars[,!apply(loopvars, MARGIN = 2,
                   FUN=function(x){col.var <- sort(x, decreasing=F)[1] == sort(x, decreasing=T)[1]
                   if(is.na(col.var)){
                     return(TRUE)
                   }else{
                     return(col.var)}
                   }), with=F]
  
  ## read in the summary metrics log from each iteration of each experiment
  for(lvid in 1:nrow(loopvars)){
    message(sprintf('--loading in summary metrics from %i of %i', lvid, nrow(loopvars)))
    
    out.dir  <- sprintf('%s/%04d', main.dir, lvid)
    
    for(iter in 1:loopvars$n.sim[1]){ ## n.sim is the same for all experiments
      if(lvid==1 & iter==1){
        summary.metrics <- fread(sprintf('%s/validation/experiment%04d_iter%04d_summary_metrics.csv', 
                                            out.dir, lvid, iter))[,lvid:=lvid]
      }else{
        summary.metrics <- rbind(summary.metrics,
                                 fread(sprintf('%s/validation/experiment%04d_iter%04d_summary_metrics.csv', 
                                               out.dir, lvid, iter))[,lvid:=lvid], fill=T)
      }
    } ## loading all iterations within
  }   ## all experiments
  
  ## save the combined metrics for the study 
  write.csv(summary.metrics, file=sprintf('%sall_summary_metrics.csv', compar.dir))
  

  ## ####################################################################
  ## ####################################################################
  ## now that should be everything!
  ## we can make a bunch of plots
  ## ####################################################################
  ## ####################################################################
  cm.all <- summary.metrics
  
  ## make new labels for INLA_EB, INLA_CCD, TMB
  cm.all[inla.int.strat == 'eb' & mean.l.model == 'inla', fit_type := 'INLA_EB']
  cm.all[inla.int.strat == 'ccd' & mean.l.model == 'inla', fit_type := 'INLA_CCD']
  cm.all[mean.l.model == 'tmb', fit_type := 'TMB']

  library(tidyr)
  library(ggplot2)
  
  ## Gather columns into key-value pairs and move them from wide to long format
  long.cov <- data.table(gather(cm.all,
                                target_cov,  ## name of NEW key col
                                obs_cov,     ## name of NEW value col
                                cov25:cov95, ## cols to convert from wide to long
                                factor_key = TRUE
                                ))
  long.cov[target_cov == 'cov25', n_target_cov := 0.25]
  long.cov[target_cov == 'cov50', n_target_cov := 0.50]
  long.cov[target_cov == 'cov80', n_target_cov := 0.80]
  long.cov[target_cov == 'cov90', n_target_cov := 0.90]
  long.cov[target_cov == 'cov95', n_target_cov := 0.95]
  
  long.cov$noise_spatial_ratio <- long.cov$clust.var / long.cov$sp.var
  long.cov$noise_spatial_ratio[is.na(long.cov$noise_spatial_ratio)] <- 0

  long.cov[mean.l.model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
  long.cov[mean.l.model == 'tmb', pred.time :=  as.numeric(pred_time)]
  long.cov[mean.l.model == 'inla', fit.time :=  as.numeric(fit_time)]
  long.cov[mean.l.model == 'inla', pred.time :=  as.numeric(pred_time)]
  long.cov[, total.time :=  pred.time + fit.time]

  ## ########
  ## COVERAGE
  ## ########
  
  ## loop through binomial and normal with different obs variances
  ## make one plot for each different data model
  
  loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                            nv=c(NA, loopvars[, unique(norm.var)]))
  
  for(ii in 1:loopvars[,.N]){
    
    dl <- loop.params[ii, dl]
    nv <- loop.params[ii, nv]
    
    ## make the plot title
    pt <- sprintf('Average pixel coverage with %s observations', toupper(dl))
    if(dl == 'normal'){
      pt <- paste(pt, sprintf('\n normal SD = %.03f', nv))
    }
      
    ## TODO make plots with correct subsetting and titles
  }
  
  ## facet by observations
  long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
                             groupvars = c('n_target_cov', 'fit_type', 'n.clust', 'data.lik'))
  long.cov.sum$fit_n_lik <- paste(long.cov.sum$fit_type, long.cov.sum$data.lik, sep='_')
  
  fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                    aes(x = n_target_cov, y = obs_cov, shape=data.lik,
                                        colour = fit_n_lik, group = fit_n_lik),
                                    position=position_jitter(w=0.02, h=0.02)) + 
  geom_errorbar(aes(ymin = l.ci, ymax = u.ci), width = .025) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(. ~ n.clust) + ggtitle(pt)

  ggsave(sprintf('%s/%s_coverage_summary_nclust.png', compar.dir, loopvars$data.lik[1]),
         plot = fit_coverage_CI_summary,
         device = 'png', units = 'in',
         width = 12, height = 12)

  ## facet by noise to spatial signal
  long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
                             groupvars=c("n_target_cov","fit_type", 'noise_spatial_ratio'))
  
  fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                    aes(x=n_target_cov, y=obs_cov, colour=fit_type, group = fit_type)) + 
  geom_errorbar(aes(ymin=obs_cov-ci, ymax=obs_cov+ci), width=.025) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(. ~ noise_spatial_ratio) +
  ggtitle(sprintf('Comparison of coverage in: %s, faceted by clust.var/sp.var',
                  loopvars$data.lik[1]))

  ggsave(sprintf('%s/%s_coverage_summary_noise_to_spatial_var.png', compar.dir, loopvars$data.lik[1]),
         plot = fit_coverage_CI_summary,
         device = 'png', units = 'in',
         width = 12, height = 12)

  ## ##########
  ## TIME
  ## ##########

  
  cm.all[model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
  cm.all[model == 'tmb', pred.time :=  as.numeric(pred_time)]
  cm.all[model == 'inla', fit.time :=  as.numeric(fit_time)]
  cm.all[model == 'inla', pred.time :=  as.numeric(pred_time)]
  cm.all[, total.time :=  pred.time + fit.time]

  long.cov <- data.table(gather(cm.all,
                                operation,
                                time_s, 
                                fit.time:total.time,
                                factor_key = TRUE))

  
  ## facet by observations
  long.cov.sum <-  summarySE(long.cov, measurevar="time_s",
                             groupvars=c("operation","fit_type", 'n.clust'))
  fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                    aes(x=n.clust, y=time_s, colour=fit_type, group = fit_type)) + 
  geom_errorbar(aes(ymin=time_s-ci, ymax=time_s+ci), width=.01) +
  geom_line() +
  geom_point() +
  facet_wrap(. ~ operation) + ggtitle(sprintf('Comparison of fit time (sec) in: %s', loopvars$data.lik[1]))

  ggsave(sprintf('%s/%s_fit_time_nclust.png', compar.dir, loopvars$data.lik[1]),
         plot = fit_coverage_CI_summary,
         device = 'png', units = 'in',
         width = 12, height = 12)

  ## facet by noise to spatial signal
  long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
                             groupvars=c("n_target_cov","fit_type", 'noise_spatial_ratio'))
  
  fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                    aes(x=n_target_cov, y=obs_cov, colour=fit_type, group = fit_type)) + 
  geom_errorbar(aes(ymin=obs_cov-ci, ymax=obs_cov+ci), width=.025) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(. ~ noise_spatial_ratio) +
  ggtitle(sprintf('Comparison of coverage in: %s, faceted by clust.var/sp.var',
                  loopvars$data.lik[1]))

  ggsave(sprintf('%s/%s_coverage_summary_noise_to_spatial_var.png', compar.dir, loopvars$data.lik[1]),
         plot = fit_coverage_CI_summary,
         device = 'png', units = 'in',
         width = 12, height = 12)

}

##


## ## average my modeling tool
## if(i == 1){
##   cm.m.tmb  <- cm[mean.l.model == 'tmb,' lapply(.SD, mean), by=mean.l.model]
##   cm.m.inla <- cm[mean.l.model == 'inla', lapply(.SD, mean), by=mean.l.model]
## }else{
##   cm.m.tmb  <- rbind(cm.m.tmb,
##                      cm[mean.l.model == 'tmb,' lapply(.SD, mean), by=mean.l.model])
##   cm.m.inla <- rbind(cm.m.inla,
##                      cm[mean.l.model == 'inla', lapply(.SD, mean), by=mean.l.model])
## }
