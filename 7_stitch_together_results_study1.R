## this script may be run in a clean env:
# ## Set core_repo location and tmb_repo loc
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

## load util funcs
## Now we can switch to the TMB repo
source(paste0(tmb_repo,'/realistic_sim_utils.R'))

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(ggplot2)

## this script pulls together some overall results from an experiment run
# main.dir.name <- main.dir.names <- c('2019_12_18_23_25_42') ## study 1
# main.dir.name <- main.dir.names <- c('2019_09_29_21_38_14') ## study 2

## source('/homes/azimmer/tmb_inla_comp/7_stitch_together_results_study1.R')

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
  
  ## if we haven't already done this, 
  ## read in the summary metrics log from each iteration of each experiment
  if(!file.exists(sprintf('%sall_summary_metrics.csv', compar.dir))){
    for(lvid in 1:nrow(loopvars)){
      message(sprintf('--loading in summary metrics from %i of %i', lvid, nrow(loopvars)))
      
      out.dir  <- sprintf('%s/%04d', main.dir, lvid)
      
      for(iter in 1:loopvars$n.sim[1]){ ## n.sim is the same for all experiments
        if(lvid==1 & iter==1){
          summary.metrics <- fread(sprintf('%s/validation/experiment%04d_iter%04d_summary_metrics.csv', 
                                           out.dir, lvid, iter))[,lvid:=lvid]
          cov.gp   <- fread(sprintf('%s/validation/experiment%04d_iter%04d_GP_magnitude_coverage_summary.csv', 
                                           out.dir, lvid, iter))[,lvid:=lvid]
          cov.dist <- fread(sprintf('%s/validation/experiment%04d_iter%04d_distance_coverage_summary.csv', 
                                           out.dir, lvid, iter))[,lvid:=lvid]
        }else{
          summary.metrics <- rbind(summary.metrics,
                                   fread(sprintf('%s/validation/experiment%04d_iter%04d_summary_metrics.csv', 
                                                 out.dir, lvid, iter))[,lvid:=lvid], fill=T)
          cov.gp   <- rbind(cov.gp,
                            fread(sprintf('%s/validation/experiment%04d_iter%04d_GP_magnitude_coverage_summary.csv', 
                                          out.dir, lvid, iter))[,lvid:=lvid], fill = T)
          cov.dist <- rbind(cov.dist,
                            fread(sprintf('%s/validation/experiment%04d_iter%04d_distance_coverage_summary.csv', 
                                          out.dir, lvid, iter))[,lvid:=lvid], fill = T)
        }
      } ## loading all iterations within
    }   ## all experiments
    
    ## save the combined metrics for the study 
    write.csv(summary.metrics, file = sprintf('%sall_summary_metrics.csv', compar.dir))
    write.csv(cov.gp, file = sprintf('%sall_gp_coverage_metrics.csv', compar.dir))
    write.csv(cov.dist, file = sprintf('%sall_dist_coverage_metrics.csv', compar.dir))

    
    message('--summary metrics from all experiments and all iterations succesfully combined and saved')
 
  }else{
    
    ## reload the prepped file
    message('--reloading in combined summary metrics')
    summary.metrics <- fread(sprintf('%sall_summary_metrics.csv', compar.dir))
    cov.gp          <- fread(sprintf('%sall_gp_coverage_metrics.csv', compar.dir))
    cov.dist        <- fread(sprintf('%sall_dist_coverage_metrics.csv', compar.dir))

  }
  
  message('summary metrics from all experiments and all iterations prepped and ready to go')
  
  ## process a few things and make some necessary columns for plotting
  
  ## TODO fix a few things...
  ## TODO fix this elsewhere!
  setnames(summary.metrics, 'gauss_prec', 'gauss_prec_mean')
  setnames(summary.metrics, 'clust_prec', 'clust_prec_mean')
  
  ## get the true logkappa and logtau params
  summary.metrics[,logkappa := log(sqrt(8) / sp.range)]
  summary.metrics[,logtau   := log(sqrt(1 / (4 * pi * exp(logkappa) ^ 2 * sp.var)))]
  
  ## make new labels for INLA_EB, INLA_CCD, TMB
  summary.metrics[inla.int.strat == 'eb' & mean.l.model == 'inla', fit_type := 'INLA_EB']
  summary.metrics[inla.int.strat == 'ccd' & mean.l.model == 'inla', fit_type := 'INLA_CCD']
  summary.metrics[mean.l.model == 'tmb', fit_type := 'TMB']

  ## ####################################################################
  ## ####################################################################
  ## now that should be everything!
  ## we can make a bunch of plots
  ## after we format the data into long for ggplot
  ## ####################################################################
  ## ####################################################################
  
  ## Gather columns into key-value pairs and move them from wide to long format
  ## for plotting pixel level covereage
  long.cov <- data.table(gather(summary.metrics,
                                target_cov,  ## name of NEW key col
                                obs_cov,     ## name of NEW value col
                                pix.cov25:pix.cov95, ## cols to convert from wide to long
                                factor_key = TRUE
                                ))
  
  ## get the coverages calculated
  nom.cov.names <- long.cov[,sort(unique(target_cov))] 
  nom.cov.vals <- as.numeric(gsub("[^0-9]", "",  nom.cov.names))/100
  
  ## and assign the numeric coverage for the rows
  for(i in 1:length(nom.cov.names)){
    long.cov[target_cov == nom.cov.names[i], n_target_cov := nom.cov.vals[i]]
  }
  
  ## calculate noise to spatial ratio
  long.cov$noise_spatial_ratio <- long.cov$clust.var / long.cov$sp.var
  long.cov$noise_spatial_ratio[is.na(long.cov$noise_spatial_ratio)] <- 0

  ## process the total fitting and predict times
  long.cov[mean.l.model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
  long.cov[mean.l.model == 'tmb', pred.time :=  as.numeric(pred_time)]
  long.cov[mean.l.model == 'inla', fit.time :=  as.numeric(fit_time)]
  long.cov[mean.l.model == 'inla', pred.time :=  as.numeric(pred_time)]
  long.cov[, total.time :=  pred.time + fit.time]
  
  ## set cluster var NA to 0. NOTE: this is fit w/o cluster var params!
  long.cov[is.na(clust.var), clust.var:=0]

  ## ##############################################################
  ## PIXEL COVERAGE: facets by number of observations
  ## ##############################################################
  
  ## loop through binomial and normal with different obs variances
  ## make one plot for each different data model
  
  message('plotting average pixel coverage plots by sample size')
  
  loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                            nv=c(NA, loopvars[, unique(norm.var)]))
  
  cov.width <- 0.5 ## set the coverage interval width
  
  for(ii in 1:loop.params[,.N]){
    
    dl <- loop.params[ii, dl]
    nv <- loop.params[ii, nv]
    
    ## make the plot title
    pt <- sprintf('Average Pixel Coverage with %s Observations | %0.0f%% E.T. Intervals', 
                  stringr::str_to_title(dl), cov.width*100)
    if(dl == 'normal'){
      pt <- paste(pt, sprintf('\n Gaussian SD = %.03f', nv))
    }
      
    if(dl == 'binom') {
      long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
    }else if(dl == 'normal') {
      long.cov.sub <- subset(long.cov, data.lik == dl & norm.var == nv)
    }
    
    ## facet by observations
    long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                               groupvars = c('n_target_cov', 'fit_type', 'n.clust', 'data.lik', 'clust.var'),
                               conf.interval = cov.width)
    ## groups to plot lines
    long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$clust.var, sep='_')
    
    ## set facet names using the labeller
    ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
    facet_labs <- as_labeller(eval(parse(text=paste0('c(', paste0('`', sort(unique(long.cov$n.clust)), '`', '=', '\'', paste0('Num. Cluster Obs: ', sort(unique(long.cov$n.clust))), '\'', sep='', collapse=','), ')'))))
    
    pd <- position_dodge(0.05)
    fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                      aes(x = n_target_cov, y = obs_cov,
                                          shape = fit_type, 
                                          linetype = fit_type,
                                          color = factor(clust.var),
                                          group = line_group),
                                      position=position_jitter(w=0.02, h=0.02)) + 
      geom_errorbar(position=pd, aes(ymin = l.ci, ymax = u.ci), width = .025) +
      geom_line(position=pd) +
      geom_point(position=pd) +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(. ~ n.clust, labeller = facet_labs) + 
      ggtitle(pt) +
      ## fix the legends a bit
      labs(color = "Cluster Var", shape='Fit Type', linetype='Fit Type') + ## legend titles
      xlab('Nominal Coverage') + 
      ylab('Observed Monte Carlo Coverage')
    
    if(interactive()){ ## then we can view
      print(fit_coverage_CI_summary)
    }
    
    ggsave(sprintf('%s/%s_%0.2fnormVar_pixel_coverage_summary_Observations.png', compar.dir, dl, nv),
           plot = fit_coverage_CI_summary,
           device = 'png', units = 'in',
           width = 12, height = 8)
    
    Sys.sleep(10) ## to help with R studio crashing?
  }

  ## ##########################################################################
  ## PIXEL COVERAGE: facets by noise to spatial signal
  ## ##########################################################################
  
  message('plotting average pixel coverage plots by spatial noise:signal ratio')
  
  loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                            nv=c(NA, loopvars[, unique(norm.var)]))
  
  cov.width <- 0.5 ## set the coverage interval width
  for(ii in 1:loop.params[,.N]){
    
    dl <- loop.params[ii, dl]
    nv <- loop.params[ii, nv]
    
    ## make the plot title
    pt <- sprintf('Average Pixel Coverage with %s Observations | %0.0f%% E.T. Intervals', 
                  stringr::str_to_title(dl), cov.width*100)
    if(dl == 'normal'){
      pt <- paste(pt, sprintf('\n Gaussian SD = %.03f', nv))
    }
    
    if(dl == 'binom') {
      long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
    }else if(dl == 'normal') {
      long.cov.sub <- subset(long.cov, data.lik == dl & norm.var == nv)
    }
    
    ## facet by noise:spatial var ratio
    long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                               groupvars=c('n_target_cov', 'fit_type', 'n.clust', 'data.lik',
                                           'noise_spatial_ratio'))
    
    ## groups to plot lines
    long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')
    
    ## set facet names using the labeller
    ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
    clust.sp.pairs <- unique(long.cov, by=c("clust.var","sp.var"))[,.(clust.var, sp.var, noise_spatial_ratio)] 
    ## sort by unique noise_spatial_ratio order
    clust.sp.pairs <- clust.sp.pairs[order(noise_spatial_ratio),]
    facet_labs <- as_labeller(eval(parse(text = 
                    paste0('c(', paste0('`', sort(unique(long.cov$noise_spatial_ratio)), '`', '=', '\'', 
                                  sprintf('(Cluster Var)/(Spatial Var): %0.2f/%0.2f = %0.2f', 
                                          clust.sp.pairs$clust.var,
                                          clust.sp.pairs$sp.var,
                                          sort(unique(long.cov$noise_spatial_ratio))), '\'', sep='', collapse=','),
                           ')'))))
    
    pd <- position_dodge(0.05)
    fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                      aes(x=n_target_cov, y=obs_cov, 
                                          shape = fit_type,
                                          linetype = fit_type,
                                          colour=factor(n.clust), 
                                          group = line_group),
                                      position=position_jitter(w=0.01, h=0.01)) + 
      geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
      geom_line(position=pd) +
      geom_point(position=pd) +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(. ~ noise_spatial_ratio, labeller = facet_labs) +
      ggtitle(pt) +
      ## fix the legends a bit
      labs(color = 'Num. Obs', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
      xlab('Nominal Coverage') + 
      ylab('Observed Monte Carlo Coverage')
    if(interactive()){ ## then we can view
      print(fit_coverage_CI_summary)
    }
  
    ggsave(sprintf('%s/%s_%0.2fnormVar_pixel_coverage_summary_NoiseSpatialRatio.png', compar.dir, dl, nv),
           plot = fit_coverage_CI_summary,
           device = 'png', units = 'in',
           width = 12, height = 8)
  }
  
  
  
  ## ################################
  ## fixed effects bias & coverage ##
  ## ################################
  
  message('plotting fixed effect bias')
  
  ## get all the fixed effects
  non.pix.eff.sd.colnames <- names(summary.metrics)[grep('_sd$', names(summary.metrics))]
  ## now, get the first part of the name
  non.pix.eff.b.colnames <- paste0(unlist(lapply( non.pix.eff.sd.colnames, 
                                            function(x){
                                              substr(x, 1, nchar(x)-3)
                                            })), '_bias')
  
  summary.metrics[,fe_int_bias := fe_int_mean - alpha]
  summary.metrics[,gauss_sd_bias := sqrt(1/gauss_prec_mean) - sqrt(norm.var)]
  summary.metrics[,matern_logtau_bias := matern_logtau_mean - logtau]
  summary.metrics[,matern_logkappa_bias := matern_logkappa_mean - logkappa]
  summary.metrics[,clust_sd_bias := sqrt(1/clust_prec_mean) - sqrt(clust.var)]
  
  ## get some coverage probs
  key.summ.metrics <- summary.metrics[,c(non.pix.eff.b.colnames, 
                                         'data.lik', 'norm.var', 'n.clust', 'fit_type', 'clust.var'), with=F]
  fe.mean.long <- melt(key.summ.metrics, id=c('data.lik', 'norm.var', 'n.clust', 'fit_type', 'clust.var'))
  
  ## drop NAs. eg there is no gauss_prec_bias in a binom model, 
  ## and there is no clust_prec_bias when clust.var==NA
  fe.mean.long <- fe.mean.long[!is.na(value),]
  
  cov.width <- 0.5 ## set the coverage interval width
  for(ii in 1:loop.params[,.N]){
    
    dl <- loop.params[ii, dl]
    nv <- loop.params[ii, nv]
    
    ## make the plot title
    pt <- sprintf('Average Param Coverage with %s Observations | %0.0f%% E.T. Intervals', 
                  stringr::str_to_title(dl), cov.width*100)
    if(dl == 'normal'){
      pt <- paste(pt, sprintf('\n Gaussian SD = %.03f', nv))
    }
    
    if(dl == 'binom') {
      fe.mean.long.sub <- subset(fe.mean.long, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
    }else if(dl == 'normal') {
      fe.mean.long.sub <- subset(fe.mean.long, data.lik == dl & norm.var == nv)
    }
    
    ## fe.mean.long.sub <- subset(fe.mean.long.sub, fit_type=='INLA_EB')
    
    ## facet by ???
    fe.mean.long.sub <-  summarySE(fe.mean.long.sub, measurevar="value",
                               groupvars=c("data.lik",'norm.var', 'n.clust', "fit_type", 'clust.var', 'variable'),
                               conf.interval = 0.5)
    
    ## groups to plot lines
    ## long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')
    
    pd <- position_dodge(.25)
    fe_bias_summary <- ggplot(fe.mean.long.sub,
                              aes(x=log(n.clust), y=med, 
                                  shape = fit_type,
                                  linetype = fit_type,
                                  colour=factor(clust.var))) + 
      geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
      geom_line(position=pd) +
      geom_point(position=pd, size=2) +
      geom_abline(intercept = 0, slope=0) +
      facet_wrap(. ~ variable, scales='free_y') +
      ggtitle(pt) +
      ## fix the legends a bit
      labs(color = 'Clust. Var', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
      xlab('Num Spatial Locations Sampled') + 
      scale_x_continuous(breaks = log(sort(fe.mean.long.sub$n.clust)), ## add log x-axis labels
                         labels = paste0('ln(', sort(fe.mean.long.sub$n.clust), ')')) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 1)) + ## rotate x-axis ticks
      ylab('Bias: Est-True')
    
    if(interactive()){ ## then we can view
      print(fe_bias_summary)
    }
    
     ggsave(sprintf('%s/%s_%0.2fnormVar_param_bias_summary.png', compar.dir, dl, nv),
           plot = fe_bias_summary,
           device = 'png', units = 'in',
           width = 12, height = 8)
  }
  
  ################
  ## time plots v1 - more below ## 
  ################

  ## get the coverages calculated
  nom.cov.names <- long.cov[,sort(unique(target_cov))] 
  nom.cov.vals <- as.numeric(gsub("[^0-9]", "",  nom.cov.names))/100
  
  ## and assign the numeric coverage for the rows
  for(i in 1:length(nom.cov.names)){
    long.cov[target_cov == nom.cov.names[i], n_target_cov := nom.cov.vals[i]]
  }
  
  ## calculate noise to spatial ratio
  long.cov$noise_spatial_ratio <- long.cov$clust.var / long.cov$sp.var
  long.cov$noise_spatial_ratio[is.na(long.cov$noise_spatial_ratio)] <- 0
  
  ## process the total fitting and predict times
  long.cov[mean.l.model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
  long.cov[mean.l.model == 'tmb', pred.time :=  as.numeric(pred_time)]
  long.cov[mean.l.model == 'inla', fit.time :=  as.numeric(fit_time)]
  long.cov[mean.l.model == 'inla', pred.time :=  as.numeric(pred_time)]
  long.cov[, total.time :=  pred.time + fit.time]
  
  ## set cluster var NA to 0. NOTE: this is fit w/o cluster var params!
  long.cov[is.na(clust.var), clust.var:=0]
  
  
  loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                            nv=c(NA, loopvars[, unique(norm.var)]))
  
  cov.width <- 0.5 ## set the coverage interval width
  for(ii in 1:loop.params[,.N]){
    
    dl <- loop.params[ii, dl]
    nv <- loop.params[ii, nv]
    
    ## make the plot title
    pt <- sprintf('Average Pixel Coverage with %s Observations | %0.0f%% E.T. Intervals', 
                  stringr::str_to_title(dl), cov.width*100)
    if(dl == 'normal'){
      pt <- paste(pt, sprintf('\n Gaussian SD = %.03f', nv))
    }
    
    if(dl == 'binom') {
      long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
    }else if(dl == 'normal') {
      long.cov.sub <- subset(long.cov, data.lik == dl & norm.var == nv)
    }
    
    ## facet by noise:spatial var ratio
    long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                               groupvars=c('n_target_cov', 'fit_type', 'n.clust', 'data.lik',
                                           'noise_spatial_ratio'))
    
    ## groups to plot lines
    long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')
    
    ## set facet names using the labeller
    ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
    clust.sp.pairs <- unique(long.cov, by=c("clust.var","sp.var"))[,.(clust.var, sp.var, noise_spatial_ratio)] 
    ## sort by unique noise_spatial_ratio order
    clust.sp.pairs <- clust.sp.pairs[order(noise_spatial_ratio),]
    facet_labs <- as_labeller(eval(parse(text = 
                                           paste0('c(', paste0('`', sort(unique(long.cov$noise_spatial_ratio)), '`', '=', '\'', 
                                                               sprintf('(Cluster Var)/(Spatial Var): %0.2f/%0.2f = %0.2f', 
                                                                       clust.sp.pairs$clust.var,
                                                                       clust.sp.pairs$sp.var,
                                                                       sort(unique(long.cov$noise_spatial_ratio))), '\'', sep='', collapse=','),
                                                  ')'))))
    
    pd <- position_dodge(0.05)
    fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                      aes(x=n_target_cov, y=obs_cov, 
                                          shape = fit_type,
                                          linetype = fit_type,
                                          colour=factor(n.clust), 
                                          group = line_group),
                                      position=position_jitter(w=0.01, h=0.01)) + 
      geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
      geom_line(position=pd) +
      geom_point(position=pd) +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(. ~ noise_spatial_ratio, labeller = facet_labs) +
      ggtitle(pt) +
      ## fix the legends a bit
      labs(color = 'Num. Obs', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
      xlab('Nominal Coverage') + 
      ylab('Observed Monte Carlo Coverage')
    if(interactive()){ ## then we can view
      print(fit_coverage_CI_summary)
    }
    
    ggsave(sprintf('%s/%s_%0.2fnormVar_pixel_coverage_summary_NoiseSpatialRatio.png', compar.dir, dl, nv),
           plot = fit_coverage_CI_summary,
           device = 'png', units = 'in',
           width = 12, height = 8)
  }
  
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
  
  if(interactive()){
    ## then we can view
    print(fit_coverage_CI_summary)
  }

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

## also need to get bias and coverage for 
## get all params by looking for things that end with _sd
non.pix.eff.sd.colnames <- names(summary.metrics)[grep('_sd$', names(summary.metrics))]
## now, get the first part of the name
non.pix.eff.sd.colnames <- unlist(lapply( non.pix.eff.sd.colnames, 
                                          function(x){
                                            substr(x, 1, nchar(x)-3)
                                          }))
#### get some coverage probs
coverage_probs <- c(25, 50, 80, 90, 95)
## it's much faster to get all quantiles at once then to iteratively calculate them
lui <- apply(all.preds, 1, quantile, 
             p = c((1 - coverage_probs/100)/2, 
                   coverage_probs/100 + (1 - coverage_probs/100) / 2), na.rm=T)
for(cc in 1:length(coverage_probs)){
  c <- coverage_probs[cc]
  message(paste0('-------- calcing ', c ,'% coverage of binomial prob'))
  li       <- lui[cc,]
  ui       <- lui[length(coverage_probs) + cc]
  d[, paste0('pix.cov.',c)] <- d[['truth']] >= li & d[['truth']] <= ui
}




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
