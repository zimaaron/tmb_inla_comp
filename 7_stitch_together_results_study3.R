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
setwd(tmb_repo)
source('./realistic_sim_utils.R')

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(ggplot2)

## this script pulls together some overall results from an experiment run
## main.dir.name <- main.dir.names <- c('2020_01_07_08_13_43') ## study 1
# main.dir.name <- main.dir.names <- c('2019_09_29_21_38_14') ## study 2

## source('/homes/azimmer/tmb_inla_comp/7_stitch_together_results_study1.R')

## message(sprintf('ON MAIN DIR %i of %i: %s', which(main.dir.names %in% main.dir.name), length(main.dir.names), main.dir.name))

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
missing.files <- data.table(lvid=integer(),
                            iter=integer())
if(!file.exists(sprintf('%sall_summary_metrics.csv', compar.dir))){
  for(lvid in 1:nrow(loopvars)){
    message(sprintf('--loading in summary metrics from %i of %i', lvid, nrow(loopvars)))
    
    out.dir  <- sprintf('%s/%06d', main.dir, lvid)
    
    for(iter in 1:loopvars$n.sim[1]){ ## n.sim is the same for all experiments
      
      ## check to see that the file is there
      if(!file.exists(sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv', 
                              out.dir, lvid, iter))){
        message(sprintf('----WARNING!! summary file for exp %06d, iter %06d does not exist!', lvid, iter))
        missing.files <- rbind(missing.files, 
                               data.table(lvid=lvid,
                                          iter=iter))
      }else{
        if(lvid==1 & iter==1){
          summary.metrics <- fread(sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv', 
                                           out.dir, lvid, iter))[,lvid:=lvid]
          cov.gp   <- fread(sprintf('%s/validation/experiment%06d_iter%06d_GP_magnitude_coverage_summary.csv', 
                                    out.dir, lvid, iter))[,lvid:=lvid]
          cov.dist <- fread(sprintf('%s/validation/experiment%06d_iter%06d_distance_coverage_summary.csv', 
                                    out.dir, lvid, iter))[,lvid:=lvid]
        }else{
          summary.metrics <- rbind(summary.metrics,
                                   fread(sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv', 
                                                 out.dir, lvid, iter))[,lvid:=lvid], fill=T)
          cov.gp   <- rbind(cov.gp,
                            fread(sprintf('%s/validation/experiment%06d_iter%06d_GP_magnitude_coverage_summary.csv', 
                                          out.dir, lvid, iter))[,lvid:=lvid], fill = T)
          cov.dist <- rbind(cov.dist,
                            fread(sprintf('%s/validation/experiment%06d_iter%06d_distance_coverage_summary.csv', 
                                          out.dir, lvid, iter))[,lvid:=lvid], fill = T)
        }
      } ## file.exists check
    }   ## loading all iterations within
  }     ## all experiments
  
  ## save the combined metrics for the study 
  write.csv(summary.metrics, file = sprintf('%sall_summary_metrics.csv', compar.dir))
  write.csv(cov.gp, file = sprintf('%sall_gp_coverage_metrics.csv', compar.dir))
  write.csv(cov.dist, file = sprintf('%sall_dist_coverage_metrics.csv', compar.dir))
  
  if(nrow(missing.files) > 0){
    message('-- WARNING! not all experiments and iters had summary metrics files. missing:')
    print(missing.files)
  }else{
    message('--summary metrics from all experiments and all iterations succesfully combined and saved')  
  }
  
}else{
  
  ## reload the prepped file
  message('--reloading in combined summary metrics')
  summary.metrics <- fread(sprintf('%sall_summary_metrics.csv', compar.dir))
  cov.gp          <- fread(sprintf('%sall_gp_coverage_metrics.csv', compar.dir))
  cov.dist        <- fread(sprintf('%sall_dist_coverage_metrics.csv', compar.dir))
  
}

message('summary metrics from all experiments and all iterations prepped and ready to go')

## process a few things and make some necessary columns for plotting

## set the order of clust.var so that NA -> None and comes first
## set factor order of clust.var
clust.var <- as.character(summary.metrics$clust.var)
non.na.cv <- as.character(sort(unique(na.omit(summary.metrics$clust.var))))
clust.var[is.na(clust.var)] <- 'None'
summary.metrics$clust.var.cat <- factor(clust.var, levels = c('None', non.na.cv))

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

## ##########################################################################
## ##########################################################################
## PIXEL COVERAGE: facets by number of observations
## ##########################################################################
## ##########################################################################

## loop through binomial and normal with different obs variances
## make one plot for each different data model

message('plotting average pixel coverage plots by sample size')

loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                          nv=c(NA, loopvars[, unique(norm.var)]))

cov.width <- 0.5 ## set the coverage interval width

for(dl in c('binom', 'normal')){
  
  ## make the plot title
  pt <- sprintf('Average Pixel Coverage given %s Observations', 
                stringr::str_to_title(dl))
  pt <- paste(pt, sprintf('\n %0.0f%% Equal Tail Intervals', cov.width*100))
  
  if(dl == 'binom') {
    long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
  }else if(dl == 'normal') {
    long.cov.sub <- subset(long.cov, data.lik == dl)
  }
  
  ## facet by observations
  long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                             groupvars = c('n_target_cov', 'fit_type', 'n.clust', 'data.lik', 'clust.var.cat', 'norm.var'),
                             conf.interval = cov.width)
  ## groups to plot lines
  long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$clust.var, sep='_')
  
  ## set facet names using the labeller
  ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
  ## facet_labs <- as_labeller(eval(parse(text=paste0('c(', paste0('`', sort(unique(long.cov$n.clust)), '`', '=', '\'', paste0('Num. Cluster Obs: ', sort(unique(long.cov$n.clust))), '\'', sep='', collapse=','), ')'))))
  
  ## rename for plotting
  setnames(long.cov.sum, 'n.clust', 'num.obs')
  setnames(long.cov.sum, 'norm.var', 'obs.var')
  
  pd <- position_dodge(0.05)
  fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                    aes(x = n_target_cov, y = obs_cov,
                                        shape = fit_type, 
                                        linetype = fit_type,
                                        color = clust.var.cat,
                                        group = line_group),
                                    position=position_jitter(w=0.02, h=0.02)) + 
    geom_errorbar(position=pd, aes(ymin = l.ci, ymax = u.ci), width = .025) +
    geom_line(position=pd) +
    geom_point(position=pd) +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle(pt) +
    ## fix the legends a bit
    labs(color = "Cluster Var", shape='Fit Type', linetype='Fit Type') + ## legend titles
    xlab('Nominal Coverage') + 
    ylab('Observed Monte Carlo Coverage')
  
  if(dl=='binom'){
    fit_coverage_CI_summary <- fit_coverage_CI_summary +
      facet_wrap(. ~ num.obs, labeller = label_both)
  }
  if(dl == 'normal'){
    fit_coverage_CI_summary <- fit_coverage_CI_summary +
      facet_grid(obs.var ~ num.obs, labeller = label_both)
  }
  
  if(interactive()){ ## then we can view
    print(fit_coverage_CI_summary)
  }
  
  ggsave(sprintf('%s/%s_plot_01_pixel_coverage_summary_by_observations.png', compar.dir, dl),
         plot = fit_coverage_CI_summary,
         device = 'png', units = 'in',
         width = 12, height = 8)
}

## ##########################################################################
## ##########################################################################
## PIXEL COVERAGE: facets by noise to spatial signal
## ##########################################################################
## ##########################################################################

message('plotting average pixel coverage plots by spatial noise:signal ratio')

loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                          nv=c(NA, loopvars[, unique(norm.var)]))

cov.width <- 0.5 ## set the coverage interval width
for(dl in c('binom', 'normal')){
  
  ## make the plot title
  pt <- sprintf('Average Pixel Coverage given %s Observations', 
                stringr::str_to_title(dl))
  pt <- paste(pt, sprintf('\n %0.0f%% Equal Tail Intervals', cov.width*100))
  
  if(dl == 'binom') {
    long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
  }else if(dl == 'normal') {
    long.cov.sub <- subset(long.cov, data.lik == dl)
  }
  
  ## facet by noise:spatial var ratio
  long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                             groupvars=c('n_target_cov', 'fit_type', 'n.clust', 'data.lik',
                                         'noise_spatial_ratio', 'norm.var'))
  
  ## groups to plot lines
  long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')
  
  ## rename for plotting
  setnames(long.cov.sum, 'noise_spatial_ratio', 'clust.var_over_sp.var')
  setnames(long.cov.sum, 'norm.var', 'obs.var')
  

  # ## set facet names using the labeller
  # ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
  # clust.sp.pairs <- unique(long.cov, by=c("clust.var","sp.var"))[,.(clust.var, sp.var, noise_spatial_ratio)] 
  # ## sort by unique noise_spatial_ratio order
  # clust.sp.pairs <- clust.sp.pairs[order(noise_spatial_ratio),]
  # facet_labs <- as_labeller(eval(parse(text = 
  #                                        paste0('c(', paste0('`', sort(unique(long.cov$noise_spatial_ratio)), '`', '=', '\'', 
  #                                                            sprintf('(Cluster Var)/(Spatial Var): %0.2f/%0.2f = %0.2f', 
  #                                                                    clust.sp.pairs$clust.var,
  #                                                                    clust.sp.pairs$sp.var,
  #                                                                    sort(unique(long.cov$noise_spatial_ratio))), '\'', sep='', collapse=','),
  #                                               ')'))))
  
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
    ## facet_wrap(. ~ noise_spatial_ratio, labeller = facet_labs) +
    ggtitle(pt) +
    ## fix the legends a bit
    labs(color = 'Num. Obs', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
    xlab('Nominal Coverage') + 
    ylab('Observed Monte Carlo Coverage')
  
  if(dl=='binom'){
    fit_coverage_CI_summary <- fit_coverage_CI_summary +
      facet_wrap(. ~ clust.var_over_sp.var, labeller = label_both)
  }
  if(dl == 'normal'){
    fit_coverage_CI_summary <- fit_coverage_CI_summary +
      facet_grid(clust.var_over_sp.var ~ obs.var, labeller = label_both)
  }
  
  if(interactive()){ ## then we can view
    print(fit_coverage_CI_summary)
  }
  
  ggsave(sprintf('%s/%s_plot_02_pixel_coverage_summary_by_noiseSpatialRatio.png', compar.dir, dl),
         plot = fit_coverage_CI_summary,
         device = 'png', units = 'in',
         width = 12, height = 8)
}


## ##########################################################################
## ##########################################################################
## effects bias & coverage 
## ##########################################################################
## ##########################################################################

message('plotting fixed effect bias')

## get all the fixed effects
non.pix.eff.sd.colnames <- names(summary.metrics)[grep('_sd$', names(summary.metrics))]
## now, get the first part of the name
non.pix.eff.b.colnames <- paste0(unlist(lapply( non.pix.eff.sd.colnames, 
                                                function(x){
                                                  substr(x, 1, nchar(x)-3)
                                                })), '_bias')

(summary.metrics[,fe_int_bias := fe_int_mean - alpha])
(summary.metrics[,gauss_var_bias := gauss_var_mean - norm.var])
(summary.metrics[,matern_range_bias := matern_range_mean - sp.range])
(summary.metrics[,matern_sigma_bias := matern_sigma_mean - sqrt(sp.var)])
(summary.metrics[,clust_var_bias := clust_var_mean - clust.var])

## get some coverage probs
key.summ.metrics <- summary.metrics[,c(non.pix.eff.b.colnames, 
                                       'data.lik', 'norm.var', 'n.clust', 'fit_type', 'clust.var.cat'), with=F]
fe.mean.long <- melt(key.summ.metrics, id=c('data.lik', 'norm.var', 'n.clust', 'fit_type', 'clust.var.cat'))

## drop NAs. eg there is no gauss_prec_bias in a binom model, 
## and there is no clust_prec_bias when clust.var==NA
fe.mean.long <- fe.mean.long[!is.na(value),]

cov.width <- 0.5 ## set the coverage interval width
for(dl in c('binom', 'normal')){
  
  ## make the plot title
  pt <- sprintf('Average Param Bias for %s Observations', 
                stringr::str_to_title(dl))
  pt <- paste(pt, sprintf('\n %0.0f%% Equal Tail Intervals', cov.width*100))
  
  if(dl == 'binom') {
    fe.mean.long.sub <- subset(fe.mean.long, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
  }else if(dl == 'normal') {
    fe.mean.long.sub <- subset(fe.mean.long, data.lik == dl)
  }
  
  ## facet by ???
  fe.mean.long.sum <-  summarySE(fe.mean.long.sub, measurevar="value",
                                 groupvars=c("data.lik",'norm.var', 'n.clust', "fit_type", 'clust.var.cat', 'variable', 'norm.var'),
                                 conf.interval = 0.5)
  
  ## groups to plot lines
  ## long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')
  
  ## rename for plotting
  setnames(fe.mean.long.sum, 'n.clust', 'num.obs')
  setnames(fe.mean.long.sum, 'norm.var', 'obs.var')
  
  pd <- position_dodge(.25)
  fe_bias_summary <- ggplot(fe.mean.long.sum,
                            aes(x=log(num.obs), y=med, 
                                shape = fit_type,
                                linetype = fit_type,
                                colour = fit_type)) + 
    geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2) +
    geom_abline(intercept = 0, slope=0) +
    ## facet_wrap(. ~ variable, scales='free_y') +
    ggtitle(pt) +
    ## fix the legends a bit
    labs(color = 'Fit', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
    xlab('Num Spatial Locations Sampled') + 
    scale_x_continuous(breaks = log(sort(fe.mean.long.sum$num.obs)), ## add log x-axis labels
                       labels = paste0('ln(', sort(fe.mean.long.sum$num.obs), ')')) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = 'plain', family='sans')) + ## rotate x-axis ticks
    ylab('Bias: Est-True')
  
  if(dl=='binom'){
    fe_bias_summary <- fe_bias_summary +
      facet_wrap(. ~ variable, labeller = label_both, scales='free_y')
  }
  if(dl == 'normal'){
    fe_bias_summary <- fe_bias_summary +
      facet_grid(obs.var ~ variable, labeller = label_both, scales='free')
  }
  
  
  if(interactive()){ ## then we can view
    print(fe_bias_summary)
  }
  
  ggsave(sprintf('%s/%s_plot_03_param_bias_summary.png', compar.dir, dl),
         plot = fe_bias_summary,
         device = 'png', units = 'in',
         width = 12, height = 8)
}

## ##########################################################################
## ##########################################################################
## PLOT COVERAGE BY GP DECILE
## ##########################################################################
## ##########################################################################

message('plotting coverage by gp decile')

## cov.gp is already long, and mostly ready for plotting
## we need to average over the right things and make a column for line groupings

cov.width <- 0.8 ## set the coverage interval width
for(ii in 1:loop.params[,.N]){
  
  d.l <- dl <- loop.params[ii, dl]
  n.v <- nv <- loop.params[ii, nv]

  ## make the plot title
  if(d.l == 'normal'){
    pt <- sprintf('Average Pixel Coverage by GP Magnitude Decile given %s Observations with Var = %.03f', 
                  stringr::str_to_title(d.l), n.v)
  }else{
    pt <- sprintf('Average Pixel Coverage by GP Magnitude Decile given %s Observations', 
                  stringr::str_to_title(d.l))
  }
  pt <- paste(pt, sprintf('\n Gridded on Number of Observations by Cluster Variance || %0.0f%% Equal Tail Intervals', cov.width*100))
  
  if(d.l == 'binom') {
    cov.gp.sub <- subset(cov.gp, dl == d.l) ## NOTE norm.var==0.01 just to make sure nothing breaks
  }else if(d.l == 'normal') {
    cov.gp.sub <- subset(cov.gp, dl == d.l & nv == n.v)
  }
  
  ## facet by ???
  cov.gp.sum <- summarySE(cov.gp.sub, measurevar="obs_cov",
                          groupvars=c("dl",'nv', 'n.clust', "fit_type", 'clust.var', 'n_target_cov', 'gp.dec'),
                          conf.interval = 0.5)
  
  ## groups to plot lines
  cov.gp.sum$line_group <- paste(cov.gp.sum$fit_type, cov.gp.sum$gp.dec, sep='_')
  
  pd <- position_dodge(0.05)
  gp_cov_summary <- ggplot(cov.gp.sum,
                           aes(x=n_target_cov, y=med, 
                               shape = fit_type,
                               linetype = fit_type,
                               colour=factor(gp.dec),
                               group = line_group)) + 
    geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2) +
    scale_colour_brewer(type='div', palette = 'RdYlBu', direction=-1) + ## RdYlBu, RdYlGn, Spectral
    geom_abline(intercept = 0, slope=1) +
    facet_grid(clust.var~n.clust, scales='free_y') +
    ggtitle(pt) +
    ## fix the legends a bit
    labs(color = 'GP Decile', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
    xlab('Nominal Coverage') + 
    ylab('Observed Monte Carlo Coverage') +
    theme_dark()
  
  if(interactive()){ ## then we can view
    print(gp_cov_summary)
  }
  
  ggsave(sprintf('%s/%s_plot_04_%0.2fnormVar_gp_coverage_decile_summary.png', compar.dir, dl, nv),
         plot = gp_cov_summary,
         device = 'png', units = 'in',
         width = 12, height = 8)
}



## ##########################################################################
## ##########################################################################
## PLOT COVERAGE BY DISTANCE TO DATA OBS
## ##########################################################################
## ##########################################################################

message('plotting coverage by dist to data')

## cov.dist is already long, and mostly ready for plotting
## we need to average over the right things and make a column for line groupings

cov.width <- 0.5 ## set the coverage interval width
for(ii in 1:loop.params[,.N]){
  
  d.l <- dl <- loop.params[ii, dl]
  n.v <- nv <- loop.params[ii, nv]
  
  ## make the plot title
  if(dl == 'normal'){
    pt <- sprintf('Average Pixel Coverage by Dist to Data given %s Observations with Var = %.03f', 
                  stringr::str_to_title(dl), nv)
  }else{
    pt <- sprintf('Average Pixel Coverage by Dist to Data given %s Observations', 
                  stringr::str_to_title(dl))
  }
  pt <- paste(pt, sprintf('\n Gridded on Number of Observations by Cluster Variance || %0.0f%% Equal Tail Intervals', cov.width*100))

  
  if(dl == 'binom') {
    cov.dist.sub <- subset(cov.dist, dl == d.l) ## NOTE norm.var==0.01 just to make sure nothing breaks
  }else if(dl == 'normal') {
    cov.dist.sub <- subset(cov.dist, dl == d.l & nv == n.v)
  }
  
  ## average across relevant features
  cov.dist.sum <- summarySE(cov.dist.sub, measurevar="obs_cov",
                            groupvars=c("dl",'nv', 'n.clust', "fit_type", 'clust.var', 'n_target_cov', 'obs.dist'),
                            conf.interval = 0.5)
  
  ## drop distances without enough observations
  cov.dist.sum <- subset(cov.dist.sum, N >= 25)
  
  ## groups to plot lines
  cov.dist.sum$line_group <- paste(cov.dist.sum$fit_type, cov.dist.sum$obs.dist, sep='_')
  
  ## set factor order of clust.var
  clust.var <- as.character(cov.dist.sum$clust.var)
  clust.var[is.na(clust.var)] <- 'None'
  cov.dist.sum$clust.var <- factor(clust.var, levels = c('None', '0.01', '0.04', '0.16'))

  pd <- position_dodge(0.05)
  gp_cov_summary <- ggplot(cov.dist.sum,
                           aes(x=n_target_cov, y=med, 
                               shape = fit_type,
                               linetype = fit_type,
                               colour = obs.dist,
                               group = line_group)) + 
    geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2) +
    scale_colour_distiller(type='seq', palette = 10) + 
    geom_abline(intercept = 0, slope=1) +
    facet_grid(clust.var~n.clust, scales='free_y') +
    ggtitle(pt) +
    ## fix the legends a bit
    labs(color = 'Pixels to data', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
    xlab('Nominal Coverage') + 
    ylab('Observed Monte Carlo Coverage') + 
    theme_dark()
  
  if(interactive()){ ## then we can view
    print(gp_cov_summary)
  }
  
  ggsave(sprintf('%s/%s_plot_05_%0.2fnormVar_gp_coverage_distance_summary.png', compar.dir, dl, nv),
         plot = gp_cov_summary,
         device = 'png', units = 'in',
         width = 12, height = 8)
}





## ##########################################################################
## ##########################################################################
## ##########################################################################
## OLDER PLOT CODE BELOW...
## ##########################################################################
## ##########################################################################
## ##########################################################################
# 
# ## ##########################################################################
# ################
# ## time plots v1 - more below ## 
# ################
# ## ##########################################################################
# 
# ## get the coverages calculated
# nom.cov.names <- long.cov[,sort(unique(target_cov))] 
# nom.cov.vals <- as.numeric(gsub("[^0-9]", "",  nom.cov.names))/100
# 
# ## and assign the numeric coverage for the rows
# for(i in 1:length(nom.cov.names)){
#   long.cov[target_cov == nom.cov.names[i], n_target_cov := nom.cov.vals[i]]
# }
# 
# ## calculate noise to spatial ratio
# long.cov$noise_spatial_ratio <- long.cov$clust.var / long.cov$sp.var
# long.cov$noise_spatial_ratio[is.na(long.cov$noise_spatial_ratio)] <- 0
# 
# ## process the total fitting and predict times
# long.cov[mean.l.model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
# long.cov[mean.l.model == 'tmb', pred.time :=  as.numeric(pred_time)]
# long.cov[mean.l.model == 'inla', fit.time :=  as.numeric(fit_time)]
# long.cov[mean.l.model == 'inla', pred.time :=  as.numeric(pred_time)]
# long.cov[, total.time :=  pred.time + fit.time]
# 
# ## set cluster var NA to 0. NOTE: this is fit w/o cluster var params!
# long.cov[is.na(clust.var), clust.var:=0]
# 
# 
# loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
#                           nv=c(NA, loopvars[, unique(norm.var)]))
# 
# cov.width <- 0.5 ## set the coverage interval width
# for(ii in 1:loop.params[,.N]){
#   
#   dl <- loop.params[ii, dl]
#   nv <- loop.params[ii, nv]
#   
#   ## make the plot title
#   pt <- sprintf('Average Pixel Coverage with %s Observations | %0.0f%% E.T. Intervals', 
#                 stringr::str_to_title(dl), cov.width*100)
#   if(dl == 'normal'){
#     pt <- paste(pt, sprintf('\n Gaussian SD = %.03f', nv))
#   }
#   
#   if(dl == 'binom') {
#     long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
#   }else if(dl == 'normal') {
#     long.cov.sub <- subset(long.cov, data.lik == dl & norm.var == nv)
#   }
#   
#   ## facet by noise:spatial var ratio
#   long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
#                              groupvars=c('n_target_cov', 'fit_type', 'n.clust', 'data.lik',
#                                          'noise_spatial_ratio'))
#   
#   ## groups to plot lines
#   long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')
#   
#   ## set facet names using the labeller
#   ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
#   clust.sp.pairs <- unique(long.cov, by=c("clust.var","sp.var"))[,.(clust.var, sp.var, noise_spatial_ratio)] 
#   ## sort by unique noise_spatial_ratio order
#   clust.sp.pairs <- clust.sp.pairs[order(noise_spatial_ratio),]
#   facet_labs <- as_labeller(eval(parse(text = 
#                                          paste0('c(', paste0('`', sort(unique(long.cov$noise_spatial_ratio)), '`', '=', '\'', 
#                                                              sprintf('(Cluster Var)/(Spatial Var): %0.2f/%0.2f = %0.2f', 
#                                                                      clust.sp.pairs$clust.var,
#                                                                      clust.sp.pairs$sp.var,
#                                                                      sort(unique(long.cov$noise_spatial_ratio))), '\'', sep='', collapse=','),
#                                                 ')'))))
#   
#   pd <- position_dodge(0.05)
#   fit_coverage_CI_summary <- ggplot(long.cov.sum,
#                                     aes(x=n_target_cov, y=obs_cov, 
#                                         shape = fit_type,
#                                         linetype = fit_type,
#                                         colour=factor(n.clust), 
#                                         group = line_group),
#                                     position=position_jitter(w=0.01, h=0.01)) + 
#     geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
#     geom_line(position=pd) +
#     geom_point(position=pd) +
#     geom_abline(intercept = 0, slope = 1) +
#     facet_wrap(. ~ noise_spatial_ratio, labeller = facet_labs) +
#     ggtitle(pt) +
#     ## fix the legends a bit
#     labs(color = 'Num. Obs', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
#     xlab('Nominal Coverage') + 
#     ylab('Observed Monte Carlo Coverage')
#   if(interactive()){ ## then we can view
#     print(fit_coverage_CI_summary)
#   }
#   
#   ggsave(sprintf('%s/%s_%0.2fnormVar_pixel_coverage_summary_NoiseSpatialRatio.png', compar.dir, dl, nv),
#          plot = fit_coverage_CI_summary,
#          device = 'png', units = 'in',
#          width = 12, height = 8)
# }
# 
# ## ##########################################################################
# ## ##########
# ## TIME
# ## ##########
# ## ##########################################################################
# 
# cm.all[model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
# cm.all[model == 'tmb', pred.time :=  as.numeric(pred_time)]
# cm.all[model == 'inla', fit.time :=  as.numeric(fit_time)]
# cm.all[model == 'inla', pred.time :=  as.numeric(pred_time)]
# cm.all[, total.time :=  pred.time + fit.time]
# 
# long.cov <- data.table(gather(cm.all,
#                               operation,
#                               time_s, 
#                               fit.time:total.time,
#                               factor_key = TRUE))
# 
# 
# ## facet by observations
# long.cov.sum <-  summarySE(long.cov, measurevar="time_s",
#                            groupvars=c("operation","fit_type", 'n.clust'))
# fit_coverage_CI_summary <- ggplot(long.cov.sum,
#                                   aes(x=n.clust, y=time_s, colour=fit_type, group = fit_type)) + 
#   geom_errorbar(aes(ymin=time_s-ci, ymax=time_s+ci), width=.01) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(. ~ operation) + ggtitle(sprintf('Comparison of fit time (sec) in: %s', loopvars$data.lik[1]))
# 
# if(interactive()){
#   ## then we can view
#   print(fit_coverage_CI_summary)
# }
# 
# ggsave(sprintf('%s/%s_fit_time_nclust.png', compar.dir, loopvars$data.lik[1]),
#        plot = fit_coverage_CI_summary,
#        device = 'png', units = 'in',
#        width = 12, height = 12)
# 
# ## facet by noise to spatial signal
# long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
#                            groupvars=c("n_target_cov","fit_type", 'noise_spatial_ratio'))
# 
# fit_coverage_CI_summary <- ggplot(long.cov.sum,
#                                   aes(x=n_target_cov, y=obs_cov, colour=fit_type, group = fit_type)) + 
#   geom_errorbar(aes(ymin=obs_cov-ci, ymax=obs_cov+ci), width=.025) +
#   geom_line() +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(. ~ noise_spatial_ratio) +
#   ggtitle(sprintf('Comparison of coverage in: %s, faceted by clust.var/sp.var',
#                   loopvars$data.lik[1]))
# 
# ggsave(sprintf('%s/%s_coverage_summary_noise_to_spatial_var.png', compar.dir, loopvars$data.lik[1]),
#        plot = fit_coverage_CI_summary,
#        device = 'png', units = 'in',
#        width = 12, height = 12)
# 
# 
# ## ## average my modeling tool
# ## if(i == 1){
# ##   cm.m.tmb  <- cm[mean.l.model == 'tmb,' lapply(.SD, mean), by=mean.l.model]
# ##   cm.m.inla <- cm[mean.l.model == 'inla', lapply(.SD, mean), by=mean.l.model]
# ## }else{
# ##   cm.m.tmb  <- rbind(cm.m.tmb,
# ##                      cm[mean.l.model == 'tmb,' lapply(.SD, mean), by=mean.l.model])
# ##   cm.m.inla <- rbind(cm.m.inla,
# ##                      cm[mean.l.model == 'inla', lapply(.SD, mean), by=mean.l.model])
# ## }
