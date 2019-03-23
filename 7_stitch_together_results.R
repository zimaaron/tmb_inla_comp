## this script pulls together some overall results from an experiment run
main.dir.names <- c("2019_03_21_05_25_32", ## binomial
                    "2019_03_21_06_11_57") ## normal

##utility function from the interwebs
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
    # New version of length which can handle NA's: if na.rm==T, don't count them
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
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- plyr::rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

for(main.dir.name in main.dir.names){

  main.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', main.dir.name)

  compar.dir <- sprintf('%s/comparisons/', main.dir)
  dir.create(compar.dir, recursive = T, showWarnings = F)

  ## read in all experiments
  loopvars <- fread(file = paste0(main.dir, '/loopvars.csv'), stringsAsFactors = F)

  ## duplicate and add tmb and inla indicator
  ## loopres <- rbind(loopvars, loopvars)
  ## loopres[, model := c(rep('tmb', nrow(loopvars)), rep('inla', nrow(loopvars)))]

  ##TODO save the loopvars with this info!
  colnames(loopvars) <- c('reg', ## 1
                          'year_list',
                          'cov_names',
                          'cov_measures',
                          'betas', ## 5
                          'alpha',
                          'sp.range',
                          'sp.var',
                          'sp.alpha',
                          'nug.var', ## 10
                          't.rho',
                          'mesh_s_params',
                          'n.clust',
                          'm.clust',
                          'sample.strat', ## 15
                          'cores',
                          'ndraws',
                          'alphaj.pri',
                          'nug.prec.pri',
                          'inla.int.strat', ## 20
                          'inla.approx', 
                          'Nsim',
                          'data.lik',
                          'norm.var',
                          'norm.prec.pri', ## 25
                          'bias.correct',
                          'sd.correct')

  ## first, for each completed experiment, we can read the experimends in
  for(i in 1:nrow(loopvars)){

    ## complete metrics across all monte carlo iterations
    cm.fn <- sprintf('%s/%i/validation/surface_metrics_complete.csv', main.dir, i)
    if(file.exists(cm.fn)){
      cm <- fread(cm.fn, stringsAsFactors = F)
      

      ## we also need to load in the param.summary.table from each iteration
      for(j in 1:loopvars[i, Nsim]){
        pst.fn <- sprintf('%s/%i/validation/iter%04d_param_summary_table.csv', main.dir, i, j)
        pst <- fread(pst.fn, stringsAsFactors = F)
        rn <- colnames(pst)
        pst <- data.table(t(pst))
        rownames(pst) <- rn
        colnames(pst) <- as.character(pst[1, ])
        pst <- pst[4:3, ]
        pst[, model := c('tmb', 'inla')]

        ## stick on the loopvars
        pst <- cbind(pst, rbind(loopvars[i, ], loopvars[i, ]))
        if(loopvars$data.lik[1] == 'binom')

          if(j == 1){
            pst.i <- pst
          }else{
            pst.i <- rbind(pst.i, pst)
          }
      }

      ##TODO remove this step! bias should be quite small, if it is
      ##large it is most likely due to convergence issues which i need
      ##to keep track of in later runs...
      cm <- cbind(cm, pst.i)
      cm <- cm[bias < 1, ]

      cm[, experiment := i]

      if(i == 1){
        cm.all <- cm
        
      }else{
        cm.all <- rbind.fill(cm.all, cm)
      }
      
    }else{
      ## skip that iteration for now
    }
  }

  cm.all <- data.table(cm.all)
  

  ## ####################################################################
  ## ####################################################################
  ## now that should be everything!
  ## we can make a bunch of plots
  ## ####################################################################
  ## ####################################################################

  ##make new labels for INLA_EB, INLA_CCD, TMB
  cm.all[inla.int.strat == 'eb' & mean.l.model == 'inla', fit_type := 'INLA_EB']
  cm.all[inla.int.strat == 'ccd' & mean.l.model == 'inla', fit_type := 'INLA_CCD']
  cm.all[mean.l.model == 'tmb', fit_type := 'TMB']

  library(tidyr)
  library(ggplot2)
  long.cov <- data.table(gather(cm.all,
                                target_cov,
                                obs_cov, 
                                cov25:cov95,
                                factor_key = TRUE
                                ))
  long.cov[target_cov == 'cov25', n_target_cov := 0.25]
  long.cov[target_cov == 'cov50', n_target_cov := 0.50]
  long.cov[target_cov == 'cov80', n_target_cov := 0.80]
  long.cov[target_cov == 'cov90', n_target_cov := 0.90]
  long.cov[target_cov == 'cov95', n_target_cov := 0.95]

  long.cov$noise_spatial_ratio <- long.cov$nug.var / long.cov$sp.var
  long.cov$noise_spatial_ratio[is.na(long.cov$noise_spatial_ratio)] <- 0

  long.cov[model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
  long.cov[model == 'tmb', pred.time :=  as.numeric(pred_time)]
  long.cov[model == 'inla', fit.time :=  as.numeric(fit_time)]
  long.cov[model == 'inla', pred.time :=  as.numeric(pred_time)]
  long.cov[, total.time :=  pred.time + fit.time]

  ## ##########
  ## COVERAGE
  ## ##########
  
  ## facet by observations
  long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
                             groupvars=c("n_target_cov","fit_type", 'n.clust'))
  fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                    aes(x=n_target_cov, y=obs_cov, colour=fit_type, group = fit_type)) + 
  geom_errorbar(aes(ymin=obs_cov-ci, ymax=obs_cov+ci), width=.025) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(. ~ n.clust) + ggtitle(sprintf('Comparison of coverage in: %s', loopvars$data.lik[1]))

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
  ggtitle(sprintf('Comparison of coverage in: %s, faceted by nug.var/sp.var',
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
                                factor_key = TRUE
                                ))

  
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
  ggtitle(sprintf('Comparison of coverage in: %s, faceted by nug.var/sp.var',
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

  











}
