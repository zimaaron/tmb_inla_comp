## #############
## #############
## VALIDATION ##
## #############
## #############

message('---- ON SCRIPT 6: running validation')

## update the tracker
write.table(x=matrix(c(sim.loop.ct, 6), ncol=2), append=T,
            file = paste0(jobtrack.dir, 
                          sprintf('exp_%04d_iter_%04d.csv', exp.lvid, exp.iter)), sep=',',
            row.names=F, col.names = F)

## 1) summarize fitted params
## 2) big plots showing difference in fits
## 3) calcualte and summarize predictive metrics  

## ###################################
## 1) summarize fitted param values ##
## ###################################
message('------ making table of true and estimated params')

## make a dt for comparing results and store some relevant computing and mesh params
res <- data.table(st_mesh_nodes = rep(nrow(epsilon_tmb_draws),2))
res[, cores           := rep(cores,2)]
res[, s_mesh_max_edge := rep(eval(parse(text = mesh_s_params))[2],2)]
res[, s_mesh_cutoff   := rep(eval(parse(text = mesh_s_params))[1],2)]
res[, draws           := c(ndraws,ndraws)]

## time variables
res[,fit_time  := c(fit_time_inla,fit_time_tmb)]
res[,pred_time := c(totalpredict_time_inla,totalpredict_time_tmb)]
res[,pt_tmb_sdreport_time := c(NA,tmb_sdreport_time)]
res[,pt_get_draws_time := c(inla_get_draws_time,tmb_get_draws_time)]

## convergence
res[, convergence := c(tmb.converge, inla.converge)]
## converge attempts
res[, convergence.fails := c(tmb.converge.fails, inla.converge.fails)]

## fe coefficients
if(!is.null(alpha) | !is.null(betas)){
  res[, paste0('fe_',res_fit$names.fixed,'_mean') := as.data.frame(rbind(res_fit$summary.fixed$mean, SD0$par.fixed[1:length(res_fit$names.fixed)]))]
  res[, paste0('fe_',res_fit$names.fixed,'_sd')   := as.data.frame(rbind(res_fit$summary.fixed$sd, sqrt(diag(SD0$cov.fixed))[1:length(res_fit$names.fixed)]))]
}

## cluster prec
if(!is.null(clust.var)){
  res[,clust_prec := unname(c(res_fit$summary.hyperpar[grep('clust.id',rownames(res_fit$summary.hyperpar)),4],
                                   SD0$value[grep('clust_prec', names(SD0$value))]))]
  res[,clust_prec_sd := c(res_fit$summary.hyperpar[grep('clust.id',rownames(res_fit$summary.hyperpar)),2],
                             sqrt(SD0$cov[grep('clust_prec', names(SD0$value)), grep('clust_prec', names(SD0$value))])) ]
}

## normal data obs prec
if(data.lik == 'normal'){
  res[,gauss_prec := unname(c(res_fit$summary.hyperpar[grep('Gaussian',rownames(res_fit$summary.hyperpar)),4],
                                   SD0$value[grep('gauss_prec', names(SD0$value))]))]
  res[,gauss_prec_sd := c(res_fit$summary.hyperpar[grep('Gaussian',rownames(res_fit$summary.hyperpar)),2],
                             sqrt(SD0$cov[grep('gauss_prec', names(SD0$value)), grep('gauss_prec', names(SD0$value))])) ]
}

## hyperparameters
res[,matern_logtau_mean := unname(c(res_fit$summary.hyperpar[grep('Theta1',rownames(res_fit$summary.hyperpar)),4],
                                      SD0$par.fixed['log_tau']))]
res[,matern_logtau_sd := c(res_fit$summary.hyperpar[grep('Theta1',rownames(res_fit$summary.hyperpar)),2],
                             sqrt(SD0$cov.fixed['log_tau','log_tau'])) ]

res[,matern_logkappa_mean := c(res_fit$summary.hyperpar[grep('Theta2',rownames(res_fit$summary.hyperpar)),4],
                                 SD0$par.fixed['log_kappa']) ]
res[,matern_logkappa_sd := c(res_fit$summary.hyperpar[grep('Theta2',rownames(res_fit$summary.hyperpar)),2],
                               sqrt(SD0$cov.fixed['log_kappa','log_kappa'])) ]

## add extra row to filled with the truth
res <- rbind(lapply(1:ncol(res), function(x){NA}), res)

## slot in the truth and also make a list of all params in the model
params <- NULL
if(!is.null(alpha)){ res[1, fe_int_mean := alpha]; params <- c(params, 'alpha')}
if(!is.null(betas) & is.null(alpha)) { res[1, grep('fe.*med', colnames(res)) := betas]; params <- c(params, rep('beta', length(betas)))}
if(!is.null(betas) & !is.null(alpha)){ res[1, grep('fe.*med', colnames(res))[-1] := betas]; params <- c(params, rep('beta', length(betas)))}
if(!is.null(clust.var)) {res[1, clust_prec := 1 / clust.var];params <- c(params, 'clust.prec')}
if(data.lik == 'normal') {res[1, gauss_prec := 1 / norm.var]; params <- c(params, 'gauss.prec')}
res[1, matern_logtau_mean := log(sp.tau)]; params <- c(params, 'logtau')
res[1, matern_logkappa_mean := log(sp.kappa)]; params <- c(params, 'logkappa')

## if(nperiods > 1){
##   res.true.params <- c(res.true.params, c(t.rho, NA))
## }

rr <- data.table(item=colnames(res))
rr <- cbind(rr, t(res))
names(rr) <- c('quantity','TRUE', 'R-INLA','TMB')
rr$diff <- rr[,3]-rr[,4]

write.table(x = rr, row.names = FALSE, sep=',', 
          file = sprintf('%s/validation/experiment%04d_iter%04d_param_summary_table.csv', out.dir, exp.lvid, exp.iter))

## we can now plot this table with: grid.table(rr)

## ###########################
## 2) make a bunch of plots ##
## ###########################
## pdf(sprintf('%s/validation/inla_tmb_summary_comparison_plot_%i_new.pdf',out.dir, exp.iter), height=15,width=30)
## TODO? one overall plot? or somehow stitch together later...

## ~~~~~~~~~~~~~~~~~~~~~~~~~
## plot the table of results
## ~~~~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting the summary table')

png(sprintf('%s/validation/experiment%04d_iter%04d_plot_01_summary_table.png', out.dir, exp.lvid, exp.iter),
    height=7, width=9, units = 'in', res = 250)
cols <- names(rr)[2:5]
rr[,(cols) := round(.SD, 3), .SDcols=cols]
grid.table(rr)
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~
## plot priors and posteriors
## ~~~~~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting priors and posteriors')

## assume:
## 1) alphas are normal
## 2) logtau and logkappa are normal
## 3) log clust precision is gamma

## NOTE: also assume that intercept is the first 'beta' and that betas are listed first!

png(sprintf('%s/validation/experiment%04d_iter%04d_plot_02_parameter_densities.png', out.dir, exp.lvid, exp.iter),
    height=9, width=9, units = 'in', res = 250)

num.dists <- length(params)

## set layout: start with a square, and drop rows is not needed
par.dim <- rep(ceiling(sqrt(num.dists)), 2)
while((par.dim[1]-1)*par.dim[2] >= num.dists){
  par.dim[1] <- par.dim[1]-1
} 
par(mfrow = par.dim)

for(ii in 1:num.dists){

  param <- params[ii]
  message(sprintf('-------- plotting prior and post for: %s', param))
  
  ## get prior curves and posterior draws
  if(param == 'alpha'){
    param.name <- names(res_fit$marginals.fixed)[ii]

    if(param.name == 'int'){
      true.val <- alpha
    }else{
      true.val <- betas[ii - 1]
    }
    
    prior.mean <- alphaj.pri[1]
    prior.sd   <- alphaj.pri[2]

    tmb.post.draws  <- alpha_tmb_draws[ii, ]
    inla.post.draws <- pred_l[ii, ]

    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
  }
  
  if(param == 'beta'){

    param.name <- names(res_fit$marginals.fixed)[ii]
    if(param.name == 'int'){
      true.val <- alpha
    }else{
      true.val <- betas[ii - 1]
    }
    
    prior.mean <- alphaj.pri[1]
    prior.sd   <- sqrt(alphaj.pri[2])
    
    tmb.post.draws  <- betas_tmb_draws[ii, ]
    inla.post.draws <- pred_l[ii, ]
    
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
  }
  
  if(param == 'logkappa'){
    
    true.val <- logkappa

    mesh.pri   <- param2.matern.orig(mesh_s) ## theta1 = logtau, theta2 = logkappa 
    prior.mean <- mesh.pri$theta.prior.mean[2]
    prior.prec <- mesh.pri$theta.prior.prec[2, 2]
    prior.sd   <- 1 / sqrt(prior.prec)

    tmb.post.draws  <- log_kappa_tmb_draws
    inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[['Theta2 for space']][, 1],
                                    size = ndraws,
                                    replace = TRUE, 
                                    res_fit$marginals.hyperpar[['Theta2 for space']][, 2])

    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

  
    tmb.post.median  <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log kappa"
  }
  if(param == 'logtau'){

    true.val <- logtau

    mesh.pri   <- param2.matern.orig(mesh_s) ## theta1 = logtau, theta2 = logkappa 
    prior.mean <- mesh.pri$theta.prior.mean[1]
    prior.prec <- mesh.pri$theta.prior.prec[1, 1]
    prior.sd   <- 1 / sqrt(prior.prec)
    
    tmb.post.draws <- log_tau_tmb_draws
    inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[['Theta1 for space']][, 1],
                                    size = ndraws,
                                    replace = TRUE, 
                                    res_fit$marginals.hyperpar[['Theta1 for space']][, 2])

    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log tau"
  }
  if(param == 'gauss.prec'){

    true.val <- 1 / norm.var

    prior.u  <- norm.prec.pri[1]
    prior.a <- norm.prec.pri[2]
    
    tmb.post.draws <- 1 / exp(log_gauss_sigma_draws * 2)
    inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[['Precision for the Gaussian observations']][, 1],
                                    size = ndraws,
                                    replace = TRUE, 
                                    res_fit$marginals.hyperpar[['Precision for the Gaussian observations']][, 2])

    if(true.val==Inf) {
      ## get a safe range for plotting, and get the prior
      xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws),
                              probs=c(.0001, .9999)) ## avoid crazy extremes
    }else{
      ## get a safe range for plotting, and get the prior
      xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                              probs=c(.0001, .9999)) ## avoid crazy extremes
    }
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    
   
    x.prior <- seq(xlim[1], xlim[2], len = 1000)
    y.prior <- dPCPriPrec(x.prior, u = prior.u, a=prior.a, give_log = 0)
    
    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log gauss prec"
  }
  if(param == 'clust.prec'){

    true.val <- 1 / clust.var

    prior.u <- clust.prec.pri[1]
    prior.a <- clust.prec.pri[2]

    tmb.post.draws <- 1 / exp(log_clust_sigma_draws * 2)
    inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[[names(res_fit$marginals.hyperpar)[grep('clust.id', names(res_fit$marginals.hyperpar))]]][, 1],
                                    size = ndraws,
                                    replace = TRUE, 
                                    res_fit$marginals.hyperpar[[names(res_fit$marginals.hyperpar)[grep('clust.id', names(res_fit$marginals.hyperpar))]]][, 2])
    
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior <- seq(xlim[1], xlim[2], len = 1000)
    y.prior <- dPCPriPrec(x.prior, u = prior.u, a=prior.a, give_log = 0)
    
    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log clust prec"
  }

  ## get posterior samples (we'll use density curves)
  tmb.dens <- density(tmb.post.draws)
  inla.dens <- density(inla.post.draws)
  xrange <- xlim ## range(c(x.prior, tmb.dens$x, inla.dens$x))
  yrange <- range(c(y.prior, tmb.dens$y, inla.dens$y))

  prior.col <- "black"
  tmb.col   <- "red"
  inla.col  <- "blue"

  ## setup plot and plot prior
  plot(x.prior, y.prior, pch = ".", col = prior.col, main = param.name,
       xlim = xrange, ylim = yrange)
  lines(x.prior, y.prior, col = prior.col)

  ## plot tmb post and data
  lines(tmb.dens$x,tmb.dens$y, col = tmb.col)
  points(tmb.post.draws, rep(0, ndraws), col = alpha(tmb.col, 0.1), cex = 2, pch = '|')
  abline(v = tmb.post.median, col = tmb.col, lwd = 2)
    
  ## plot inla post and data
  lines(inla.dens$x,inla.dens$y, col = inla.col)
  points(inla.post.draws, rep(0, ndraws), col = alpha(inla.col, 0.1), cex = 2, pch = '|')
  abline(v = inla.post.median, col = inla.col, lwd = 2)

  ## plot truth
  abline(v = true.val, col = prior.col, lwd = 2)
  
  ## add legend
  legend("topright", legend = c("prior/truth", "tmb", "inla"), col = c(prior.col, tmb.col, inla.col), lwd = rep(2, 3))
}
dev.off()



## plot results in logit space or in prevalence space?
## TODO (?) plot.in.logit.space <- FALSE 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## plot comparisons of INLA and TMB fields
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting the comparison of INLA and TMB fields')

## randomly select pixels for plotting tmb v inla scatter
samp <- sample(cellIdx(ras_med_inla[[1]]),1e4) 

for(sum.meas in c('median','stdev')){ 
  
  if(sum.meas=='median'){
    
    if(data.lik == 'binom'){
      true  <- invlogit(true.rast[[1]])
      rinla <- ras_med_inla_p[[1]]
      rtmb  <- ras_med_tmb_p[[1]]
    }else if(data.lik == 'normal'){
      true <- true.rast[[1]]
      rinla <- ras_med_inla[[1]]
      rtmb  <- ras_med_tmb[[1]]
    }
    all.vec <- c(as.vector(rtmb), as.vector(rinla), as.vector(true))
    all.diff.vec <- as.vector(rinla - rtmb)
    diff.truth.zrange <- range(c( as.vector(rtmb-true), as.vector(rinla-true) ), na.rm=T)
    rast.list <- list('TRUE' = true,
                      'TMB' = rtmb,
                      'INLA' = rinla)
    
    png(sprintf('%s/validation/experiment%04d_iter%04d_plot_03_median_rasters.png', out.dir, exp.lvid, exp.iter),
        height=12, width=12, units = 'in', res = 250)
  } ## sum.meas==median
  
  if(sum.meas=='stdev'){
    if(data.lik == 'binom'){
      rinla <- ras_sdv_inla_p[[1]]
      rtmb  <- ras_sdv_tmb_p[[1]]
    }else if(data.lik == 'normal'){
      rinla <- ras_sdv_inla[[1]]
      rtmb  <- ras_sdv_tmb[[1]]
    }
    all.vec <- c(as.vector(rtmb), as.vector(rinla))
    all.diff.vec <- as.vector(rinla - rtmb)
    rast.list <- list('TMB' = rtmb,
                      'INLA' = rinla)
    
    png(sprintf('%s/validation/experiment%04d_iter%04d_plot_04_stdev_rasters.png', out.dir, exp.lvid, exp.iter),
        height=8, width=8, units = 'in', res = 250)
  } ## sum.meas==stdev
  
  layout(matrix(1:length(rast.list) ^ 2, byrow = T, ncol = length(rast.list)))
  
  tmp <- subset(dt, period_id==1) ## for s-t

  ## get limits
  rast.zrange <- range(all.vec, na.rm = T)
  diff.zrange <- range(all.diff.vec, na.rm = T)
  
  ## set some plot args for mean and stdev separately
  if(sum.meas=='median'){
    ## set the legend.width to get nice plots
    lw <- 5
    tick.space.val  <- 0.5
    tick.space.truth.diff <- 0.25
    tick.space.diff <- 0.1
    mar <- c(0, 0, 1.4, 8)
  }
  if(sum.meas=='stdev'){
    ## set the legend.width to get nice plots
    lw <- 2
    tick.space.val  <- 0.1
    tick.space.diff <- 0.05
    mar <- c(0, 0, 1.4, 4)
  }
  
 

  for(i in 1:length(rast.list)){
    for(j in 1:length(rast.list)){

      if(i == j){
        ## plot raster
        par(mar = mar, bty='n')
        plot(rast.list[[i]],  maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend.width = lw,
             axis.args=list(at=c(seq(rast.zrange[1], rast.zrange[2], by=tick.space.val), rast.zrange[2]),
                            labels=round(c(seq(rast.zrange[1], rast.zrange[2], by=tick.space.val), rast.zrange[2]), 3), 
                            cex.axis=0.6),
             legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main=paste0(names(rast.list)[i], ': ', sum.meas),
             zlim=rast.zrange)
      }

      if(i < j){
        ## plot scatter
        par(mar = c(4, 4, 2, 2),bty='n')
        plot(x=as.vector(rast.list[[j]])[samp],
             y=as.vector(rast.list[[i]])[samp],
             xlab='',
             ylab='',
             cex=.05, pch=19,
             main=paste0('(sub)SCATTER (', names(rast.list)[i], ' vs ', names(rast.list)[j], '): '))
        lines(x=rast.zrange, y=rast.zrange, col='red')
        title(xlab=names(rast.list)[j], line=-1, cex.lab=1.2)
        title(ylab=names(rast.list)[i], line=-1, cex.lab=1.2)
      }

      if(i > j){
        ## plot difference rasters
        
        if(sum.meas=='median' & j == 1){
          ## these are for truth - *_estiamte
          diff.col <- rev(plasma(100))
          dz  <- diff.truth.zrange
          tsd <- tick.space.truth.diff
          
        }else{
          ## these are for inla_est - tmb_est
          diff.col <- rev(cividis(100))
          dz  <- diff.zrange
          tsd <- tick.space.diff
        }
        par(mar = mar, bty='n')
        plot(rast.list[[i]] - rast.list[[j]],  maxpixel=1e7, col=diff.col, axes=FALSE, legend.width = lw,
             axis.args=list(at=sort(c(0, seq(dz[1], dz[2], by=tsd), dz[2])),
                            labels=round(sort(c(0, seq(dz[1], dz[2], by=tsd), dz[2])), 3),
                            cex.axis=0.6),
             legend.args=list(text='', side=2, font=1, line=0, cex=0.1),
             main=paste0('Diff: ', names(rast.list)[i], ' - ', names(rast.list)[j]),
             zlim=dz)
        points( x=tmp$long,y=tmp$lat, pch=19, cex=.1 )

      }
      
    } ## j
  } ## i

  dev.off()
  
} ## sum.meas

## ~~~~~~~~~~~~~~~~~~~~~~
## plot caterpillar plots
## ~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting caterpillar plots')

png(sprintf('%s/validation/experiment%04d_iter%04d_plot_05_spatial_re_caterpillars.png', out.dir, exp.lvid, exp.iter),
      height=8, width=12, units = 'in', res = 250)

layout(matrix(1, 1, 1, byrow = TRUE))

## Compare mean and distribution of random effects
summ_gp_tmb  <- t(cbind((apply(epsilon_tmb_draws, 1, quantile,probs=c(.1,.5,.9)))))
summ_gp_inla <- t(cbind((apply(pred_s, 1, quantile,probs=c(.1,.5,.9)))))

## all time-space random effects
plot_d <- data.table(tmb_median = summ_gp_tmb[, 2],inla_median = summ_gp_inla[, 2],
                     tmb_low    = summ_gp_tmb[, 1],inla_low    = summ_gp_inla[, 1],
                     tmb_up     = summ_gp_tmb[, 3],inla_up     = summ_gp_inla[, 3])

plot_d$period <- factor(rep(1:nperiods, each=nrow(plot_d)/nperiods))
plot_d$loc    <- rep(1:(nrow(plot_d)/nperiods), rep=nperiods)
plot_d[, absdiff := abs(tmb_median-inla_median)]

## get the node locations and add them on
nodelocs <- mesh_s$loc
plot_d <- cbind(plot_d, nodelocs)

if(nrow(plot_d) > 2500)
  plot_d <- plot_d[sample(nrow(plot_d), 2500, replace=F), ]


## ## plot spatial RE differences
## ggplot(plot_d, aes(x=tmb_median,y=inla_median,col=period)) + theme_bw() +
## geom_point() + geom_line(aes(group=loc)) + geom_abline(intercept=0,slope=1,col='red') +
## ggtitle('Posterior Medians of Random Effects at Mesh Nodes, TMB v R-INLA. Connected dots same location different periods. ')

## plot locations where they are relatively different, are they near or far from data?
biggdiff <- plot_d[which(plot_d$absdiff > quantile(plot_d$absdiff, prob=0.80)), ]

plot(simple_polygon, main='Mesh nodes sized by abs difference TMB and R-INLA')
points(x=dt$long, y=dt$lat, pch=19, cex=(dt$N / max(dt$N))*.1)
points(x=plot_d$V1, y=plot_d$V2, pch=1, cex=plot_d$absdiff*5, col='red')

## catterpillar plot
plot_d <- plot_d[order(period,tmb_median)]
plot_d[ ,i := seq(1,.N), by = period]
gg_cat <- ggplot(plot_d, aes(i, tmb_median, col=i)) + theme_bw() + # [seq(1, nrow(plot_d), 5)]
          geom_linerange(aes(ymin = tmb_low, ymax = tmb_up), col='red', size=.8, alpha=.3) +
          geom_linerange(aes(x=i,ymin = inla_low, ymax = inla_up), col='blue', size=.8, alpha=.3) +
          facet_wrap(~period) +
          ggtitle('Comparison of random effects (10% to 90% quantiles) | BLUE == R-INLA | RED == TMB')
print(gg_cat)

dev.off()

## ###############################################
## 3) generate and summarize predictive metrics ##  
## ###############################################
message('------ making metrics to compare true and estimated surfaces')

## NOTE: these are in latent space for all models!
all.preds <- data.table(rbind(pred_tmb, pred_inla))

## make a data.table with prediction draws, model type, and truth
## NOTE! the truth and the median.fit are in logit-space!
non.na.idx <- which(!is.na(values(simple_raster)))
d <- data.table(truth        = rep(values(true.rast)[non.na.idx], 2),
                model        = c(rep('tmb', nrow(pred_tmb)),
                                 rep('inla', nrow(pred_inla))), 
                median.fit   = apply(all.preds, 1, median))

## get some coverage probs
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

## get error and pred var
d[, error := truth - median.fit]
d[, var := apply(all.preds, 1, var)]

if(data.lik == 'binom'){
  ## NOTE: these are in latent space for all models!
  all.preds <- data.table(rbind(pred_tmb_p, pred_inla_p))

  ## make a data.table with prediction draws, model type, and truth
  ## NOTE! the truth and the median.fit are in logit-space!
  non.na.idx <- which(!is.na(values(true.rast.p)))
  d[,`:=`(truth.p  = rep(values(true.rast.p)[non.na.idx], 2),
          model.p  = c(rep('tmb', nrow(pred_tmb_p)),
                       rep('inla', nrow(pred_inla_p))), 
          median.fit.p   = apply(all.preds, 1, median))]
  
  ## get some coverage probs
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
    d[, paste0('p.pix.cov.',c)] = d[['truth.p']] >= li & d[['truth.p']] <= ui
  }

  ## get error and pred var
  d[, error.p := truth.p - median.fit.p]
  d[, var.p := apply(all.preds, 1, var)]
}

message('-------- making final table of summary metrics and scores')
## summarize across pixels
surface.metrics <- data.table(cbind(mean.l      = d[, .(truth = mean(truth, na.rm = T)), by = c('model')], 
                                    mean.l.est  = d[, .(est = mean(median.fit, na.rm = T)), by = c('model')]$est, 
                                    bias        = d[, .(bias = mean(error, na.rm = T)), by = c('model')]$bias,
                                    rmse        = d[, .(rmse = sqrt(mean(error ^ 2, na.rm = T))), by = c('model')]$rmse,
                                    cor         = d[, .(cor = my.cor(truth, median.fit)), by = c('model')]$cor,
                                    CoV         = d[, .(cov = mean(error / var, na.rm = T)), by = c('model')]$cov,
                                    crps        = d[, .(crps = mean(crpsNormal(truth, median.fit, var), na.rm = T)), by = c('model')]$crps,
                                    cov25       = d[, .(cov25 = mean(pix.cov.25, na.rm = T)), by = c('model')]$cov25,
                                    cov50       = d[, .(cov50 = mean(pix.cov.50, na.rm = T)), by = c('model')]$cov50,
                                    cov80       = d[, .(cov80 = mean(pix.cov.80, na.rm = T)), by = c('model')]$cov80,
                                    cov90       = d[, .(cov90 = mean(pix.cov.90, na.rm = T)), by = c('model')]$cov90,
                                    cov95       = d[, .(cov95 = mean(pix.cov.95, na.rm = T)), by = c('model')]$cov95))

if(data.lik == 'binom'){
  surface.metrics.p <- data.table(cbind(
                                    mean.p      = d[, .(truth.p = mean(truth.p, na.rm = T)), by = c('model')]$truth.p, 
                                    mean.p.est  = d[, .(est.p = mean(median.fit.p, na.rm = T)), by = c('model')]$est.p, 
                                    bias.p      = d[, .(bias.p = mean(error.p, na.rm = T)), by = c('model')]$bias.p,
                                    rmse.p      = d[, .(rmse.p = sqrt(mean(error.p ^ 2, na.rm = T))), by = c('model')]$rmse.p,
                                    cor.p       = d[, .(cor.p = my.cor(truth.p, median.fit.p)), by = c('model')]$cor.p,
                                    CoV.p       = d[, .(cov.p = mean(error.p / var.p, na.rm = T)), by = c('model')]$cov.p
                                    ))
  surface.metrics <- cbind(surface.metrics, surface.metrics.p)
}

## add on things from the res tab
res.addon <- t(rr[, c('R-INLA', 'TMB')])
colnames(res.addon) <- rr[, quantity]

## and addon loovars for easy summary later
lv.addon <- rbind(loopvars[exp.lvid, ], loopvars[exp.lvid, ])
summary.metrics <- cbind(surface.metrics, res.addon, lv.addon)

## note iteration
summary.metrics[, iter := exp.iter]

## save
write.csv(summary.metrics, sprintf('%s/validation/experiment%04d_iter%04d_summary_metrics.csv', out.dir, exp.lvid, exp.iter))