## #############
## #############
## VALIDATION ##
## #############
## #############

## 1) summarize fitted params
## 2) big plots showing difference in fits
## 3) calcualte and summarize predictive metrics  

## ###################################
## 1) summarize fitted param values ##
## ###################################

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
res[, convergence := c(tmb.pd.converge, inla.converge)]

## fe coefficients
if(!is.null(alpha) | !is.null(betas)){
  res[, paste0('fe_',res_fit$names.fixed,'_med') := rbind(res_fit$summary.fixed$mean, SD0$par.fixed[1:length(res_fit$names.fixed)])]
  res[, paste0('fe_',res_fit$names.fixed,'_sd')   := rbind(res_fit$summary.fixed$sd, sqrt(diag(SD0$cov.fixed))[1:length(res_fit$names.fixed)])]
}

## nugget
if(!is.null(nug.var)){
  res[,nug_prec := unname(c(res_fit$summary.hyperpar[grep('nug.id',rownames(res_fit$summary.hyperpar)),4],
                                   SD0$value[grep('nugget_prec', names(SD0$value))]))]
  res[,nug_prec_sd := c(res_fit$summary.hyperpar[grep('nug.id',rownames(res_fit$summary.hyperpar)),2],
                             sqrt(SD0$cov[grep('nugget_prec', names(SD0$value)), grep('nugget_prec', names(SD0$value))])) ]
}

## normal var
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

## add extra row to fille with the truth
res <- rbind(lapply(1:ncol(res), function(x){NA}), res)

## slot in the truth and also make a list of all params in the model
params <- NULL
if(!is.null(alpha)){ res[1, fe_int_mean := alpha]; params <- c(params, 'alpha')}
if(!is.null(betas) & is.null(alpha)) { res[1, grep('fe.*med', colnames(res)) := betas]; params <- c(params, rep('beta', length(betas)))}
if(!is.null(betas) & !is.null(alpha)){ res[1, grep('fe.*med', colnames(res))[-1] := betas]; params <- c(params, rep('beta', length(betas)))}
if(!is.null(nug.var)) {res[1, nugget_prec := 1 / nug.var];params <- c(params, 'nug.prec')}
if(data.lik == 'normal') {res[1, gauss_prec := 1 / norm.var]; params <- c(params, 'gauss.prec')}
res[1, matern_logtau_mean := log(sp.tau)]; params <- c(params, 'logkappa')
res[1, matern_logkappa_mean := log(sp.kappa)]; params <- c(params, 'logtau')

## if(nperiods > 1){
##   res.true.params <- c(res.true.params, c(t.rho, NA))
## }

rr <- data.table(item=colnames(res))
rr <- cbind(rr, t(res))
names(rr) <- c('quantity','TRUE', 'R-INLA','TMB')
rr$diff <- rr[,3]-rr[,4]

write.csv(x = rr, row.names = FALSE, 
          file = sprintf('%s/validation/iter%04d_param_summary_table.csv', out.dir, iii))

## we can now plot this table with: grid.table(rr)

## ####################
## 2) setup big plot ##
## ####################
## pdf(sprintf('%s/validation/inla_tmb_summary_comparison_plot_%i_new.pdf',out.dir, iii), height=15,width=30)
## TODO? one overall plot? or somehow stitch together later...


## ~~~
## plot the table of results
## ~~~
png(sprintf('%s/validation/iter%04d_plot_01_sumary_table.png', out.dir, iii),
    height=7, width=9, units = 'in', res = 250)
cols <- names(rr)[2:5]
rr[,(cols) := round(.SD, 3), .SDcols=cols]
grid.table(rr)
dev.off()

## ~~~
## plot priors and posteriors
## ~~~

## assume:
## 1) alphas are normal
## 2) logtau and logkappa are normal
## 3) log nugget precision is gamma

## NOTE: also assume that intercept is the first 'beta' and that betas are listed first!

png(sprintf('%s/validation/iter%04d_plot_02_parameter_densities.png', out.dir, iii),
    height=9, width=9, units = 'in', res = 250)

num.dists <- length(params)

par(mfrow = rep(ceiling(sqrt(num.dists)), 2))

for(ii in 1:num.dists){

  param <- params[ii]
  
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

    xlim       <- range(c(tmb.post.draws, inla.post.draws, true.val)) 
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
    
    tmb.post.draws  <- alpha_tmb_draws[ii, ]
    inla.post.draws <- pred_l[ii, ]
    
    xlim       <- range(c(tmb.post.draws, inla.post.draws, true.val)) 
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

    xlim       <- range(c(tmb.post.draws, inla.post.draws, true.val)) 

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

    xlim       <- range(c(tmb.post.draws, inla.post.draws, true.val)) 

    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log tau"
  }
  if(param == 'gauss.prec'){

    true.val <- 1 / norm.var

    prior.shape  <- norm.prec.pri[1]
    prior.iscale <- norm.prec.pri[2]
    
    tmb.post.draws <- 1 / exp(log_gauss_sigma_draws * 2)
    inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[['Precision for the Gaussian observations']][, 1],
                                    size = ndraws,
                                    replace = TRUE, 
                                    res_fit$marginals.hyperpar[['Precision for the Gaussian observations']][, 2])

    xlim       <- range(c(tmb.post.draws, inla.post.draws, true.val)) 
    
    x.prior <- seq(xlim[1], xlim[2], len = 1000)
    y.prior <- dgamma(x.prior, shape = prior.shape, scale = 1 / prior.iscale)
    
    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log gauss prec"
  }
  if(param == 'nug.prec'){

    true.val <- 1 / nug.var

    prior.shape  <- nug.prec.pri[1]
    prior.iscale <- nug.prec.pri[2]

    tmb.post.draws <- 1 / exp(log_nugget_sigma_draws * 2)
    inla.post.draws <- base::sample(x = res_fit$marginals.hyperpar[[names(res_fit$marginals.hyperpar)[grep('nug.id', names(res_fit$marginals.hyperpar))]]][, 1],
                                    size = ndraws,
                                    replace = TRUE, 
                                    res_fit$marginals.hyperpar[[names(res_fit$marginals.hyperpar)[grep('nug.id', names(res_fit$marginals.hyperpar))]]][, 2])
    
    xlim       <- range(c(tmb.post.draws, inla.post.draws, true.val)) 
    
    x.prior <- seq(xlim[1], xlim[2], len = 1000)
    y.prior <- dgamma(x.prior, shape = prior.shape, scale = 1 / prior.iscale)
   
    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log nugget prec"
  }

  ## get posterior samples (we'll use density curves)
  tmb.dens <- density(tmb.post.draws)
  inla.dens <- density(inla.post.draws)
  xrange <- range(c(x.prior, tmb.dens$x, inla.dens$x))
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
  points(tmb.post.draws, rep(0, ndraws), col = alpha(tmb.col, 0.25), cex = 2, pch = '|')
  abline(v = tmb.post.median, col = tmb.col, lwd = 2)
    
  ## plot inla post and data
  lines(inla.dens$x,inla.dens$y, col = inla.col)
  points(inla.post.draws, rep(0, ndraws), col = alpha(inla.col, 0.25), cex = 2, pch = '|')
  abline(v = inla.post.median, col = inla.col, lwd = 2)

  ## plot truth
  abline(v = true.val, col = prior.col, lwd = 2)
  
  ## add legend
  legend("topright", legend = c("prior/truth", "tmb", "inla"), col = c(prior.col, tmb.col, inla.col), lwd = rep(2, 3))
}
dev.off()


## plot results in logit space or in prevalence space?
## plot.in.logit.space <- FALSE ## TODO ?

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
    all.diff.vec <- c( as.vector(true- rtmb), as.vector(true - rinla), as.vector(rtmb - rinla))
    rast.list <- list('TRUE' = true,
                      'TMB' = rtmb,
                      'INLA' = rinla)
    
    png(sprintf('%s/validation/iter%04d_plot_03_median_rasters.png', out.dir, iii),
        height=12, width=12, units = 'in', res = 250)
  }

   if(sum.meas=='stdev'){
    if(data.lik == 'binom'){
      rinla <- ras_sdv_inla_p[[1]]
      rtmb  <- ras_sdv_tmb_p[[1]]
    }else if(data.lik == 'normal'){
      rinla <- ras_sdv_inla[[1]]
      rtmb  <- ras_sdv_tmb[[1]]
    }
    all.vec <- c(as.vector(rtmb), as.vector(rinla))
    all.diff.vec <- as.vector(rtmb - rinla)
    rast.list <- list('TMB' = rtmb,
                      'INLA' = rinla)
    
    png(sprintf('%s/validation/iter%04d_plot_04_stdev_rasters.png', out.dir, iii),
        height=12, width=12, units = 'in', res = 250)
  }

  layout(matrix(1:length(rast.list) ^ 2, byrow = T, ncol = length(rast.list)))
  
  tmp <- subset(dt, period_id==1) ## for s-t

  ## get limits
  rast.zrange <- range(all.vec, na.rm = T)
  diff.zrange <- range(all.diff.vec, na.rm = T)

  for(i in 1:length(rast.list)){
    for(j in 1:length(rast.list)){

      if(i == j){
        ## plot raster
        par(mar = c(0, 0, 1.4, 8), bty='n')
        plot(rast.list[[i]],  maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend.width = 5,
             legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main=paste0(names(rast.list)[i], ': ', sum.meas),
             zlim=rast.zrange)
      }

      if(i < j){
        ## plot scatter
        par(mar = c(4, 4, 2, 2),bty='n')
        plot(x=as.vector(rast.list[[j]])[samp],
             y=as.vector(rast.list[[i]])[samp],
             xlab=names(rast.list)[j],
             ylab=names(rast.list)[i],
             cex=.05, pch=19,
             main=paste0('(sub)SCATTER (', names(rast.list)[i], ' vs ', names(rast.list)[j], '): '))
        lines(x=rast.zrange, y=rast.zrange, col='red')
      }

      if(i > j){
        ## plot difference rasters
        par(mar = c(0, 0, 1.4, 8), bty='n')
        plot(rast.list[[i]] - rast.list[[j]],  maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend.width = 5,
             legend.args=list(text='', side=2, font=1, line=0, cex=0.1),
             main=paste0('Diff: ', names(rast.list)[i], ' - ', names(rast.list)[j]),
             zlim=diff.zrange)
        points( x=tmp$long,y=tmp$lat, pch=19, cex=.1 )

      }
      
    } ## j
  } ## i

  dev.off()
  
} ## sum.meas

## ~~~
## now make caterpillar plots
## ~~~

png(sprintf('%s/validation/iter%04d_plot_05_spatial_re_caterpillars.png', out.dir, iii),
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

if(nrow(plot_d)>2500)
  plot_d <- plot_d[sample(nrow(plot_d),2500,replace=F),]


## ## plot spatial RE differences over time
## ggplot(plot_d, aes(x=tmb_median,y=inla_median,col=period)) + theme_bw() +
## geom_point() + geom_line(aes(group=loc)) + geom_abline(intercept=0,slope=1,col='red') +
## ggtitle('Posterior Medians of Random Effects at Mesh Nodes, TMB v R-INLA. Connected dots same location different periods. ')

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
gg_cat <- ggplot(plot_d, aes(i, tmb_median, col=i)) + theme_bw() + # [seq(1, nrow(plot_d), 5)]
          geom_linerange(aes(ymin = tmb_low, ymax = tmb_up), col='red', size=.8, alpha=.3) +
          geom_linerange(aes(x=i,ymin = inla_low, ymax = inla_up), col='blue', size=.8, alpha=.3) +
          facet_wrap(~period) +
          ggtitle('Comparison of random effects (10% to 90% quantiles) ... BLUE == R-INLA ... RED == TMB')
print(gg_cat)

dev.off()

## ###############################################
## 3) generate and summarize predictive metrics ##  
## ###############################################

## there are two types of predictive metrics we might look at:
##  (i) metrics comparing true surfaces: in latent space and in inv-link space if needed
## (ii) metrics comparing data to estimated surfaces

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## (i) compare true surface to fitted surface
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NOTE: these are in latent space for all models!
all.preds <- data.table(rbind(pred_tmb, pred_inla))

## make a data.table with prediction draws, model type, and truth
## NOTE! the truth and the median.fit are in logit-space!
non.na.idx <- which(!is.na(values(true.rast)))
d <- data.table(truth        = rep(values(true.rast)[non.na.idx], 2),
                model        = c(rep('tmb', nrow(pred_tmb)),
                                 rep('inla', nrow(pred_inla))), 
                median.fit   = apply(all.preds, 1, median))

## get some coverage probs
coverage_probs <- c(25, 50, 80, 90, 95)
for(c in coverage_probs){
  message(paste0('For ',c,'% coverage.'))
  coverage <- c / 100
  li       <- apply(all.preds, 1, quantile, p = (1 - coverage)/2, na.rm=T)
  ui       <- apply(all.preds, 1, quantile, p = coverage + (1 - coverage) / 2, na.rm=T)
  d[, paste0('pix.cov.',c)] = d[['truth']] >= li & d[['truth']] <= ui
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
  for(c in coverage_probs){
    message(paste0('For ',c,'% coverage.'))
    coverage <- c / 100
    li       <- apply(all.preds, 1, quantile, p = (1 - coverage)/2, na.rm=T)
    ui       <- apply(all.preds, 1, quantile, p = coverage + (1 - coverage) / 2, na.rm=T)
    d[, paste0('p.pix.cov.',c)] = d[['truth.p']] >= li & d[['truth.p']] <= ui
  }

  ## get error and pred var
  d[, error.p := truth.p - median.fit.p]
  d[, var.p := apply(all.preds, 1, var)]
}

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

## note iteration
surface.metrics[, iter := iii]

## save
write.csv(surface.metrics, sprintf('%s/validation/iter%04d_surface_metrics.csv',out.dir, iii))

## append into overall metrics for assessing monte carlo variance of metrics
if(iii == 1){
  complete.surface.metrics <- surface.metrics
}else{
  complete.surface.metrics <- rbind(complete.surface.metrics, surface.metrics)
}
