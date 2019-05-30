## #######
## #######
## INLA ##
## #######
## #######

## ########
## SETUP ##
## ########

A <- inla.spde.make.A(
  mesh = mesh_s,
  loc = dt.coords,
  group = dt.pers
)
space   <- inla.spde.make.index("space",
                                n.spde = spde$n.spde,
                                n.group = nperiods)

inla.covs <- covs$name
design_matrix <- data.frame(int = rep(1, nrow(dt)), dt[, inla.covs, with=F], nug.id = 1:nrow(dt))
stack.obs <- inla.stack(tag='est',
                        data=list(Y=dt$Y), ## response
                        A=list(A,1), ## proj matrix, not sure what the 1 is for
                        effects=list(
                          space,
                          design_matrix)
                        )

formula <- formula(paste('Y ~ -1',
                         ifelse(is.null(alpha), '', ' + int'), 
                         ifelse(is.null(betas), '', paste0(' + ', (paste(inla.covs, collapse = ' + ')))),
                         ifelse(is.null(nug.var), '',
                                paste0(' + f(nug.id, model = \'iid\', hyper = list(theta = list(prior = \'loggamma\', param = c(',
                                       nug.prec.pri[1],', ', nug.prec.pri[2], '))))')), 
                         ' + f(space, model = spde, group = space.group, control.group = list(model = \'ar1\'))',
                         sep = ''))

## function to convert from data lik string to integer
## allows easily adding more options even though overkill for just 2
inla.lik.dict <- function(x){
  dict <- list(normal = 'normal',
               binom = 'binomial')
  dict[[x]]
}
inla.setOption("enable.inla.argument.weights", TRUE)



## ######
## FIT ##
## ######

ptm <- proc.time()[3]
if(data.lik == 'normal'){
  res_fit <- inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           ## link = 1, ## removed after looking at NMM
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       prec = list(default = 1 / alphaj.pri[2] ^ 2)),
                  control.inla = list(strategy = inla.approx,
                                      int.strategy = inla.int.strat ##,
                                      ## h = 1e-3, ## removed after looking at NMM
                                      ## tolerance = 1e-6 ## removed after looking at NMM
                                      ),
                  control.compute=list(config = TRUE),
                  control.family = list(hyper = list(prec = list(prior = "loggamma", 
                                                                 param = norm.prec.pri))),
                  family = inla.lik.dict(data.lik),
                  num.threads = cores, #
                  Ntrials = dt$N,
                  scale = dt$N, 
                  ## weights = rep(1, nrow(dt)),
                  verbose = TRUE,
                  keep = FALSE)
}else if(data.lik == 'binom'){
  res_fit <- inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           ## link = 1, ## removed after looking at NMM
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       prec = list(default = 1 / alphaj.pri[2] ^ 2)),
                  control.inla = list(strategy = inla.approx,
                                      int.strategy = inla.int.strat ##,
                                      ## h = 1e-3, ## removed after looking at NMM
                                      ## tolerance = 1e-6 ## removed after looking at NMM
                                      ),
                  control.compute=list(config = TRUE),
                  family = inla.lik.dict(data.lik),
                  num.threads = cores, #
                  Ntrials = dt$N,
                  scale = dt$N, 
                  ## weights = rep(1, nrow(dt)),
                  verbose = TRUE,
                  keep = FALSE)
}
fit_time_inla <- proc.time()[3] - ptm

## check to see if INLA converged nicely
## using the check suggested here: https://groups.google.com/forum/#!topic/r-inla-discussion-group/LsCpuCsr-Qo
## and noting that failing this check may not be terrible
inla.mode.converge <- ifelse(res_fit$mode$mode.status == 0, TRUE, FALSE)

## ##########
## PREDICT ##
## ##########

ptm <- proc.time()[3]
inla_draws <- inla.posterior.sample(ndraws, res_fit, use.improved.mean = bias.correct)
inla_get_draws_time <- proc.time()[3] - ptm

## get parameter names
par_names <- rownames(inla_draws[[1]]$latent)

## index to spatial field and linear coefficient samples
s_idx <- grep('^space.*', par_names)
l_idx <- which(!c(1:length(par_names)) %in% grep('^space.*|Predictor|[*:*]', par_names))

## get spatial draws as matrices and project to deaws at locations 
pred_s <- sapply(inla_draws, function (x) x$latent[s_idx])
pred_inla <- as.matrix(A.pred %*% pred_s)

## get intercept and coef draws and convert to covariate effects
if(!is.null(alpha) | !is.null(betas)){
  pred_l <- sapply(inla_draws, function (x) x$latent[l_idx])
  if(!is.matrix(pred_l)){
    pred_l <- matrix(pred_l, ncol = length(pred_l))
  }
  rownames(pred_l) <- res_fit$names.fixed

  ## extract cell values  from covariates, deal with timevarying covariates here
  non_space_names <- par_names[l_idx]
  cov_effs <- list()
  for(p in 1:nperiods)  cov_effs[[p]] <- cov_vals[[p]][, non_space_names] %*% pred_l
  
  cov_effs <- do.call(rbind, cov_effs)
  
  pred_inla <- pred_inla + cov_effs
}

if(!is.null(nug.var)){
  ## get draws of nugget precision
  pred_n <- sapply(inla_draws, function(x) {
    nug.idx <- which(grepl('nug.id', names(inla_draws[[1]]$hyper)))
    x$hyperpar[[nug.idx]]}) ## this gets the precision for the nugget

  ## simulate nugget noise
  cell_n <- sapply(pred_n, function(x){rnorm(n = nrow(pred_inla),
                                             sd = 1 / sqrt(x),
                                             mean = 0)})
  ## add it on
  pred_inla <- cell_nug + pred_inla
}

## make them into time bins
len = nrow(pred_inla)/nperiods
totalpredict_time_inla <- proc.time()[3] - ptm

## #######
## SAVE ##
## #######

## summarize the latent field
summ_inla <- cbind(median = (apply(pred_inla, 1, median)),
                   sd     = (apply(pred_inla, 1, sd)))

ras_med_inla <- insertRaster(simple_raster, matrix(summ_inla[, 1], ncol = nperiods))
ras_sdv_inla <- insertRaster(simple_raster, matrix(summ_inla[, 2], ncol = nperiods))

saveRDS(file = sprintf('%s/modeling/outputs/inla/experiment%04d_iter%04d_inla_preds_median_raster.rds', out.dir, par.iter, iii), object = ras_med_inla)
saveRDS(file = sprintf('%s/modeling/outputs/inla/experiment%04d_iter%04d_inla_preds_stdev_raster.rds', out.dir, par.iter, iii), object = ras_sdv_inla)

if(data.lik == 'binom'){
  ## convert to prevalence space and summarize, rasterize, and save again
  pred_inla_p <- plogis(pred_inla)
  
  summ_inla_p <- cbind(median = (apply(pred_inla_p, 1, median)),
                       sd     = (apply(pred_inla_p, 1, sd)))

  ras_med_inla_p <- insertRaster(simple_raster, matrix(summ_inla_p[, 1], ncol = nperiods))
  ras_sdv_inla_p <- insertRaster(simple_raster, matrix(summ_inla_p[, 2], ncol = nperiods))

  saveRDS(file = sprintf('%s/modeling/outputs/inla/experiment%04d_iter%04d_inla_preds_median_raster_PREV.rds', out.dir, par.iter, iii),
          object = ras_med_inla_p)
  saveRDS(file = sprintf('%s/modeling/outputs/inla/experiment%04d_iter%04d_inla_preds_stdev_raster_PREV.rds', out.dir, par.iter, iii),
          object = ras_sdv_inla_p)
}