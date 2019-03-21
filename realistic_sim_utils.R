## some functions to help with realistic inla/tmb simulation comparison
## written by a0z 5/17/18

## qsub_sim: function to launch sims (and sim comparisons) on the cluster
qsub_sim <- function(iter, ## if looping through muliple models, used to give different names
                     main.dir.nm, ## head dir to store all results
                     codepath,
                     slots, 
                     singularity = 'default',
                     singularity_opts = NULL,
                     extra_name = '', 
                     logloc = NULL ## defaults to input/output dir in main.dir/iter/
                     ){

  ## some early checks
  if(length(unique(c(length(cov_names), length(cov_measures), length(betas)))) != 1){
    messge('cov_names, cov_measrures, and betas lengths do not match! fix and rerun')
    stop()
  }

  ## if(sp.range < 0 | sp.var < 0 | sp.alpha < 0 | nug.var < 0){
  ##   message('sp range, var, alpha, or nugget var is not positive. fix and rerun!')
  ##   stop()
  ## }

  ## setup some stuff that I don't think I'll want to change
  proj <- 'proj_geo_nodes'
  node.flag <- '-l geos_node=TRUE'
  shell <- '/share/code/geospatial/lbd_core/mbg_central/share_scripts/shell_sing.sh'
  sing_image <- get_singularity(image = singularity)
  if(is.null(singularity_opts)) singularity_opts <- list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)
  if(is.null(logloc)) logloc <- sprintf('/homes/azimmer/tmb_inla_sim/%s/%s/logs', main.dir.nm, iter)
  
  ## Piece together lengthy `qsub` command
  qsub <- paste0("qsub",
                 " -e ", logloc, "/errors/",
                 " -o ", logloc, "/output/",
                 " -pe multi_slot ", slots,
                 " -P ", proj, " ", node.flag)

  ## add on stuff to launch singularity
  qsub <- qsub_sing_envs(qsub, singularity_opts,
                         sing_image)

  ## append job name, shell, and code to run 
  qsub <- paste0(qsub,
                 sprintf(" -N sim_job_%s_%i", extra_name, iter), ## job name
                 " ", shell, " ", codepath) ## shell and code path

  ## add on all remaining arguments 
  qsub <- paste(qsub,
                iter, ## which row in loopvars
                main.dir.nm, ## which dir to load from
                sep = " ")

  return(qsub)
  
}

## a function to simulate from a GF using INLA
## this function is taken from the spde tutorial
rspde <- function (coords, kappa, variance = 1, alpha = 2, n = 1, mesh,
                   verbose = FALSE, seed, return.attributes = FALSE)
{

  if(is.null(seed)) seed = 0 ## If ‘seed=0L’ then GMRFLib will set the seed
                             ## intelligently/at 'random'
  
    t0 <- Sys.time()
    theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
    if (verbose)
        cat("theta =", theta, "\n")
    if (missing(mesh)) {
        mesh.pars <- c(0.5, 1, 0.1, 0.5, 1) * sqrt(alpha - ncol(coords)/2)/kappa
        if (verbose)
            cat("mesh.pars =", mesh.pars, "\n")
        attributes <- list(mesh = inla.mesh.2d(, coords[chull(coords),
                                                        ],
                                               max.edge = mesh.pars[1:2],
                                               cutoff = mesh.pars[3],
                                               offset = mesh.pars[4:5]))
        if (verbose)
            cat("n.mesh =", attributes$mesh$n, "\n")
    }
    else attributes <- list(mesh = mesh)
    attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
    attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
    attributes$A <- inla.mesh.project(mesh = attributes$mesh,
                                      loc = coords)$A
    if (n == 1)
        result <- drop(attributes$A %*% inla.qsample(Q = attributes$Q,
                                                     constr = attributes$spde$f$extraconstr))
    t1 <- Sys.time()
    result <- inla.qsample(n, attributes$Q, seed = ifelse(missing(seed),
                                                          0, seed),
                           constr = attributes$spde$f$extraconstr)
    if (nrow(result) < nrow(attributes$A)) {
        result <- rbind(result, matrix(NA, nrow(attributes$A) -
                                           nrow(result), ncol(result)))
        dimnames(result)[[1]] <- paste("x", 1:nrow(result), sep = "")
        for (j in 1:ncol(result)) result[, j] <- drop(attributes$A %*%
                                                      result[1:ncol(attributes$A),
                                                             j])
    }
    else {
        for (j in 1:ncol(result)) result[1:nrow(attributes$A),
                                         j] <- drop(attributes$A %*% result[, j])
        result <- result[1:nrow(attributes$A), ]
    }
    t2 <- Sys.time()
    attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2 -
                                                              t0)
    if (return.attributes)
        attributes(result) <- c(attributes(result), attributes)
    return(drop(result))
}


sim.realistic.data <- function(reg,
                               year_list,
                               data.lik, ## either 'binom' or 'normal'
                               sd.norm = NULL, ## sd of normal observations 
                               betas = NULL, ## if null, use no covs
                               sp.kappa,
                               sp.alpha,
                               t.rho,
                               nug.var = NULL, 
                               n.clust,
                               m.clust,
                               covs = NULL,
                               cov_layers = NULL, ## if supplied, use this instead of reloading covs
                               simple_raster,
                               simple_polygon, 
                               out.dir,
                               pop_raster = NULL,
                               obs.loc.strat = 'rand', ## either 'rand' or 'pop.strat'. NOTE: random is proportional to population!
                               urban.pop.pct = 1, ## percent (in percent = alpha*100% space - i.e. 1 for 1%) of population that comprises urban
                               urban.strat.pct = 40, ## percent of sample locations that should come from urban pixels
                               sp.field.sim.strat = 'RF', ## one of RF or SPDE ## TODO add t-dist, extremal
                               seed = NULL,
                               exp.iter = 1){ ## exp.iter, used for file saving

  ## make some checks and set things
  if(is.null(pop_raster) & obs.loc.strat == 'pop.strat'){
    stop("You need to supply a pop raster to use obs.loc.strat=='pop.strat'")
  }
  if(obs.loc.strat == 'pop.strat' & (is.null(urban.pop.pct) | is.null(urban.strat.pct))){
      stop("You must pass in urban.pop.pct and urban.strat.pct values when you set obs.loc.strat=='pop.strat'")
  }
  if(data.lik == 'normal' & is.null(sd.norm)){
        stop("You must pass in sd.norm if setting data.lik=='normal'")
  }
  if(!is.null(betas)){
    if(!is.null(cov_layers)){
      if(length(cov_layers) != length(betas)){
        stop('The supplied cov_layers object does not match the beta arg in length. Something has gone wrong')
      }
    }else if(!is.null(covs)){
      if(nrow(covs) != length(betas)){
        stop('nrow(covs) does not match the beta arg in length. Something has gone wrong')
      }
    }else{
      stop('You must pass in either but no covs of cov_layers when you pass in betas.')
    }
  }

  ## set seed if required
  if(!is.null(seed)) set.seed(seed)
  
  ## create dir for simulated objects
  dir.create(sprintf('%s/simulated_obj/', out.dir), recursive = T, showWarnings = F)
  
  ## #####################################
  ## load and prepare covariate rasters ##
  ## #####################################
  
##  if(!is.null(betas)){
    if(is.null(cov_layers)){
      message('\n\nLOADING COVS\n')
      cov_layers <- load_and_crop_covariates_annual(covs            = covs[, name],
                                                    measures        = covs[, meas],
                                                    simple_polygon  = simple_polygon,
                                                    start_year      = min(year_list),
                                                    end_year        = max(year_list),
                                                    interval_mo     = 12) ## always grab annual, then subset if need be
      ## subset to years in yearlist
      for(cc in 1:length(cov_layers)){
        if(dim(cov_layers[[cc]])[3] > 1){
          cov_layers[[cc]] <- cov_layers[[cc]][[which( min(year_list):max(year_list) %in% year_list )]]
        }
        
        ## center-scale covs. (TODO by year or across years?)
        cov_layers[[cc]] <- (cov_layers[[cc]] - mean(values(cov_layers[[cc]]), na.rm = T)) / sd(values(cov_layers[[cc]]), na.rm = T)
      }
      
      ## we also want our cov_layers to align with simple_raster
      for(l in 1:length(cov_layers)) {
        message(sprintf("On cov %i out of %i", l, length(cov_layers)))
        cov_layers[[l]]  <- crop(cov_layers[[l]], extent(simple_raster))
        cov_layers[[l]]  <- setExtent(cov_layers[[l]], simple_raster)
        cov_layers[[l]]  <- mask(cov_layers[[l]], simple_raster)
      }
    }else{
      message('\n\nUSING PRE-SUPPLIED AND PREPPED COVS\n')
    }  
    
    ## plot center-scaled covariates
    pdf(sprintf('%s/simulated_obj/cov_plot.pdf', out.dir), width = 16, height = 16)
    for(cc in 1:length(cov_layers)){
      message(sprintf('Plotting covariate: %s\n', names(cov_layers)[[cc]]))
      if(dim(cov_layers[[cc]])[3] == 1){
        message('--plotting single synoptic map\n')
        par(mfrow = c(1, 1))
      }else{
        message('--plotting time series\n')
        par(mfrow = rep( ceiling(sqrt( dim(cov_layers[[cc]])[3] )), 2))
      }
      for(yy in 1:dim(cov_layers[[cc]])[3]){
        raster::plot(cov_layers[[cc]][[yy]],
                     main = paste(names(cov_layers)[cc],
                                  ifelse(dim(cov_layers[[cc]])[3] == 1,
                                         'synoptic',
                                         year_list[yy]), 
                                  sep = ': '))
      }
    }
    dev.off()
##  }else{
##    message('\n\nNO COVS\n')
##  }
  
  ## now we can simulate our true surface
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## ##########################
  ## simulate space-time gp ##
  ## ##########################
  
  message('SIMULATE GP\n')
  
  ## FIRST, get the pixel coords- these are useful later too
  
  ## convert simple raster of our region to spatialpolygonsDF
  pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)
  
  ## reproject sp obj
  geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
  pix.pts <- spTransform(pix.pts, CRS(geo.prj)) 
  proj4string(pix.pts)
  
  ## to simulate, we need lat-lon locs for the entire raster
  ## get coords
  pix.pts@data <- data.frame(pix.pts@data, long=coordinates(pix.pts)[,1],
                             lat=coordinates(pix.pts)[,2])
  pix.pts.numeric <- as.data.frame(pix.pts@data)
  
  ## sim using SPDE mesh ## TODO export this mesh to use in fitting
  if(sp.field.sim.strat == 'SPDE'){ 
    
    ## now we can use these coords to simulate GP from rspde()
    reg.mesh <- inla.mesh.2d(boundary = inla.sp2segment(simple_polygon),
                             loc = pix.pts@data[, 2:3],
                             max.edge = c(0.25, 5),
                             offset = c(1, 5),
                             cutoff = 0.25)
    
    
    ## get spatial fields that are GPs across space and are indep in time
    sf.iid <- rspde(coords = cbind(pix.pts.numeric[, 2], pix.pts.numeric[, 3]),
                    kappa = sp.kappa,
                    variance = sp.var,
                    alpha = sp.alpha,
                    mesh = reg.mesh,
                    n = length(year_list),
                    seed = seed)
  }else{
    reg.mesh <- NULL
  }
  
  ## use random fields package on simple_raster to simulate GP for spatial field
  if(sp.field.sim.strat == 'RF'){ 
    model <- RMmatern(nu = sp.alpha - 1, ## from INLA book
                      scale = sqrt(2 * (sp.alpha - 1)) / sp.kappa, 
                      var = sp.var)      
    ## sf.iid <- geostatsp::RFsimulate(model, x = simple_raster, n = length(year_list))
    sf.iid <- RFsimulate(model, x = pix.pts.numeric[, 2], y = pix.pts.numeric[, 3], n = length(year_list), spConform = FALSE)
  }
  
  ## simulate t dist with low DOF
  if(sp.field.sim.strat == 't'){
    message("sp.field.sim.strat=='t' is not yet implemented")
  }
  
  ## simulate extremal dist
  if(sp.field.sim.strat == 'ext'){
    message("sp.field.sim.strat=='ext' is not yet implemented")
  }
  
  ## ---------
  ## introduce temporal ar1 correlation at the pixel level
  if(length(year_list) > 1){ ## then, correlate gp draws
    sf.cor <- sf.iid
    for(ii in 2:ncol(sf.cor)){
      sf.cor[, ii] <- t.rho * sf.cor[, ii - 1] + sqrt(1 - t.rho ^ 2) * sf.iid[, ii]
    }
    
    
    ## convert them to rasters
    for(cc in 1:ncol(sf.iid)){
      if(cc == 1){
        sf.rast <- rasterize(x = pix.pts@data[, 2:3],
                             y = simple_raster,
                             field = sf.cor[, cc])
      }else{
        sf.rast <- addLayer(sf.rast,
                            rasterize(x = pix.pts@data[, 2:3],
                                      y = simple_raster,
                                      field = sf.cor[, cc]))
      }
    }
  }else{ ## we have a single year, no time corr needed
    sf.rast <- rasterize(x = pix.pts@data[, 2:3],
                         y = simple_raster,
                         field = sf.iid)
  }
  
  
  ## plot gp
  pdf(sprintf('%s/simulated_obj/iter%04d_st_gp_plot.pdf', out.dir, exp.iter), width = 16, height = 16)
  par(mfrow = rep( ceiling(sqrt( dim(sf.rast)[3] )), 2))
  for(yy in 1:dim(sf.rast)[3]){
    raster::plot(sf.rast[[yy]],
                 main = paste('GP',
                              year_list[yy],
                              sep = ': '))
  }
  dev.off()
  
  ## TODO alternatively, simulate non-Gaussian field!
  ## using RandomFields, e.g.
  ## model <- RPopitz(RMexp(), alpha=2)
  ## z1 <- RFsimulate(model, seq(0,2,0.01),seq(0,2,0.01),grid=T)
  ## plot(z1, type = 'l')
  
  ## now we can make the iid nugget
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## simulate IID normal draws to add as nugget
  
  if(!is.null(nug.var) & data.lik != 'normal'){
    ## normal plus nugget means add nugget in normal draws, not to every pixel!
    ## TODO is this right to set it up and limit it this way??
    
    message('SIMULATE PIXEL NUGGET\n')
    
    ## take rnorm() draws and convert them to rasters
    for(cc in 1:length(year_list)){
      if(cc == 1){ ## initialize
        nug.rast <- rasterize(x = pix.pts@data[, 2:3],
                              y = simple_raster,
                              field = rnorm(n = nrow(pix.pts@data),
                                            mean = 0,
                                            sd = sqrt(nug.var)))
      }else{ ## add a layer
        nug.rast <- addLayer(nug.rast,
                             rasterize(x = pix.pts@data[, 2:3],
                                       y = simple_raster,
                                       field = rnorm(n = nrow(pix.pts@data),
                                                     mean = 0,
                                                     sd = sqrt(nug.var))))        
      } ## else
    } ## year loop
    
    ## plot nugget
    pdf(sprintf('%s/simulated_obj/nugget_plot.pdf', out.dir), width = 16, height = 16)
    par(mfrow = rep( ceiling(sqrt( dim(sf.rast)[3] )), 2))
    for(yy in 1:dim(nug.rast)[3]){
      raster::plot(nug.rast[[yy]],
                   main = paste('GP',
                                year_list[yy],
                                sep = ': ')) 
    }
    dev.off()
    
  }else{
    message('NO NUGGET\n')
  }
  
  
  ## now we can make the true surface
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## #####################################################
  ## make true surface by combining cov effects and gp ##
  ## #####################################################
  
  ## finally, we combine the gp and the covariate effects to get our surface in link (e.g. logit if binomial) space
  
  ## start with spatial field
  true.rast <- sf.rast
  
  ## add on covs
  if(!is.null(betas)){
    for(cc in 1:length(cov_layers)){ ## loop though and add on coefficients*covariates to gp raster layers
      true.rast <- true.rast + betas[cc] * cov_layers[[cc]] ## should work for both stationary and time-varying
    }
  }
  
  ## we append the gp to the cov_layers
  ## TODO this seems like a bad idea... adjust this to keep gp out of covs 
  cov_layers[['gp']] <- sf.rast
  
  ## and, add nugget if desired
  if(!is.null(nug.var)) true.rast <- true.rast + nug.rast
  
  pdf(sprintf('%s/simulated_obj/iter%04d_true_surface_plot.pdf', out.dir, exp.iter), width = 16, height = 16)
  par(mfrow = rep( ceiling(sqrt( dim(true.rast)[3] )), 2))
  for(yy in 1:dim(true.rast)[3]){
    raster::plot(true.rast[[yy]],
                 main = paste('GP',
                              year_list[yy], 
                              sep = ': '))
  }
  dev.off()
  
  saveRDS(sprintf('%s/simulated_obj/iter%04d_true_surface_raster.rds', out.dir, exp.iter), object = true.rast)
  

  ## now the surface simulation is done and all we need to do is simulate the data

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #####################################
  ## simulate data from true surface ##
  #####################################

  message('SIMULATE DATA\n')
  
  ## randomly (for now) select data boservation locations across time
  
  ## to do this, we sample, with replacement, from the lat-longs that we used to sim the GP
  if(obs.loc.strat == 'rand'){ ## select locations totally at random
    sim.rows <- sample(x = 1:nrow(pix.pts.numeric), prob = pix.pts.numeric[, 1],
                       size = n.clust * length(year_list), replace = TRUE)
  } else{ ## stratify by "urban"/"rural"
    ## given the % of population you want to be urban, find the population value cutoff
    urban_thresh <- quantile(probs = (1 - urban.pop.pct), na.omit(values(pop_raster)))
    
    ## make a binary urban rural raster and get the lat-longs of the pixels
    u_r_raster <- pop_raster[[1]] ## urban is 1, rural is 0
    u_r_raster[pop_raster[[1]] < urban_thresh] <- 0
    u_r_raster[pop_raster[[1]] >= urban_thresh] <- 1
    
    ## convert pixels to a data frame
    u_r.pts <- rasterToPoints(u_r_raster, spatial = TRUE)
    u_r.pts@data <- data.frame(u_r.pts@data, long=coordinates(u_r.pts)[,1],
                               lat=coordinates(u_r.pts)[,2])
    u_r.pts.numeric <- as.data.frame(u_r.pts@data)
    
    ## sample stratified locations
    u.rows <- sample(x = which(u_r.pts.numeric[, 1] == 1), size = round(n.clust * urban.strat.pct),
                     replace = TRUE)
    r.rows <- sample(x = which(u_r.pts.numeric[, 1] == 0), size = round(n.clust * (1 - urban.strat.pct)),
                     replace = TRUE)
    sim.rows <- c(u.rows, r.rows)
  }
  
  ## generate a table of simulated data at the selected locations
  sim.dat <- as.data.table(pix.pts.numeric[, -1])
  sim.dat <- sim.dat[sim.rows, ]
  
  ## add in years
  sim.dat[, year := rep(year_list, each = n.clust)]
  
  ## extract the value of the true surface at data locations
  true_p_logit <- numeric(nrow(sim.dat))
  for(yy in unique(year_list)){
    true_p_logit[which(sim.dat[, year] == yy)] <- raster::extract(x = true.rast[[ which(year_list %in% yy) ]],
                                                                  y = sim.dat[year == yy, .(long, lat)])
  }
  
  ## sim sample size of observations
  sim.dat[, N := rpois(n = nrow(sim.dat), lambda = m.clust)]
  
  if(data.lik == 'binom'){
    sim.dat[, p_true := inv.logit(true_p_logit)]
    
    ## and now we simulate binomial observations from the true surface
    sim.dat[, Y := rbinom(n = nrow(sim.dat), size = sim.dat[, N], prob = sim.dat[, p_true])]
    
    ## and get empirical p_obs
    sim.dat[, p_obs := Y / N]
    
    ## lastly, we tag on weight. TODO flesh out weighting options
    sim.dat[, weight := 1]
  }else if(data.lik == 'normal'){
    sim.dat[, p_true := true_p_logit] ## no link to transform
    
    ## and now we simulate binomial observations from the true surface
    sim.dat[, Y := apply(.SD, 1, function(x) mean(rnorm(n = x[1], mean = x[2], sd = sd.norm))), .SDcols = c('N', 'p_true')]
    
    ## and get empirical p_obs
    sim.dat[, p_obs := Y]
    
    ## lastly, we tag on weight. TODO flesh out weighting options
  }
  
  ## now we just finish by making some convenience objects and saving everything
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## #########################################################################################
  ## for convenience, we also extract covariate values to the same df (and the true gp val) ##
  ## ##########################################################################################
  
  message('PREPARE AND SAVE OBJECTS\n')
  
  cov.mat <- matrix(ncol = length(cov_layers),
                    nrow = nrow(sim.dat))
  for( cc in 1:length(cov_layers) ){
    tmp <- numeric(nrow(sim.dat))
    
    if(dim(cov_layers[[cc]])[3] > 1){ ## check for time-varying
      for(ll in 1:dim(cov_layers[[cc]])[3]){
        tmp[which(sim.dat[, year] == year_list[ll])] <- raster::extract(x = cov_layers[[cc]],
                                                                        y = sim.dat[year == year_list[ll], .(long, lat)],
                                                                        layer = ll)[, 1]
      }
    }else{ ## space only rasters
      tmp <- raster::extract(x = cov_layers[[cc]], y = sim.dat[, .(long, lat)])
    }
    
    cov.mat[, cc] <- tmp  
  }
  cov.mat <- as.data.table(cov.mat)
  setnames(cov.mat, names(cov_layers))
  
  ## ###################################
  ## combine into a single master df ##
  ## ###################################
  sim.dat <- cbind(sim.dat, cov.mat)
  
  
  ## #################################
  ## save everything we might want ##
  ## #################################
  saveRDS(object = sim.dat,
          file = sprintf('%s/simulated_obj/iter%04d_sim_data.rds', out.dir, exp.iter))
  
  saveRDS(object = cov_layers,
          file = sprintf('%s/simulated_obj/iter%04d_cov_gp_rasters.rds', out.dir, exp.iter))

  if(sp.field.sim.strat == 'SPDE'){
    saveRDS(object = reg.mesh,
            file = sprintf('%s/simulated_obj/iter%04d_region_mesh.rds', out.dir, exp.iter))
  }
  
  #########################
  ## return a named list ##
  #########################
  return(list(sim.dat = sim.dat,
              cov.gp.rasters = cov_layers,
              true.rast = true.rast, 
              mesh_s = reg.mesh))
}



#######
## Helper function for turning an xyzt table into a raster
rasterFromXYZT <- function(table,
                           z,t){
  require(data.table)
  n_periods = length(unique(table[,t]))
  table$t=table[,t]
  res=  stack(rasterFromXYZ(as.matrix(table[t==1,c('x','y',z),with=F])))
  if(n_periods>1)
    for(r in 2:n_periods)
      res=addLayer(res, rasterFromXYZ(as.matrix(table[t==r,c('x','y',z),with=F])))
  return(res)
}
