## some functions to help with realistic inla/tmb simulation comparison
## written by a0z 5/17/18

## qsub_sim: function to launch sims (and sim comparisons) on the cluster
qsub_sim <- function(exp.lvid, ## if looping through multiple experiments - i.e. row of loopvars
                     exp.iter, ## if monte carlo iteration within an experiment
                     exp.hash, ## unique 6char string to identify all jobs belonging to the same loopvar submission
                     main.dir.nm, ## head dir to store all results
                     codepath,
                     singularity = 'default',
                     singularity_opts = list(SET_OMP_THREADS = cores, SET_MKL_THREADS = cores),
                     extra_name = '',
                     mem = '20G',
                     time = '01:00:00:00', 
                     queue = 'geospatial.q',
                     priority = 0, ## defaults 0, can be as low as -1023 to reduce priority
                     hold.jid = NULL, ## jobid to hold on. default no hold
                     logloc = NULL ## defaults to input/output dir in main.dir/exp.lvid/
                     ){

  ## some early checks
  if(length(unique(c(length(cov_names), length(cov_measures), length(betas)))) != 1){
    messge('cov_names, cov_measrures, and betas lengths do not match! fix and rerun')
    stop()
  }

  ## set correct project based on queue
  proj <- ifelse(queue=='geospatial.q', 'proj_geo_nodes', 'proj_geospatial')
  
  ## make sure we have access to J if running on non geo nodes
  node.flag <- ifelse(queue != 'geospatial.q', ' -l archive=TRUE ', '')
  
  ## grab the shell script we want
  shell <- '/share/code/geospatial/azimmer/lbd_core/mbg_central/share_scripts/shell_sing.sh'
  sing_image <- get_singularity(image = singularity)
  
  ## set the loglocation for output/error files
  if(is.null(logloc)){
    logloc <- sprintf('/ihme/scratch/users/azimmer/tmb_inla_sim/%s/%04d/logs', main.dir.nm, exp.lvid)
  }
  error_log_dir <- paste0(logloc, '/errors/')
  output_log_dir <- paste0(logloc, '/output/')
  
  ## Piece together lengthy `qsub` command
  qsub <- paste0("qsub",
                 " -e ", logloc, "/errors/",
                 " -o ", logloc, "/output/",
                 " -q ", queue, 
                 " -P ", proj,
                 " -p ", priority)

  ## add on job resource requests
  qsub <- paste0(qsub,
                ' -l m_mem_free=', mem,
                ' -l fthread=1',
                ' -l h_rt=', time, 
                node.flag ## either ' -l geos_node=TRUE', or ''
  )
  
  ## add on stuff to launch singularity
  qsub <- qsub_sing_envs(qsub, singularity_opts,
                         sing_image)
  
  ## add hold on jobid flag
  if(!is.null(hold.jid)) {
    qsub <- paste(qsub, 
                  "-hold_jid", 
                  hold.jid)
  }
  
  ## append job name, shell, and code to run 
  qsub <- paste0(qsub,
                 sprintf(" -N sim_job_%s_hash_%s_exp%04d_iter%04d", 
                         extra_name, exp.hash, exp.lvid, exp.iter), ## job name
                 " ", shell, " ", codepath) ## shell and code path

  ## add on all remaining arguments 
  qsub <- paste(qsub,
                exp.lvid, ## which row in loopvars
                exp.iter, ## which iteration of experiment
                main.dir.nm, ## which dir to load from
                sep = " ")

  return(qsub)
}

## used to monitor a set of jobs given amatrix of job.ids
track.exp.iter <- function(jid.dt, main.dir) {
  
  ## query the system to see all jobs
  ## first two rows are fluff
  q.ret <- system("qstat", intern = TRUE)
  
  ## get all jobids and their status in a nice data.table
  q.jobs <- do.call('rbind',
                    lapply(q.ret[-(1:2)], 
                    FUN = function(x) { q.ret.spl <- strsplit(x, split = ' ')[[1]];
                                       data.table(jid    = q.ret.spl[3],
                                                  j.status = q.ret.spl[12])
                                     }))
  
  ## merge status onto jid.dt
  jid.dt <- merge(jid.dt, q.jobs, by ='jid', all.x=T, all.y=F)
  
  ## check the csvs
  jt.dir <- paste(main.dir, 'common', 'job_tracking', sep = '/')
  jt.csvs <- list.files(path = jt.dir)
  if(length(jt.csvs) == 0) { ## no files yet, 
    jt.info <- data.table(exp     = character(),
                          iter    = character(),
                          sim_loop_ct = numeric(), ## sim.loop.ct from 1_run_space_sim.R
                          script_num  = numeric())
  } else { 
    jt.info <- do.call('rbind',
                       lapply(jt.csvs,
                              FUN = function(x) {
                                csv.spl <- strsplit(x, split='.')[[1]];
                                csv.spl <- strsplit(x, split='_')[[1]];
                                csv.dat <- read.csv(paste(jt.dir, x, sep='/'));
                                data.table(exp     = csv.spl[2],
                                           iter    = substr(csv.spl[4], start=1, stop=4), ## the status is appended, use last row
                                           sim_loop_ct = csv.dat[nrow(csv.dat), 1], ## sim.loop.ct from 1_run_space_sim.R
                                           script_num  = csv.dat[nrow(csv.dat), 2]) ## script number
                              } ## func
                       ) ## lapply
    ) ## do.call
  } ## else: length(jt.csvs) > 0
  
  ## merge on csv jt info
  jid.dt <-  merge(jid.dt, jt.info, by = c('exp', 'iter'), all.x=T)
  
  ## determine status (not started, running, failed, completed) for each job
  
  ## number of jobs in queue
  jid.dt[(grepl('qw', j.status) | grepl('r', j.status)), in_q := 1]
  ## not started means job in qw
  jid.dt[grepl('qw', j.status), not_started := 1]
  ## running
  jid.dt[grepl('r', j.status), running := 1]
  ## errored
  jid.dt[(is.na(j.status) & (script_num != 0 | is.na(script_num))), errored := 1]
  ## not errored
  jid.dt[is.na(errored), on_track := 1]
  ## completed
  jid.dt[(is.na(j.status) & script_num==0), completed := 1]
  
  ## swap NA for 0 for averaging
  jid.dt[is.na(jid.dt)] <- 0
  
  ## count number of iterations in q by experiment
  exp.sum <- jid.dt[, lapply(.SD, sum), by=exp, .SDcols=c('in_q')]
  
  ## take the mean of running status
  exp.sum <- merge(exp.sum, 
                   jid.dt[, lapply(.SD, FUN=function(x){round(100*mean(x, na.rm=T), 2)}), 
                             by=exp, 
                    .SDcols=c('on_track', 'not_started', 'running', 
                              'errored', 'completed')],
                   all.x=T, all.y=F, by='exp')
  
  return(list(full.tracker = jid.dt,
              summ.tracker = exp.sum))
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
                               pixel.iid.var = NULL, ## pixel iid variance (spatial discontinuity)
                               n.clust,
                               m.clust,
                               clust.re.var = NULL, ## iid RE variance for each cluster observed 
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
                               verbose = FALSE,
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
      if(verbose) message('\n\nLOADING COVS\n')
      cov_layers <- load_and_crop_covariates_annual(covs            = covs[, name],
                                                    measures        = covs[, meas],
                                                    simple_polygon  = simple_polygon,
                                                    start_year      = min(year_list),
                                                    end_year        = max(year_list),
                                                    interval_mo     = 12) ## always grab annual, then subset if need be
      
      ## loop through covs, subset to years, align with simple raster, center-scale
      for(cc in 1:length(cov_layers)) {
        
        ## subset to years in yearlist
        if(dim(cov_layers[[cc]])[3] > 1){
          cov_layers[[cc]] <- cov_layers[[cc]][[which( min(year_list):max(year_list) %in% year_list )]]
        }
        
        ## align covs with simple_raster
        if(verbose) message(sprintf("On cov %i out of %i", cc, length(cov_layers)))
        cov_layers[[cc]]  <- crop(cov_layers[[cc]], extent(simple_raster))
        cov_layers[[cc]]  <- setExtent(cov_layers[[cc]], simple_raster)
        cov_layers[[cc]]  <- mask(cov_layers[[cc]], simple_raster)
        
        ## center-scale covs. (TODO by year or across years?)
        cov_layers[[cc]] <- (cov_layers[[cc]] - mean(values(cov_layers[[cc]]), na.rm = T)) / sd(values(cov_layers[[cc]]), na.rm = T)
      }
      
    
      
    }else{
      if(verbose) message('\n\nUSING PRE-SUPPLIED AND PREPPED COVS\n')
    }  
    
    ## plot center-scaled covariates
    pdf(sprintf('%s/simulated_obj/cov_plot.pdf', out.dir), width = 16, height = 16)
    for(cc in 1:length(cov_layers)){
      if(verbose) message(sprintf('Plotting covariate: %s\n', names(cov_layers)[[cc]]))
      if(dim(cov_layers[[cc]])[3] == 1){
        if(verbose) message('--plotting single synoptic map\n')
        par(mfrow = c(1, 1))
      }else{
        if(verbose) message('--plotting time series\n')
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
##    if(verbose) message('\n\nNO COVS\n')
##  }
  
  ## now we can simulate our true surface
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## ##########################
  ## simulate space-time gp ##
  ## ##########################
  
  if(verbose) message('SIMULATE SPATIAL FIELD\n')
  
  ## FIRST, get the pixel coords- these are useful later too
  
  ## convert simple raster of our region to spatialpolygonsDF
  ## these are only the NON-NA pixels!
  pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)
  
  ## reproject sp obj to default used in covariates
  geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
  pix.pts <- spTransform(pix.pts, CRS(geo.prj)) 
  ## proj4string(pix.pts) ## check proj is what we wanted
  
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
                             max.edge = c(0.1, 5),
                             offset = c(1, 5),
                             cutoff = 0.1)
    # FOR PAPER FIGURE:  
    # png('~/tmb_inla_paper/spatial_matern_varying_ranges.png', width=8, height=8, units='in', res=300)
    # set.seed(413)
    # par(mfrow=c(2, 2), mai=c(.5, .5, .5, .85))
    # for(sp.range in c(.5, 1, sqrt(8), 5)){
    #   sp.kappa = sqrt(8)/sp.range
    
    ## get spatial fields that are GPs across space and are indep in time
    sf.iid <- rspde(coords = cbind(pix.pts.numeric[, 2], pix.pts.numeric[, 3]),
                    kappa = sp.kappa,
                    variance = sp.var,
                    alpha = sp.alpha,
                    mesh = reg.mesh,
                    n = length(year_list),
                    seed = seed)
    
    # FOR PAPER FIGURE:  
    #   ## rasterize and plot (for testing params)
    #   sf.rast <- rasterize(x = pix.pts@data[, 2:3],
    #                        y = simple_raster_full,
    #                        field = sf.iid)
    #   if(sp.range==.5) plot(sf.rast, col=viridis(256),
    #                         main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = 0.5')))
    #   if(sp.range==1) plot(sf.rast, col=viridis(256),
    #                        main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = 1')))
    #   if(sp.range==sqrt(8)) plot(sf.rast, col=viridis(256),
    #                              main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = ', sqrt(8))))
    #   if(sp.range==5) plot(sf.rast, col=viridis(256),
    #                        main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = 5')))
    # }
    # dev.off()
    
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
    stop("sp.field.sim.strat=='t' is not yet implemented")
  }
  
  ## simulate extremal dist
  if(sp.field.sim.strat == 'ext'){
    stop("sp.field.sim.strat=='ext' is not yet implemented")
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
  
  if(!is.null(pixel.iid.var)){ #} & data.lik != 'normal'){
    ## normal plus nugget means add nugget in normal draws, not to every pixel!
    ## TODO is this right to set it up and limit it this way??
    
    if(verbose) message('--simulate discontinuous pixel RE\n')
    
    ## take rnorm() draws and convert them to rasters
    for(cc in 1:length(year_list)){
      if(cc == 1){ ## initialize
        nug.rast <- rasterize(x = pix.pts@data[, 2:3],
                              y = simple_raster,
                              field = rnorm(n = nrow(pix.pts@data),
                                            mean = 0,
                                            sd = sqrt(pixel.iid.var)))
      }else{ ## add a layer
        nug.rast <- addLayer(nug.rast,
                             rasterize(x = pix.pts@data[, 2:3],
                                       y = simple_raster,
                                       field = rnorm(n = nrow(pix.pts@data),
                                                     mean = 0,
                                                     sd = sqrt(pixel.iid.var))))        
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
    if(verbose) message('--no spatial discontinuity RE\n')
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
  if(!is.null(pixel.iid.var)) true.rast <- true.rast + nug.rast
  
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

  if(verbose) message('SIMULATE DATA\n')
  
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
  ## this is the true value of the linear predictor! not necessarily in logit=space...
  true_lin_pred <- numeric(nrow(sim.dat))
  for(yy in unique(year_list)){
    true_lin_pred[which(sim.dat[, year] == yy)] <- raster::extract(x = true.rast[[ which(year_list %in% yy) ]],
                                                                  y = sim.dat[year == yy, .(long, lat)])
  }
  
  ## add in the iid cluste RE
  if(!is.null(clust.re.var)){
    if(verbose) message('-- adding in survey cluster RE')
    cluster_re <- rnorm(n=nrow(sim.dat), mean=0, sd=sqrt(clust.re.var))
  } else{
    cluster_re <- rep(0, nrow(sim.dat))
  }
  
  ## sim sample size of observations
  sim.dat[, N := rpois(n = nrow(sim.dat), lambda = m.clust)]
  
  if(data.lik == 'binom'){
    ## p_true is spatial + cluster!
    ## not the same as what the values from the truth raster
    sim.dat[, p_true := inv.logit(true_lin_pred + cluster_re)] 
           
    ## and now we simulate binomial observations from the true surface
    sim.dat[, Y := rbinom(n = nrow(sim.dat), size = sim.dat[, N], prob = sim.dat[, p_true])]
    
    ## and get empirical p_obs
    sim.dat[, p_obs := Y / N]
    
    ## lastly, we tag on weight. TODO flesh out weighting options
    sim.dat[, weight := 1]
    
  }else if(data.lik == 'normal'){
    sim.dat[, p_true := true_lin_pred + cluster_re] ## identity transform
    
    ## and now we simulate normal observations from the true surface
    ## each individual observation has sd(obs) = sd.norm. 
    ## the sd of the mean of the cluster is sd.norm/sqrt(n)
    sim.dat[, Y := apply(.SD, 1, function(x) mean(rnorm(n = x[1], mean = x[2], sd = sd.norm))), 
            .SDcols = c('N', 'p_true')]
    
    ## and get empirical p_obs
    sim.dat[, p_obs := Y]
  }
  
  ## now we just finish by making some convenience objects and saving everything
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## #########################################################################################
  ## for convenience, we also extract covariate values to the same df (and the true gp val) ##
  ## ##########################################################################################
  
  if(verbose) message('PREPARE AND SAVE OBJECTS\n')
  
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

######
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

#######
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
