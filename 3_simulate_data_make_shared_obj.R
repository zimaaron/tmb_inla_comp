## ################
## simulate data ##
## ################

message('---- ON SCRIPT 3: simulating data and prepping objects for fitting')

## update the tracker
write.table(x=matrix(c(sim.loop.ct, 3), ncol=2), append=T,
            file = paste0(jobtrack.dir, 
                          sprintf('exp_%04d_iter_%04d.csv', exp.lvid, exp.iter)), sep=',', 
            row.names=F, col.names = F)

if(exp.iter == 1){ ## first time, must load covs, after that, we can reuse them
  sim.obj <- sim.realistic.data(reg = reg,
                                year_list = year_list,
                                data.lik = data.lik,
                                sd.norm = sqrt(norm.var), 
                                alpha = alpha, 
                                betas = betas,
                                sp.kappa = sp.kappa,
                                sp.alpha = sp.alpha,
                                t.rho = t.rho,
                                pixel.iid.var = NULL,
                                n.clust = n.clust,
                                m.clust = m.clust,
                                clust.re.var = clust.var, 
                                covs = covs,
                                cov_layers = NULL, 
                                simple_raster = simple_raster,
                                simple_polygon = simple_polygon,
                                out.dir = out.dir,
                                pop_raster = pop_raster, 
                                obs.loc.strat = obs.loc.strat,
                                urban.pop.pct = urban.pop.pct,
                                urban.strat.pct = urban.strat.pct, 
                                sp.field.sim.strat = 'SPDE', 
                                seed = seed,
                                verbose = FALSE,
                                exp.iter = exp.iter)
  
  ## save the cov_list for future iterations to speed things up
  covs.gp <- sim.obj$cov.gp.rasters   ## rasters of covs and true simulated gp field
  cov_list <- covs.gp[!grepl('gp',names(covs.gp))]
  ## save the cov_list to reload in future iterations of this experiment
  saveRDS(object = cov_list,
          file = sprintf('%s/cov_list.rds', common.dir))

}else{
  
  ## reuse covs
  cov_list <- readRDS(sprintf('%s/cov_list.rds', common.dir))
  
  sim.obj <- sim.realistic.data(reg = reg,
                                year_list = year_list,
                                data.lik = data.lik,
                                sd.norm = sqrt(norm.var),
                                alpha = alpha, 
                                betas = betas,
                                sp.kappa = sp.kappa,
                                sp.alpha = sp.alpha,
                                t.rho = t.rho,
                                clust.re.var = clust.var, 
                                n.clust = n.clust,
                                m.clust = m.clust,
                                covs = covs,
                                cov_layers = cov_list, ## which is created from each sim.obj after stripping GP from the list. ~line64
                                simple_raster = simple_raster,
                                simple_polygon = simple_polygon,
                                pop_raster = pop_raster, 
                                obs.loc.strat = obs.loc.strat,
                                urban.pop.pct = urban.pop.pct,
                                urban.strat.pct = urban.strat.pct, 
                                out.dir = out.dir,
                                sp.field.sim.strat = 'SPDE', 
                                seed = seed,
                                exp.iter = exp.iter)

saveRDS(file = sprintf('%s/simulated_obj/experiment%04d_iter%04d_sim_obj.rds', 
                       out.dir, exp.lvid, exp.iter),
        object = sim.obj)

}

## process parts of the returned sim obj list into pieces we need for model fitting
dt <- sim.obj$sim.dat ## simulated data, lat-long, year, covs, true surface
covs.gp <- sim.obj$cov.gp.rasters   ## rasters of covs and true simulated gp field
true.gp <- covs.gp[['gp']]
cov_list <- covs.gp[!grepl('gp',names(covs.gp))]
true.rast <- sim.obj$true.rast
if(data.lik == 'binom'){
  true.rast.p <- true.rast
  values(true.rast.p) <- plogis(values(true.rast))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ##################################
## ##################################
## setup for tmb and INLA modeling ##
## ##################################
## ##################################

## ############################
## SETUP SOME SHARED OBJECTS ##
## ############################

## ~~~~~~~~~~~~~~~~~~~~~~
## required for model fit
## ~~~~~~~~~~~~~~~~~~~~~~

## setup space time data locs
dt[, id := 1:.N]
dt[, period_id := as.numeric(as.factor(dt[, year]))]
dt.coords <- as.matrix(dt[, .(long, lat)])
dt.pers   <- dt[, period_id]
nperiods  <- length(year_list)

## get triangulation params 
mesh.params    <- eval(parse(text=mesh_s_params))

## make delaunay triangulation from all pixel centroids in domain

## first, get the centroids
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

## make the triangulation
mesh_s <- inla.mesh.2d(loc.domain = as.matrix(pix.pts.numeric[,2:3]), 
                       max.e = mesh.params)

## plot the triangulation
pdf(sprintf('%s/modeling/inputs/experiment%04d_iter%04d_mesh.pdf', out.dir, exp.lvid, exp.iter))
plot(mesh_s)
plot(simple_raster, add = TRUE) ## just to show loc of simple_raster under mesh for scale
plot(mesh_s, add = TRUE)
points(dt.coords, col = 'red', pch = '.')
dev.off()

nodes <- mesh_s$n ## get number of mesh nodes
## spde <- inla.spde2.matern(mesh_s, alphac = 2)
## TODO - if we move to pc priors, need to adjust draws, plots, and validation
##        bc inla produce range and variance (not theta1,2) by default when using PC
# spde <- inla.spde2.pcmatern(mesh = mesh_s, 
#                             prior.range = c(0.05, 0.01), ## P(range < 0.05) = 0.01
#                             prior.sigma = c(2, 0.01))    ## P(sigma > 2) = 0.01
## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
## ^ this gives us a linear reduction of \Sigma^{-1} as:
## \Sigma^{-1} = \kappa^4 M_0 + 2\kappa^2M_1 + M_2
## M_2 = M_1M_0^{-1}M_1
## Where the Ms are all sparse matrices stored as "dgTMatrix"
## names(spde$param.inla)
spde <- inla.spde2.pcmatern(mesh=mesh_s, alpha =2,
                            ## constr = TRUE, integrate-to-zero constraint
                            prior.range = matern.pri[1:2],
                            prior.sigma = matern.pri[3:4])

## use inla helper functions to project the spatial effect from mesh points to data points
A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = dt.coords,
                           group = dt.pers)

## save relevant objects. NOTE: meshes are not deterministic!?
saveRDS(file = sprintf('%s/modeling/inputs/experiment%04d_iter%04d_mesh.rds', out.dir, exp.lvid, exp.iter), mesh_s)
saveRDS(file = sprintf('%s/modeling/inputs/experiment%04d_iter%04d_spde.rds', out.dir, exp.lvid, exp.iter), spde)

## now that the mesh is made, we can grab the default priors that it generates
## only needed for non-pc matern priors
mesh.info <- param2.matern.orig(mesh_s)
## theta1 = log(tau). The prior on theta1 is:   log(tau)   ~ Normal(mean, var)
spde.theta1.pri <- c(mesh.info$theta.prior.mean[1], mesh.info$theta.prior.prec[1, 1])
## theta2 = log(kappa). The prior on theta2 is: log(kappa) ~ Normal(mean, var)
spde.theta2.pri <- c(mesh.info$theta.prior.mean[2], mesh.info$theta.prior.prec[2, 2])

## ~~~~~~~~~~~~~~~~~~~~
## required for predict
## ~~~~~~~~~~~~~~~~~~~~

## get space-locs grid to predict onto
pcoords <- xyFromCell(simple_raster, which(!is.na(values(simple_raster))))

## replicate in time if needed
if(nperiods > 1){
  pcoords <- do.call(rbind,
                     replicate(nperiod,
                               pcoords,
                               simplify = FALSE))
}

## get time groupings
groups_periods <- rep(1:nperiods, each = nrow(pcoords))

## use inla helper functions to project the spatial effect.
A.pred <- inla.spde.make.A(
  mesh = mesh_s,
  loc = pcoords,
  group = groups_periods)

## extract cell values from covariates, deal with timevarying covariates here

## pull out covariates in format we expect them
## a list of length periods with a brick of named covariates inside
new_cl <- list()
for(p in 1:nperiods){
  new_cl[[p]] <- list()
  for(n in names(cov_list)){
    if(dim(cov_list[[n]])[3] == 1){
      new_cl[[p]][[n]] <- cov_list[[n]]
    }else{
      new_cl[[p]][[n]] <- cov_list[[n]][[p]]
    }
  }
  new_cl[[p]] <- brick(new_cl[[p]])
}

cov_vals <- list()
for(p in 1:nperiods){
  cov_vals[[p]] <- raster::extract(new_cl[[p]], pcoords[1:(sum(!is.na(values(simple_raster)))/nperiods),])
  colnames(cov_vals[[p]]) <- cov_names
  cov_vals[[p]] <- (cbind(int = 1, cov_vals[[p]]))
}