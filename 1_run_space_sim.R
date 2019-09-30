## this script simulates some realistic datasets for comparison between INLA and TMB
## it leverages existing architecture that the LBD team at IHME has already created
## written by AOZ
## last editted Aug 2019

## options(error = recover)

## ## for working on local laptop
## load('/homes/azimmer/scratch/tmb_space_debug.RData')
## load('~/Desktop/tmb_inla_comp/scratch/tmb_space_debug.RData')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

exp.lvid      <- as.numeric(commandArgs()[4]) ## all we need is to grab which (i.e. row of loopvars),
exp.iter      <- as.numeric(commandArgs()[5])  ## which monte carlo iteration,
main.dir.name <- as.character(commandArgs()[6]) ## the run_date folder name

message(sprintf('ON EXPERIMENT LV ID: %04d', exp.lvid))
message(sprintf('-- ON SIM ITER: %04d', exp.iter))

## set the seed for the entire experiment using minutes, 
## the experiemnt loopvar id, and the iteration
## this seed will rarely repeat across experiments and give unique starting values within each experiment
## NOTE: the max seed integer is (2^31 - 1) = 2147483647
## so, we use the last 3 digits for iter
##     we use the 3 prior digits for exp.lvid
##     that leaves 4 digits, with a max value of 2147, first 4 digits
##     we take the minutes and seconds of run_date, modulo 2146 to get these first 4 digits
seed <- as.numeric(sprintf('%i%03d%03d',
                           (as.numeric(paste0(strsplit(main.dir.name, '_')[[1]][5:6], collapse='')) %% 2146),
                           exp.lvid,   ## experiment loopvar id
                           exp.iter))  ## experiment iteration for given lvid
set.seed(seed)
 

## #####################################################
## setup the environment for singularity R            ##
## load in the loopvars and setup params for this run ##
## setup the spatial modeling domain                  ##
## #####################################################
## Set core_repo location and tmb_repo loc
user      <- Sys.info()['user']
tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)
setwd(tmb_repo)

## set the counter how many loops we've gone through
sim.loop.ct <- -1 ## initial setup
source('./2_setup_experiment_settings.R')

## ######################################################################
## setup convergence checks and loop until met. record number of times ##
## ######################################################################
tmb.converge        <- 0 ## set to 1 when tmb converges
tmb.converge.fails  <- 0 ## counts how many times TMB didn't converge
inla.converge       <- 0 ## set to 1 when inla converges
inla.converge.fails <- 0 ## counts how many times INLA didn't converge

## loop until both methods have converged or one method has failed fifth time
sim.loop.ct <- 0 ## ready to start looping
while( (tmb.converge != 1 | inla.converge != 1) & !(tmb.converge.fails >= 15 | inla.converge.fails >= 15) ) {
  
  ## update the counter
  sim.loop.ct <- sim.loop.ct + 1
  
  ## #####################################################
  ## simulate data and setup obj shared by tmb and inla ##
  ## #####################################################
  source('./3_simulate_data_make_shared_obj.R')
  
  ## ######
  ## TMB ##
  ## ######  
  source('./4_setup_run_predict_tmb.R')
  
  ## update convergence args for while loop
  if(tmb.pd.converge){
    tmb.converge <- 1
  }else{
    tmb.converge <- 0
    tmb.converge.fails <- tmb.converge.fails + 1
  }
  
  ## #######
  ## INLA ##
  ## #######
  source('./5_setup_run_predict_inla.R')
  
  ## update convergence args for while loop
  if(inla.mode.converge){
    inla.converge <- 1
  }else{
    inla.converge <- 0
    inla.converge.fails <- inla.converge.fails + 1
  }
  
  ## update the seed so that a new GP is simulated 
  ## 1e6 is added to update the 'seconds' of the 
  ## minute-second timestamp used to create the seed
  ## this ensures that every iteration within this lvid e
  ## experiment will still have a unique seed
  seed <- (seed + 1e6) %% (2^31 - 1) ## modulo max allowed seed
}

## 
if(tmb.converge==1 & inla.converge==1){
  message('-- BOTH INLA AND TMB CONVERGED!')
}else{
  message('-- AT LEAST ONE METHOD DID NOT CONVERGED!')
}

## #############
## VALIDATION ##
## #############
source('./6_run_validation.R')

## update the tracker to compelted (0)
write.table(x=matrix(c(sim.loop.ct, 0), ncol=2), append=TRUE, 
            file = paste0(jobtrack.dir, 
                          sprintf('exp_%04d_iter_%04d.csv', exp.lvid, exp.iter)), sep=',',
            row.names = F, col.names = F)

## write to bottom of outputs and errors files for ease of checking
message('DONE');print('DONE')

# ## save results from all Nsim monte carlo simulations in this run
# write.csv(complete.summary.metrics, sprintf('%s/validation/summary_metrics_complete.csv', out.dir))


#############
## SCRATCH ##
#############


## ######################################
## load in some real data (if desired) ##
## ######################################

## if(use_real_data){

##   reg <- 'sssa'
##   indicator = 'hiv_test'
##   indicator_group = 'hiv'
##   age <- holdout <- test <- 0
##   yearload <- 'annual'
##   withtag <- TRUE
##   datatag <- '_survey'
##   use_share <- 'FALSE'
##   year_list <- 2000:2016

##   pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)


##   ## load in the region shapefile and prep the boundary
##   gaul_list           <- get_gaul_codes(reg)
##   simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
##   subset_shape        <- simple_polygon_list[[1]]
##   simple_polygon      <- simple_polygon_list[[2]]

##   ## Load list of raster inputs (pop and simple)
##   raster_list        <- build_simple_raster_pop(subset_shape)
##   simple_raster      <- raster_list[['simple_raster']]
##   pop_raster         <- raster_list[['pop_raster']]

##   run_date <- make_time_stamp(TRUE)
##   dt <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
##                         simple      = simple_polygon,
##                         agebin      = age,
##                         removeyemen = TRUE,
##                         pathaddin   = pathaddin,
##                         years       = yearload,
##                         withtag     = as.logical(withtag),
##                         datatag     = datatag,
##                         use_share   = as.logical(use_share),
##                         yl          = year_list)

##   ## just to make sure everything goes smoothly, add in single datapoints to missing years
##   missing.yrs <- setdiff(year_list, unique(dt[, year]))

##   if(length(missing.yrs) > 0){
##     for(yy in missing.yrs){
##       new.row <- dt[1, ]
##       new.row[, weight := 0]
##       new.row[, year := yy]
##       dt <- rbind(dt, new.row)
##     }
##   }


##   dt[, Y:= hiv_test]

## }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## !! moved into realistic_sim_utils.R/sim.realistic.data() on 0ct 3, 2018 !!
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## top_pop_urban <- 0.01
## urban_strat <- 0.4
## urban_thresh <- quantile(probs = (1 - top_pop_urban), na.omit(values(pop_raster)))
## u_r_raster <- pop_raster[[1]] ## urban is 1, rural is 0
## u_r_raster[pop_raster[[1]] < urban_thresh] <- 0
## u_r_raster[pop_raster[[1]] >= urban_thresh] <- 1

## pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)
## u_r.pts <- rasterToPoints(u_r_raster, spatial = TRUE)

## ## reproject sp obj
## geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
## pix.pts <- spTransform(pix.pts, CRS(geo.prj))
## u_r.pts <- spTransform(u_r.pts, CRS(geo.prj))
## proj4string(pix.pts)
## proj4string(u_r.pts)

## ## get coords
## pix.pts@data <- data.frame(pix.pts@data, long=coordinates(pix.pts)[,1],
##                            lat=coordinates(pix.pts)[,2])
## pix.pts.numeric <- as.data.frame(pix.pts@data)

## u_r.pts@data <- data.frame(u_r.pts@data, long=coordinates(u_r.pts)[,1],
##                            lat=coordinates(u_r.pts)[,2])
## u_r.pts.numeric <- as.data.frame(u_r.pts@data)

## sim.rows <- sample(x = 1:nrow(pix.pts.numeric), size = n.clust * length(year_list),
##                    replace = TRUE)
## sim.dat <- as.data.table(pix.pts.numeric[, -1])
## sim.dat <- sim.dat[sim.rows, ]


## u.rows <- sample(x = which(u_r.pts.numeric[, 1] == 1), size = round(n.clust * urban_strat),
##                  replace = TRUE)
## r.rows <- sample(x = which(u_r.pts.numeric[, 1] == 0), size = round(n.clust * (1 - urban_strat)),
##                  replace = TRUE)

## sim.dat <- as.data.table(pix.pts.numeric[, -1])
## sim.dat <- sim.dat[c(u.rows, r.rows), ]


## pdf(sprintf('%s/validation/pop_strat.pdf',out.dir), height=10,width=20)
## par(mfrow = c(1, 2))
## raster::plot(pop_raster[[1]])
## points(sim.dat)
## raster::plot(u_r_raster)
## dev.off()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~