## this script simulates some realistic datasets for comparison between INLA and TMB
## it leverages existing architecture that the LBD team at IHME has already created
## written by AOZ
## last editted Oct 3, 2018

## options(error = recover)

## ## for working on local laptop
## load('/homes/azimmer/scratch/tmb_space_debug.RData')
## load('~/Desktop/tmb_inla_comp/scratch/tmb_space_debug.RData')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par.iter      <- as.numeric(  commandArgs()[4]) ## all we need is to grab the (parallel) iteration of this run
main.dir.name <- as.character(commandArgs()[5]) ## and the main directory for all experiments in this run so we know where to load from

## #####################################################
## setup the environment for singularity R            ##
## load in the loopvars and setup params for this run ##
## setup the spatial modeling domain                  ##
## #####################################################
source('./2_setup_experiment_settings.R')

## ####################################
## perform the experiment Nsim times ##
## ####################################

for(iii in 1:Nsim){ ## repeat Nsim times

  ## ######################################################################
  ## setup convergence checks and loop until met. record number of times ##
  ## ######################################################################
  tmb.converge <- 0 ## set to 1 when tmb converges
  tmb.converge.fails <- 0 ## counts how many times TMB didn't converge
  inla.converge <- 0 ## set to 1 when inla converges
  inla.conerge.fails <- 0 ## counts how many times INLA didn't converge

  ## loop until both methods have converged or one method has failed fifth time
  while(tmb.converge != 1 | inla.converge != 1 | tmb.converge.fails < 5 | inla.converge.fails < 5) {
    
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
    
  }

  ## #############
  ## VALIDATION ##
  ## #############
  source('./6_run_validation.R')


} ## end iii loop repeating iterations over 1:Nsim

## save results from all Nsim monte carlo simulations in this run
write.csv(complete.summary.metrics, sprintf('%s/validation/summary_metrics_complete.csv',out.dir))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## (i) compare data to fitted surface
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




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

