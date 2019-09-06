## plots made for TMB-INLA paper
## AOZ - SEP 2018

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## varying meshes used in study 3
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## set triangulation params: c(cutoff, max.edge[1], max.edge[2])
## default: c(0.1 1.0 5.0)

## in order to add title...
custom.plot.inla.mesh <- function (x, col = "white", t.sub = 1:nrow(mesh$graph$tv), add = FALSE, 
                            lwd = 1, xlim = range(mesh$loc[, 1]), 
                            ylim = range(mesh$loc[, 2]), main = NULL, rgl = FALSE, size = 2, 
                            draw.vertices = FALSE, 
                            vertex.color = "black", draw.edges = TRUE, 
                            edge.color = rgb(0.3, 0.3, 0.3), draw.segments = draw.edges, ...) 
{
  INLA:::inla.require.inherits(x, "inla.mesh", "'mesh'")
  mesh = x
  if (rgl) {
    stopifnot(inla.require("rgl"))
    if (!add) {
      dev = rgl::open3d()
      rgl::view3d(0, 0, fov = 0)
    }
    else {
      dev = NULL
    }
    tv = mesh$graph$tv[t.sub, , drop = FALSE]
    if (draw.vertices) {
      idx = intersect(unique(as.vector(tv)), mesh$idx$loc)
      rgl::points3d(mesh$loc[idx, , drop = FALSE], size = 2 * 
                      size, lwd = lwd, color = "blue", ...)
    }
    if (draw.segments) {
      if (!is.null(mesh$segm$bnd)) 
        lines(mesh$segm$bnd, mesh$loc, lwd = lwd + 1, 
              rgl = TRUE, add = TRUE, ...)
      if (!is.null(mesh$segm$int)) 
        lines(mesh$segm$int, mesh$loc, lwd = lwd + 1, 
              rgl = TRUE, add = TRUE, ...)
    }
    plot.inla.trimesh(tv, mesh$loc, color = col, size = size, 
                      lwd = lwd, draw.vertices = draw.vertices, vertex.color = vertex.color, 
                      draw.edges = draw.edges, edge.color = edge.color, 
                      ...)
    return(invisible(dev))
  }
  else {
    idx = cbind(mesh$graph$tv[t.sub, c(1:3, 1), drop = FALSE], 
                NA)
    x = mesh$loc[t(idx), 1]
    y = mesh$loc[t(idx), 2]
    if (!add) {
      plot.new()
      plot.window(xlim = xlim, ylim = ylim, ...)
    }
    if (draw.edges) {
      lines(x, y, type = "l", col = edge.color, lwd = lwd)
    }
    tv = mesh$graph$tv[t.sub, , drop = FALSE]
    if (draw.vertices) {
      idx = unique(as.vector(tv))
      points(mesh$loc[idx, , drop = FALSE], pch = 20, col = vertex.color, 
             ...)
      idx = intersect(idx, mesh$idx$loc)
      points(mesh$loc[idx, , drop = FALSE], pch = 20, col = "blue", 
             ...)
    }
    if (draw.segments) {
      if (!is.null(mesh$segm$bnd)) 
        lines(mesh$segm$bnd, mesh$loc, lwd = lwd + 1, 
              ...)
      if (!is.null(mesh$segm$int)) 
        lines(mesh$segm$int, mesh$loc, lwd = lwd + 1, 
              ...)
    }
    if (!add && missing(main)) {
      if (mesh$meta$is.refined) {
        title("Constrained refined Delaunay triangulation", 
              ...)
      }
      else {
        title("Constrained Delaunay triangulation", ...)
      }
    }else {
      title(main, ...)
    }
    return(invisible())
  }
}

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


mesh_list <- list(c(0.15,  5.0),
                  c(0.20,  5.0),
                  c(0.40,  5.0),
                  c(0.80,  5.0)
                  )
png('~/tmb_inla_paper/mesh_2x2_varying_density3.png', width=8, height=8, units='in', res=350)
par(mfrow=c(2,2), mai=c(0.1, 0.2, 0.6, 0.2))
for(mm in 1:length(mesh_list)){
  
  ## an alternate way to get decent mesh - but the number of points 
  ## on the boundary can be too dense/odd
  # bound <- inla.nonconvex.hull(as.matrix(pix.pts.numeric[,2:3]),
  #                               -0.05,-0.05, resolution=c(100,100))
  # plot(bound$loc);plot(simple_polygon, add=T)
  # mesh_s <- inla.mesh.2d(boundary=bound, max.e=mesh_list[[mm]])

  ## make delaunay triangulation
  mesh_s <- inla.mesh.2d(loc.domain = as.matrix(pix.pts.numeric[,2:3]), 
                         max.e=mesh_list[[mm]])
  
  ## plot mesh
  if(TRUE){## use custom plot to get title correctly added...
    custom.plot.inla.mesh(mesh_s, 
                          main = sprintf('Delaunay FEM Triangulation \n interior edge size ~ %.02f degrees\n number of nodes: %i',
                                         mesh_list[[mm]][1], mesh_s$n))
    
    plot(simple_raster, legend=FALSE, add=TRUE,
         xlim=range(mesh_s$loc[,1]),
         ylim=range(mesh_s$loc[,2]))
    plot(mesh_s, add = TRUE)
  }
  
}
dev.off()

