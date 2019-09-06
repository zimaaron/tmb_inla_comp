################################################################################
## Testing priors with normalize()
################################################################################

code.dir <- '/Users/azimmer/Desktop/tmb_normalize_examples/'

### Load Libraries
library(TMB)
library(INLA)

## Set up example data (use locations from PRprec)
data(PRprec)
coords <- as.matrix(PRprec[,1:2])
pts.bound <- inla.nonconvex.hull(coords, 0.3, 0.3)

## Set up mesh, projection from mesh --> data locations, spde
mesh <- inla.mesh.2d(coords, boundary=pts.bound,
                     max.edge=c(0.3,1), offset=c(1e-5,1.5), cutoff=0.1)
A <- inla.spde.make.A(mesh, loc=coords)
spde <- inla.spde2.matern(mesh=mesh, alpha=2)

## Sample S ~ N(0, Q^{-1})
Q <- inla.spde2.precision(spde, theta= c(-1.5, 0.4))
S <- inla.qsample(n=1L, Q=Q, seed = 111L)

## Sample data
set.seed(1002)
alpha <- -3
expit <- function(x) { exp(x) /(1+exp(x))}
prob <- as.vector(expit(alpha + A %*% S))
n <- rep(50, length(prob))
y <- rbinom(length(n), n, prob)

## TMB inputs
Data <- list( flag = 1,
              y = y,
              n = n,
              M0 = spde$param.inla$M0,
              M1 = spde$param.inla$M1,
              M2 = spde$param.inla$M2,
              Aproj = A)

Parameters <- list(alpha = -3.5,
                   logtau = spde$param.inla$theta.initial[1],
                   logkappa = spde$param.inla$theta.initial[2],
                   S = rep(0, nrow(Q)))

Random <- c("S")

## Uses normalize flag
compile(paste0(code.dir, "normalize_example_01.cpp")) 
dyn.load(dynlib(paste0(code.dir, "normalize_example_01")))

## Does not use normalize flag (GMRF with normalizing constant already added in)
compile(paste0(code.dir, "nonnormalize_example_01.cpp")) 
dyn.load(dynlib(paste0(code.dir, "nonnormalize_example")))

##### EXAMPLE 1
# alpha ~ N(-3, 1)
# logtau ~ N(-1.6, 1)
# logkappa ~ N(0.4, 1)
Data$alphamean <- -3
Data$logtaumean <- -1.6
Data$logkappamean <- 0.4

# normalize
obj1a <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="normalize_example")
obj1a <- normalize(obj1a, flag = "flag")
opt1a <- do.call("optim", list(par=obj1a$par, fn=obj1a$fn, gr=obj1a$gr, method="BFGS"))
sd1a <- sdreport(obj1a)

# nonnormalize
obj1b <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="nonnormalize_example")
opt1b <- do.call("optim", list(par=obj1b$par, fn=obj1b$fn, gr=obj1b$gr, method="BFGS"))
sd1b <- sdreport(obj1b)


##### EXAMPLE 2: change alphamean
# alpha ~ N(0, 1)
# logtau ~ N(-1.6, 1)
# logkappa ~ N(0.4, 1)
Data$alphamean <- 0
Data$logtaumean <- -1.6
Data$logkappamean <- 0.4

# normalize
obj2a <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="normalize_example")
obj2a <- normalize(obj2a, flag = "flag")
opt2a <- do.call("optim", list(par=obj2a$par, fn=obj2a$fn, gr=obj2a$gr, method="BFGS"))
sd2a <- sdreport(obj2a)

# nonnormalize
obj2b <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="nonnormalize_example")
opt2b <- do.call("optim", list(par=obj2b$par, fn=obj2b$fn, gr=obj2b$gr, method="BFGS"))
sd2b <- sdreport(obj2b)



##### EXAMPLE 3: change logtaumean & logkappamean
# alpha ~ N(-3, 1)
# logtau ~ N(-5, 1)
# logkappa ~ N(3, 1)
Data$alphamean <- -3
Data$logtaumean <- -3
Data$logkappamean <- 3

# normalize
obj3a <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="normalize_example")
obj3a <- normalize(obj3a, flag = "flag")
opt3a <- do.call("optim", list(par=obj3a$par, fn=obj3a$fn, gr=obj3a$gr, method="BFGS"))
sd3a <- sdreport(obj3a)

# nonnormalize
obj3b <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="nonnormalize_example")
opt3b <- do.call("optim", list(par=obj3b$par, fn=obj3b$fn, gr=obj3b$gr, method="BFGS"))
sd3b <- sdreport(obj3b)


#### Compare
sd1a
sd2a
sd3a

sd1b
sd2b
sd3b
