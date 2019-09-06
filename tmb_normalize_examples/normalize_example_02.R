################################################################################
## Testing priors with normalize()
################################################################################

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

## Sample nu ~ N(0, sigma^2)
sigma <- exp(-1)
set.seed(4343)
nu <- rnorm(nrow(A), 0, sigma)

## Sample data
set.seed(1002)
alpha <- -3
expit <- function(x) { exp(x) /(1+exp(x))}
prob <- as.vector(expit(alpha + A %*% S + nu))
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
                   logsigma = -0.9,
                   S = rep(0, nrow(Q)),
                   nu = rep(0, length(y)))

Random <- c("S", "nu")

## Uses normalize flag
compile("normalize_example.cpp") 
dyn.load(dynlib("normalize_example"))

## Does not use normalize flag (GMRF with normalizing constant already added in)
compile("nonnormalize_example.cpp") 
dyn.load(dynlib("nonnormalize_example"))

##### EXAMPLE 1
# alpha ~ N(-3, 1)
# logtau ~ N(-1.6, 1)
# logkappa ~ N(0.4, 1)
# logsigma ~ N(-1, 1)
Data$alphamean <- -3
Data$logtaumean <- -1.6
Data$logkappamean <- 0.4
Data$logsigmamean <- -1

# normalize
obj1a <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="normalize_example")
obj1a <- normalize(obj1a, flag = "flag")
opt1a <- do.call("optim", list(par=obj1a$par, fn=obj1a$fn, gr=obj1a$gr, method="BFGS"))
sd1a <- sdreport(obj1a)

# nonnormalize
obj1b <- MakeADFun(Data, Parameters, random = Random, hessian = T, DLL="nonnormalize_example")
opt1b <- do.call("optim", list(par=obj1b$par, fn=obj1b$fn, gr=obj1b$gr, method="BFGS"))
sd1b <- sdreport(obj1b)

#### Compare
sd1a

sd1b
