## Fit SCR+move model in JAGS and using custom Gibbs sampler
## Must include model for capture of telemetered guys
## Otherwise, some of them would have all zero encounter histories, which
## isn't consistent with the (S)CR likelihood


## First we have to build the R package on the cluster
print(.libPaths())
.libPaths("./")                         ## Add local library for new package
print(.libPaths())

install.packages("../../scrmove", repos=NULL, lib=.libPaths()[1])


## Load the data
load("../deer_scr_telem.RData")


nThreads <- 2L


## Sys.setenv("OMP_PLACES" = "threads")
Sys.setenv(OMP_PROC_BIND = "true")
## Sys.setenv(OMP_PROC_BIND = "spread")
## Sys.setenv(OMP_SCHEDULE = "static")
Sys.setenv(OMP_SCHEDULE = "dynamic")
Sys.setenv(OMP_NUM_THREADS = nThreads)
## print(Sys.getenv())
## .libPaths("./")
print(.libPaths())
library(scrmove, lib.loc=.libPaths()[1])


gc()

## Do a short run to evaluate performance and burnin.
system.time({
    fm4.1 <- scrMoveDA(ycap=y.cap,
                      ydet=y.det,
                      u=u,
                      M=300,
                      x=cam.locs,
                      oper=oper,
                      plotit=FALSE, 
                      n.iters=1000, n.mci=100L,
                      buffer=5000, trim=100, nthreads=nThreads,
                      report=10, verbose=FALSE,
                      block.rhosigma=TRUE,
                      block.lam0kappa=TRUE,
                      ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                      tune=c(0.29, 6.9e-4, 5.7, 0.15, 1.2, 0.05, 300, 300, 3.01e-3, -0.05))
}) 



save(fm4.1, file="fm4_1-1.gzip")



## Do a longer run without storing posterior samples of s and u
system.time({
    fm4.2 <- scrMoveDA(ycap=y.cap,
                       ydet=y.det,
                       u=u,
                       M=300,
                       x=cam.locs,
                       oper=oper,
                       plotit=FALSE,
                       n.iters=10000, n.mci=100L,
                       buffer=5000, trim=100, nthreads=nThreads,
                       report=100, verbose=FALSE,
                       inits=fm4.1$final.state,
                       block.rhosigma=TRUE,
                       block.lam0kappa=TRUE,
                       ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                       ## tune=c(0.29, 0.0007, 4, 0.15, 1.1, 0.11, 300, 300, -0.05))
                       tune=c(0.29, 6.6e-4, 5.5, 0.15, 1.2, 0.11, 300, 300, 3.01e-3, -0.05))
}) ## 


save(fm4.2, file="fm4_2-1.gzip")



## Keep every 20th sample of s and u
## TODO: Run this with increasing values of n.mci
system.time({
    fm4.3 <- scrMoveDA(ycap=y.cap,
                       ydet=y.det,
                       u=u,
                       M=300,
                       x=cam.locs,
                       oper=oper,
                       plotit=FALSE,
                       n.iters=100000, n.mci=100L,
                       buffer=5000, trim=100, nthreads=nThreads,
                       su.post.keep=200,
                       report=100, verbose=FALSE,
                       inits=fm4.2$final.state,
                       block.rhosigma=TRUE,
                       block.lam0kappa=TRUE,
                       ## tune=c(0.29, 0.0007, 4, 0.15, 1.1, 0.11, 300, 300, -0.05))
                       tune=c(0.29, 6.6e-4, 5.5, 0.15, 1.2, 0.11, 300, 300, 3.01e-3, -0.05))
}) ##

save(fm4.3, file="fm4_3-1.gzip")




