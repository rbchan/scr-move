## Fit SCR model using custom Gibbs sampler
## Must include model for capture of telemetered guys
## Otherwise, some of them would have all zero encounter histories, which
## isn't consistent with the (S)CR likelihood


## First we have to build the R package on the cluster
print(.libPaths())
.libPaths("./")                         ## Add local library for new package
print(.libPaths())
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp") ## Flags might not be necessary
Sys.setenv("PKG_LIBS" = "-fopenmp")

install.packages("../../scrmove", repos=NULL, lib=.libPaths()[1])

## library(scrmove)


## Load the data
load("../deer_scr_telem.RData")


## Set up the workers (one core for each chain)
library(parallel)

nChains <- 4
(nCores <- detectCores())
nThreads <- nCores/nChains


if(file.exists("parallel_outfile_DA_nomove.Rout")) unlink("parallel_outfile_DA_nomove.Rout")

cl1 <- makeCluster(nChains, outfile="parallel_outfile_DA_nomove.Rout")


clusterExport(cl1, c("y.cap", "y.det", "cam.locs", "oper", "nThreads"))


## Set environment variables and load package on each core
invisible(clusterEvalQ(cl1, {
    Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
    Sys.setenv("PKG_LIBS" = "-fopenmp")
    ## Sys.setenv("OMP_PLACES" = "threads")
    Sys.setenv(OMP_PROC_BIND = "true")
    ## Sys.setenv(OMP_PROC_BIND = "spread")
    Sys.setenv(OMP_SCHEDULE = "static")
    Sys.setenv(OMP_NUM_THREADS = nThreads)
    ##    print(Sys.getenv())
    print(system("env | grep ^OMP*"))
##    source("../scr_move_semi_capdet_mcmc.R")
    .libPaths("./")
    library(scrmove, lib.loc=.libPaths()[1])
}))

## plot(rmvnorm(1000, c(0.001, 300), matrix(c(0.0002^2, -.001, -.001, 20^2),2)))

## Do a short run to evaluate performance and burnin.
system.time({
    fm6p.1 <- clusterEvalQ(cl1, {
        fmp1 <- scrDA(ycap=y.cap,
                      ydet=y.det,
                      M=300,
                      x=cam.locs,
                      oper=oper,
                      plotit=FALSE, 
                      n.iters=1000, 
                      buffer=5000, trim=3000,
                      nthreads=nThreads,
                      report=10, verbose=FALSE,
                      ## block.lam0sigma=TRUE,
                      ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                      tune=c(0.0002, 25, 500, -0.001))
        return(fmp1)
    })
}) ## 710 it/hr

## xx <- as.mcmc.list(lapply(fm6p.3, function(x) as.mcmc(x$samples)))



## save(fmp1, file="fmp4_1-1.gzip")

save(fm6p.1, file="fm6p_1.gzip")



## Do a longer run without storing posterior samples of s and u
system.time({
    fm6p.2 <- clusterEvalQ(cl1, {
        fmp2 <- scrDA(ycap=y.cap,
                      ydet=y.det,
                      M=300,
                      x=cam.locs,
                      oper=oper,
                      plotit=FALSE,
                      n.iters=10000, 
                      buffer=5000, trim=3000,
                      nthreads=nThreads,
                      report=100, verbose=FALSE,
                      inits=fmp1$final.state,
                      ## block.lam0sigma=TRUE,
                      tune=c(0.0002, 25, 500, -0.001))
        return(fmp2)
    })
}) ## 1990 it/hr


save(fm6p.2, file="fm6p_2.gzip")



## Keep every 20th sample of s and u
## TODO: Run this with increasing values of n.mci
system.time({
    fm6p.3 <- clusterEvalQ(cl1, {
        fmp3 <- scrDA(ycap=y.cap,
                      ydet=y.det,
                      M=300,
                      x=cam.locs,
                      oper=oper,
                      plotit=FALSE,
                      n.iters=100000, 
                      buffer=5000, trim=3000,
                      nthreads=nThreads,
                      s.post.keep=20,
                      report=100, verbose=FALSE,
                      inits=fmp2$final.state,
                      ## block.lam0sigma=TRUE,
                      ## tune=c(0.29, 0.0007, 4, 0.15, 1.1, 0.11, 300, 300, -0.05))
                      tune=c(0.0002, 25, 500, -0.001))
        return(fmp3)
    })
}) ##

save(fm6p.3, file="fm6p_3.gzip")





stopCluster(cl1)


