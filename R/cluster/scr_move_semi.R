## Fit SCR+move model in JAGS and using custom Gibbs sampler
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


if(file.exists("parallel_semi_outfile.Rout")) unlink("parallel_semi_outfile.Rout")

cl1 <- makeCluster(nChains, outfile="parallel_semi_outfile.Rout")


clusterExport(cl1, c("y.cap", "y.det", "u", "cam.locs", "oper", "nThreads"))


## Set environment variables and load package on each core
invisible(clusterEvalQ(cl1, {
    Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
    Sys.setenv("PKG_LIBS" = "-fopenmp")
    ## Sys.setenv("OMP_PLACES" = "threads")
    ## Sys.setenv(OMP_PROC_BIND = "true")
    ## Sys.setenv(OMP_PROC_BIND = "spread")
    Sys.setenv(OMP_NUM_THREADS = nThreads)
##    print(Sys.getenv())
##    source("../scr_move_semi_capdet_mcmc.R")
    .libPaths("./")
    library(scrmove, lib.loc=.libPaths()[1])
}))


## Do a short run to evaluate performance and burnin.
system.time({
    fm1p.1 <- clusterEvalQ(cl1, {
        fmp1 <- scrMoveSemi(ycap=y.cap,
                            ydet=y.det,
                            u=u,
                            x=cam.locs,
                            oper=oper,
                            plotit=FALSE,
                            n.iters=5000, n.mci=2000,
                            buffer=5000, trim=600, nthreads=nThreads,
                            report=100, verbose=FALSE,
                            block.rhosigma=TRUE,
                            block.lam0kappa=TRUE,
                            ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                            tune=c(0.29, 6.6e-4, 5.5, 0.15, 1.2, 0.11, 300, 300, 3.01e-3, -0.05))
        return(fmp1)
    })
}) 



save(fm1p.1, file="fm1p_1.gzip")



## Do a longer run without storing posterior samples of s and u
system.time({
    fm1p.2 <- clusterEvalQ(cl1, {
        fmp2 <- scrMoveSemi(ycap=y.cap,
                            ydet=y.det,
                            u=u,
                            x=cam.locs,
                            oper=oper,
                            plotit=FALSE,
                            n.iters=10000, n.mci=2000,
                            buffer=5000, trim=600, nthreads=nThreads,
                            report=100, verbose=FALSE,
                            inits=fmp1$final.state,
                            block.rhosigma=TRUE,
                            block.lam0kappa=TRUE,
                            ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                            ## tune=c(0.29, 0.0007, 4, 0.15, 1.1, 0.11, 300, 300, -0.05))
                            tune=c(0.29, 6.6e-4, 5.5, 0.15, 1.2, 0.11, 300, 300, 3.01e-3, -0.05))
        return(fmp2)
    })
}) ## 1990 it/hr


save(fm1p.2, file="fm1p_2.gzip")



## Keep every 20th sample of s and u
## TODO: Run this with increasing values of n.mci
system.time({
    fm1p.3 <- clusterEvalQ(cl1, {
        fmp3 <- scrMoveSemi(ycap=y.cap,
                            ydet=y.det,
                            u=u,
                            x=cam.locs,
                            oper=oper,
                            plotit=FALSE,
                            n.iters=200000, n.mci=2000,
                            buffer=5000, trim=600, nthreads=nThreads,
                            su.post.keep=200,
                            report=100, verbose=FALSE,
                            inits=fmp2$final.state,
                            block.rhosigma=TRUE,
                            block.lam0kappa=TRUE,
                            ## tune=c(0.29, 0.0007, 4, 0.15, 1.1, 0.11, 300, 300, -0.05))
                            tune=c(0.29, 6.6e-4, 5.5, 0.15, 1.2, 0.11, 300, 300, 3.01e-3, -0.05))
        return(fmp3)
    })
}) ##

save(fm1p.3, file="fm1p_3.gzip")





stopCluster(cl1)


