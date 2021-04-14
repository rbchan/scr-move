## Fit SCR+move model in JAGS and using custom Gibbs sampler
## Must include model for capture of telemetered guys
## Otherwise, some of them would have all zero encounter histories, which
## isn't consistent with the (S)CR likelihood
## This script ignores the telemetry data


## First we have to build the R package on the cluster
print(.libPaths())
.libPaths("./")                         ## Add local library for new package
print(.libPaths())
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp") ## Flags might not be necessary
Sys.setenv("PKG_LIBS" = "-fopenmp")

install.packages("../../scrmove", repos=NULL, lib=.libPaths()[1])

library(scrmove)


## Load the data
load("../deer_scr_telem.RData")


## Set up the workers (one core for each chain)
library(parallel)

nChains <- 4
(nCores <- detectCores())
(nThreads <- min(nCores/nChains, 4))


if(file.exists("parallel_nou_outfile.Rout"))
    unlink("parallel_nou_outfile.Rout")

cl2 <- makeCluster(nChains, outfile="parallel_nou_outfile.Rout")


clusterExport(cl2, c("y.cap", "y.det",
                     ## "u", ## ignore telemetry data
                     "cam.locs", "oper", "nThreads"))


## Set environment variables and load package on each core
invisible(clusterEvalQ(cl2, {
    Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
    Sys.setenv("PKG_LIBS" = "-fopenmp")
    Sys.setenv("OMP_PLACES" = "threads")
    ## Sys.setenv(OMP_PROC_BIND = "true")
    Sys.setenv(OMP_PROC_BIND = "spread")
    Sys.setenv(OMP_NUM_THREADS = nThreads)
##    print(Sys.getenv())
##    source("../scr_move_semi_capdet_mcmc.R")
    .libPaths("./")
    library(scrmove, lib.loc=.libPaths()[1])
}))


## Do a short run to evaluate performance and burnin.
system.time({
    fm2p.1 <- clusterEvalQ(cl2, {
        fmp1 <- scrMoveSemi(ycap=y.cap,
                            ydet=y.det,
                            ## u=u,  ## ignore telemetry data
                            x=cam.locs,
                            oper=oper,
                            plotit=FALSE,
                            n.iters=5000, n.mci=2000,
                            buffer=5000, trim=600, nthreads=nThreads,
                            report=100, verbose=FALSE,
##                            block.rhosigma=TRUE,
                            block.lam0kappa=TRUE,
                            pu.prob=1,
                            ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                            tune=c(0.29, 0.002, 5, 0.15, 1.1, 0.10, 300, 300, -0.02, -0.05))
        return(fmp1)
    })
}) 


## save(fm2p.1, file="fm2p_1_desktop.gzip")

save(fm2p.1, file="fm2p_1.gzip")



## Do a longer run without storing posterior samples of s and u
system.time({
    fm2p.2 <- clusterEvalQ(cl2, {
        fmp2 <- scrMoveSemi(ycap=y.cap,
                            ydet=y.det,
##                            u=u, ## ignore telemetry data
                            x=cam.locs,
                            oper=oper,
                            plotit=FALSE,
                            n.iters=10000, n.mci=2000,
                            buffer=5000, trim=600, nthreads=nThreads,
                            report=100, verbose=FALSE,
##                            block.rhosigma=TRUE,
                            block.lam0kappa=TRUE,
                            inits=fmp1$final.state,
                            pu.prob=1,
                            tune=c(0.29, 0.002, 5, 0.15, 1.1, 0.10, 300, 300, -0.02, -0.05))
        return(fmp2)
    })
}) ## 

## save(fm2p.2, file="fm2p_2_desktop.gzip")

save(fm2p.2, file="fm2p_2.gzip")



## Keep every 20th sample of s and u
## TODO: Run this with increasing values of n.mci
system.time({
    fm2p.3 <- clusterEvalQ(cl2, {
        fmp3 <- scrMoveSemi(ycap=y.cap,
                            ydet=y.det,
##                            u=u, ## ignore telemetry data
                            x=cam.locs,
                            oper=oper,
                            plotit=FALSE,
                            n.iters=200000, n.mci=2000,
                            buffer=5000, trim=600, nthreads=nThreads,
                            su.post.keep=200,
                            report=1000, verbose=FALSE,
##                            block.rhosigma=TRUE,
                            block.lam0kappa=TRUE,
                            inits=fmp2$final.state,
                            pu.prob=1,
                            tune=c(0.29, 0.002, 5, 0.15, 1.1, 0.10, 300, 300, -0.02, -0.05))
        return(fmp3)
    })
}) ##

## save(fm2p.3, file="fm2p_3_desktop.gzip")

save(fm2p.3, file="fm2p_3.gzip")





stopCluster(cl2)




## rat = p(y|u.cand)p(u.cand|s.cand)p(s.cand)q(su) / 
##       p(y|u)p(u|s)p(s)q(su.cand)

## rat = p(y|u.cand) / 
##       p(y|u)
