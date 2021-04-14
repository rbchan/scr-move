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

## library(scrmove)


## Load the data
load("../deer_scr_telem.RData")


## Set up the workers (one core for each chain)
library(parallel)

nChains <- 4
## (nCores <- detectCores())
## nThreads <- nCores/nChains
nThreads <- 3  ## Fast enough


if(file.exists("parallel_nomove_outfile.Rout"))
    unlink("parallel_nomove_outfile.Rout")

cl3 <- makeCluster(nChains, outfile="parallel_nomove_outfile.Rout")


clusterExport(cl3, c("y.cap", "y.det",
                     "cam.locs", "oper", "nThreads"))


## Set environment variables and load package on each core
invisible(clusterEvalQ(cl3, {
    Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
    Sys.setenv("PKG_LIBS" = "-fopenmp")
    ## Sys.setenv("OMP_PLACES" = "threads")
    ## Sys.setenv(OMP_PROC_BIND = "true")
    ## Sys.setenv(OMP_PROC_BIND = "spread")
    Sys.setenv(OMP_NUM_THREADS = nThreads)
##    print(Sys.getenv())
    .libPaths("./")
    library(scrmove, lib.loc=.libPaths()[1])
}))


## Do a short run to evaluate performance and burnin.
system.time({
    fm3p.1 <- clusterEvalQ(cl3, {
        fmp1 <- scrSemi(ycap=y.cap,
                        ydet=y.det,
                        x=cam.locs,
                        oper=oper,
                        aggregate.ydet=TRUE,
                        plotit=FALSE,
                        n.iters=5000, n.mci=3000,
                        buffer=5000, trim=5000, nthreads=nThreads,
                        s.post.keep=100,
                        report=100, verbose=FALSE,
                        tune=c(0.29, 0.0005, 30, 0.11, 600))
        return(fmp1)
    })
}) 



save(fm3p.1, file="fm3p_1.gzip")



## Do a longer run without storing posterior samples of s and u
system.time({
    fm3p.2 <- clusterEvalQ(cl3, {
        fmp2 <- scrSemi(ycap=y.cap,
                        ydet=y.det,
                        x=cam.locs,
                        oper=oper,
                        aggregate.ydet=TRUE,
                        plotit=FALSE,
                        n.iters=10000, n.mci=3000,
                        buffer=5000, trim=5000, nthreads=nThreads,
                        report=100, verbose=FALSE,
                        inits=fmp1$final.state,
                        tune=c(0.29, 0.0005, 30, 0.11, 600))
        return(fmp2)
    })
}) 


save(fm3p.2, file="fm3p_2.gzip")



## Keep every 20th sample of s
## TODO: Run this with increasing values of n.mci
system.time({
    fm3p.3 <- clusterEvalQ(cl3, {
        fmp3 <- scrSemi(ycap=y.cap,
                        ydet=y.det,
                        x=cam.locs,
                        oper=oper,
                        aggregate.ydet=TRUE,
                        plotit=FALSE,
                        n.iters=20000, n.mci=3000,
                        buffer=5000, trim=5000, nthreads=nThreads,
                        s.post.keep=20,
                        report=100, verbose=FALSE,
                        inits=fmp2$final.state,
                        tune=c(0.29, 0.0005, 30, 0.11, 600))
        return(fmp3)
    })
}) ##

save(fm3p.3, file="fm3p_3.gzip")





stopCluster(cl3)


