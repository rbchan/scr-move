## Fit SCR+move model in JAGS and using custom Gibbs sampler
## Must include model for capture of telemetered guys
## Otherwise, some of them would have all zero encounter histories, which
## isn't consistent with the (S)CR likelihood


## First we have to build the R package on the cluster
print(olp <- .libPaths())
.libPaths("./")                         ## Add local library for new package
print(.libPaths())

install.packages("../../scrmove", repos=NULL, lib=.libPaths()[1])

## library(scrmove)


## Load the data
load("../deer_scr_telem.RData")


## Set up the workers (one core for each chain)
library(parallel)

nChains <- 4
(nCores <- detectCores())
nThreads <- as.integer(nCores/nChains)


if(file.exists("parallel_outfile_DA.Rout")) unlink("parallel_outfile_DA.Rout")

cl1 <- makeCluster(nChains, outfile="parallel_outfile_DA.Rout")


clusterExport(cl1, c("y.cap", "y.det", "u", "cam.locs", "oper", "nThreads"))


## Set environment variables and load package on each core
invisible(clusterEvalQ(cl1, {
    ## Sys.setenv("OMP_PLACES" = "threads")
    Sys.setenv(OMP_PROC_BIND = "true")
    Sys.setenv(OMP_PROC_BIND = "spread")
    Sys.setenv(OMP_SCHEDULE = "dynamic")
    Sys.setenv(OMP_NUM_THREADS = nThreads)
##    print(Sys.getenv())
    print(system("env | grep ^OMP_*"))
    .libPaths("./")
    library(scrmove, lib.loc=.libPaths()[1])
}))


## Do a short run to evaluate performance and burnin.
system.time({
    fm4p.1 <- clusterEvalQ(cl1, {
        fmp1 <- scrMoveDA(ycap=y.cap,
                          ydet=y.det,
                          u=u,
                          M=300,
                          x=cam.locs,
                          oper=oper,
                          plotit=FALSE, 
                          n.iters=1000, n.mci=100L,
                          buffer=5000, trim=100, nthreads=nThreads,
                          report=100, verbose=FALSE,
                          block.rhosigma=TRUE,
                          block.lam0kappa=TRUE,
                          ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                          tune=c(0.29, 6.9e-4, 5.7, 0.15, 1.2, 0.05, 300, 300, 3.01e-3, -0.05))
        return(fmp1)
    })
}) 



save(fm4p.1, file="fm4p_1.gzip")



## Do a longer run without storing posterior samples of s and u
system.time({
    fm4p.2 <- clusterEvalQ(cl1, {
        fmp2 <- scrMoveDA(ycap=y.cap,
                          ydet=y.det,
                          u=u,
                          M=300,
                          x=cam.locs,
                          oper=oper,
                          plotit=FALSE,
                          n.iters=10000, n.mci=100L,
                          buffer=5000, trim=100, nthreads=nThreads,
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


save(fm4p.2, file="fm4p_2.gzip")



## Keep every 20th sample of s and u
## TODO: Run this with increasing values of n.mci
system.time({
    fm4p.3 <- clusterEvalQ(cl1, {
        fmp3 <- scrMoveDA(ycap=y.cap,
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
                          inits=fmp2$final.state,
                          block.rhosigma=TRUE,
                          block.lam0kappa=TRUE,
                          ## tune=c(0.29, 0.0007, 4, 0.15, 1.1, 0.11, 300, 300, -0.05))
                          tune=c(0.29, 6.6e-4, 5.5, 0.15, 1.2, 0.11, 300, 300, 3.01e-3, -0.05))
        return(fmp3)
    })
}) ##

save(fm4p.3, file="fm4p_3.gzip")





stopCluster(cl1)


