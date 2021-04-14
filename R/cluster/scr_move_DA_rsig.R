## Fit SCR+move model in JAGS and using custom Gibbs sampler
## Must include model for capture of telemetered guys
## Otherwise, some of them would have all zero encounter histories, which
## isn't consistent with the (S)CR likelihood


## First we have to build the R package on the cluster
print(olibs <- .libPaths())
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
nThreads <- nCores/nChains


if(file.exists("parallel_outfile_DA_rsig.Rout")) unlink("parallel_outfile_DA_rsig.Rout")

cl1 <- makeCluster(nChains, outfile="parallel_outfile_DA_rsig.Rout")


clusterExport(cl1, c("y.cap", "y.det", "u", "cam.locs", "oper", "nThreads"))


## Set environment variables and load package on each core
invisible(clusterEvalQ(cl1, {
    Sys.setenv(OMP_PLACES = "threads")
    Sys.setenv(OMP_PROC_BIND = "true")
    Sys.setenv(OMP_PROC_BIND = "spread")
    Sys.setenv(OMP_SCHEDULE = "dynamic")
    Sys.setenv(OMP_NUM_THREADS = nThreads)
    ##    print(Sys.getenv())
    print(system("env | grep ^OMP*"))
##    source("../scr_move_semi_capdet_mcmc.R")
    .libPaths("./")
    library(scrmove, lib.loc=.libPaths()[1])
}))


## if(1==2) {
##     source("../../scrmove/R/scr_move_DA_rsig_capdet_mcmc.R")
## }

## Do a short run to evaluate performance and burnin.
system.time({
    fm5p.1 <- clusterEvalQ(cl1, {
        fmp1 <- scrMoveDArsig(ycap=y.cap,
                              ydet=y.det,
                              u=u,
                              M=300,
                              x=cam.locs,
                              oper=oper,
                              random.sig=TRUE,
                              plotit=FALSE, 
                              n.iters=1000, n.mci=200L,
                              buffer=5000, trim=100, nthreads=nThreads,
                              report=100, verbose=FALSE,
                              ##                              block.rhosigma=TRUE,
                              block.lam0kappa=TRUE,
                              ## tune order: rho, log.sigma, log.sigma.mu, log.sigma.sd
                              ##             lam0, kappa, p, s, u, cov(lam0,kappa)
                              tune=c(6.9e-4, 0.4, 0.04, 0.02, 0.15, 1.5, 0.05, 300, 300, 3.01e-3))
        return(fmp1)
    })
}) 


## save(fmp1, file="fmp4_1-1.gzip")

save(fm5p.1, file="fm5p_1.gzip")



## Do a longer run without storing posterior samples of s and u
system.time({
    fm5p.2 <- clusterEvalQ(cl1, {
        fmp2 <- scrMoveDArsig(ycap=y.cap,
                              ydet=y.det,
                              u=u,
                              M=300,
                              x=cam.locs,
                              oper=oper,
                              random.sig=TRUE,
                              plotit=FALSE,
                              n.iters=10000, n.mci=200L,
                              buffer=5000, trim=100, nthreads=nThreads,
                              report=100, verbose=FALSE,
                              inits=fmp1$final.state,
                              ##                          block.rhosigma=TRUE,
                              block.lam0kappa=TRUE,
                              ## tune order: rho, log.sigma, log.sigma.mu, log.sigma.sd
                              ##             lam0, kappa, p, s, u, cov(lam0,kappa)
                              tune=c(6.9e-4, 0.4, 0.04, 0.02, 0.15, 1.5, 0.05, 300, 300, 3.01e-3))
        return(fmp2)
    })
}) ## 1990 it/hr


save(fm5p.2, file="fm5p_2.gzip")



## Keep every 20th sample of s and u
## TODO: Run this with increasing values of n.mci
system.time({
    fm5p.3 <- clusterEvalQ(cl1, {
        fmp3 <- scrMoveDArsig(ycap=y.cap,
                              ydet=y.det,
                              u=u,
                              M=300,
                              x=cam.locs,
                              oper=oper,
                              random.sig=TRUE,
                              plotit=FALSE,
                              n.iters=100000, n.mci=200L,
                              buffer=5000, trim=100, nthreads=nThreads,
                              su.post.keep=20,
                              report=100, verbose=FALSE,
                              inits=fmp2$final.state,
##                              block.rhosigma=TRUE,
                              block.lam0kappa=TRUE,
                              ## tune order: rho, log.sigma, log.sigma.mu, log.sigma.sd
                              ##             lam0, kappa, p, s, u, cov(lam0,kappa)
                              tune=c(6.9e-4, 0.4, 0.04, 0.02, 0.15, 1.5, 0.05, 300, 300, 3.01e-3))
        return(fmp3)
    })
}) ##

save(fm5p.3, file="fm5p_3.gzip")





stopCluster(cl1)


