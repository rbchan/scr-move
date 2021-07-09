## Install and load package

list.files("../../")

install.packages("../../scrmove", repos=NULL, type="source")

library(scrmove)
library(coda)


##  Load data
list.files("../")

load("../deer_scr_telem.RData")


ls()



## Simulator

sim.scr.move <- function(x, buffer, K, ##M,
                         N, sigma, rho,
                         p.cap, lam0, sigma.det,
                         thin.telem, verbose=FALSE) {

    if(!exists(".Random.seed"))
        runif(1)

    iseed <- .Random.seed

    J <- nrow(x) ## nTraps
    xlim <- range(x[,1])+c(-buffer, buffer)
    ylim <- range(x[,2])+c(-buffer, buffer)

    ## Latent activity centers
    s <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]))
    u <- array(NA_real_, c(N, K, 2))
    lambda <- dist2.ux <- array(NA_real_, c(N, J, K))
    y.det.N <- array(NA_integer_, c(N, J, K))
    sigma.step <- sqrt(sigma^2-rho^2*sigma^2)

    ## Latent movement paths
    ## Initial location drawn from stationary dist
    u[,1,] <- cbind(rnorm(N, s[,1], sigma),
                    rnorm(N, s[,2], sigma))
    for(k in 2:K) {
        mu.u1 <- s[,1]+rho*(u[,k-1,1]-s[,1])
        mu.u2 <- s[,2]+rho*(u[,k-1,2]-s[,2])
        u[,k,] <- cbind(rnorm(N, mu.u1, sigma.step),
                        rnorm(N, mu.u2, sigma.step))
    }
    ## Detections conditional on location
    for(j in 1:J) {
        for(k in 1:K) {
            dist2.ux[,j,k] <- (u[,k,1]-x[j,1])^2 + (u[,k,2]-x[j,2])^2
            lambda[,j,k] <- lam0*exp(-dist2.ux[,j,k]/(2*sigma.det^2))
            y.det.N[,j,k] <- rpois(N, lambda[,j,k])
            if(any(test <- y.det.N[,j,k]>0) & verbose) {
                cat("Detection distances =", sqrt(dist2.ux[test,j,k]), "\n")
            }
        }
    }

    y.cap.N <- rbinom(N, 1, p.cap)

    if(verbose) cat("\n")
    for(j in 1:J) {
        for(k in 1:K) {
            test <- (y.det.N[,j,k]>0) & (y.cap.N==1)
            if(any(test) & verbose) {
                ## browser()
                dists <- sqrt((u[test,k,1]-x[j,1])^2 +
                              (u[test,k,2]-x[j,2])^2)
                cat("Observed detection distances =", dists, "\n")
            }
        }
    }

    ## Augmented encounter histories
    y.det.guys <- rowSums(y.det.N)>0
    y.cap.guys <- y.cap.N>0
    ## y.cap guys must come first!
    all.guys <- unique(c(which(y.cap.guys), which(y.det.guys)))
    n.guys <- length(all.guys)
    y.det <- y.det.N[all.guys,,]
    y.cap <- y.cap.N[all.guys]

    ## Observed telemetry locations
    u.data <- array(NA_real_, c(n.guys, K, 2))
    observed.telem <- seq(1, K, thin.telem)
    u.data[1:sum(y.cap.guys),observed.telem,] <-
        u[y.cap.guys,observed.telem,]

    out <- list(pars=c(N, sigma, rho, p.cap, lam0, sigma.det),
                latent=list(s=s, u=u, y.det=y.det.N, y.cap=y.cap.N),
                x=x, xlim=xlim, ylim=ylim,
                y.det=y.det, y.cap=y.cap,
                u=u.data,
                iseed=iseed)

    return(out)
}



if(1==2) {

    ## debugonce(sim.scr.move)


    tmp <- sim.scr.move(x=cam.locs, buffer=5000, K=90, ##M=250,
                        N=100, rho=0.95, sigma=600, p.cap=0.25,
                        lam0=2, sigma.det=50, thin.telem=3, verbose=TRUE)

    str(tmp)


    table(tmp$y.det)
    rowSums(tmp$y.det)
    table(rowSums(tmp$y.det))
    sum(tmp$y.cap)

    plot(tmp$x, asp=1, pch=3, xlim=tmp$xlim, ylim=tmp$ylim)
    utmp <- tmp$latent$u
    for(i in 1:nrow(utmp)) {
        arrows(utmp[i,-dim(utmp)[2],1], utmp[i,-dim(utmp)[2],2],
               utmp[i,-1,1], utmp[i,-1,2],
               length=0.05, col=ifelse(tmp$latent$y.cap[i]==1,
                                       "blue", gray(0.8)))
    }

    (guys <- which(tmp$y.cap==1 & rowSums(tmp$y.det)>0))

    guy <- guys[2]
    (guy.trap <- which(rowSums(tmp$y.det[guy,,])>0))
    (guy.det.time <- which(colSums(tmp$y.det[guy,,])>0))

    plot(tmp$u[guy,,], asp=1)##, xlim=range(tmp$x[,1]), ylim=range(tmp$x[,2]))
    arrows(tmp$u[guy,-90,1], tmp$u[guy,-90,2],
           tmp$u[guy,-1,1], tmp$u[guy,-1,2], len=0.05)
    points(tmp$x, pch=3)
    points(tmp$x[guy.trap,,drop=FALSE], pch=3, col="red", cex=3, lwd=2)
    segments(tmp$u[guy,guy.det.time,1], tmp$u[guy,guy.det.time,2],
             tmp$x[guy.trap,1], tmp$x[guy.trap,2], col="blue", lty=3, lwd=3)
    
    ls()

##    ini <- if(exists("fm1")) fm1$final.state else NULL

    ## source("../../scrmove/R/scr_move_DA_capdet_mcmc.R")

    system.time({
    fm1 <- scrMoveDA(ycap=as.matrix(tmp$y.cap),
                     ydet=tmp$y.det,
                     u=tmp$u,
                     M=225,
                     x=tmp$x,
                     ## oper=oper,
                     plotit=FALSE, 
                     n.iters=1000, n.mci=100L,
                     buffer=5000, trim=1000, nthreads=4,
                     report=100, verbose=FALSE,
                     block.rhosigma=TRUE,
                     block.lam0kappa=TRUE,
                     ## inits=ini,
                     ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                     tune=c(0.29, 3e-3, 15, 0.29, 4, 0.1, 300, 300, 3.01e-3, -0.05))
    }) ## 8500 it/hr

    mc1 <- as.mcmc(fm1$samples)

    plot(mc1, ask=TRUE)

    rejectionRate(mc1)

    effectiveSize(mc1)



    system.time({
    fm2 <- scrDA(ycap=as.matrix(tmp$y.cap),
                 ydet=tmp$y.det,
                 M=225,
                 x=tmp$x,
                 ## oper=oper,
                 plotit=FALSE, 
                 n.iters=1000, ##n.mci=100L,
                 buffer=5000, trim=1000, nthreads=4,
                 report=100, verbose=FALSE,
                 ## block.rhosigma=TRUE,
                 ## block.lam0kappa=TRUE,
                 ## inits=ini,
                 ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                 tune=c(0.01, 20, 300, -0.001))
    }) ## 8500 it/hr

    mc2 <- as.mcmc(fm2$samples)

    plot(mc2, ask=TRUE)

    rejectionRate(mc2)

    effectiveSize(mc2)
    

}








## Simulate datasets


## rho=0.55

sqrt(600^2 - 600^2*0.55^2)

set.seed(3405)
sims0.55 <- replicate(n=100, sim.scr.move(x=cam.locs, buffer=5000, K=90, 
                                          N=100, rho=0.55, sigma=600, p.cap=0.25,
                                          lam0=2, sigma.det=50, thin.telem=3),
                      simplify=FALSE)
save(sims0.55, file="sims055.gzip")


## rho=0.65

sqrt(600^2 - 600^2*0.65^2)

set.seed(340975)
sims0.65 <- replicate(n=100, sim.scr.move(x=cam.locs, buffer=5000, K=90, 
                                          N=100, rho=0.65, sigma=600, p.cap=0.25,
                                          lam0=2, sigma.det=50, thin.telem=3),
                      simplify=FALSE)
save(sims0.65, file="sims065.gzip")


## rho=0.75

sqrt(600^2 - 600^2*0.75^2)

set.seed(368405)
sims0.75 <- replicate(n=100, sim.scr.move(x=cam.locs, buffer=5000, K=90, 
                                          N=100, rho=0.75, sigma=600, p.cap=0.25,
                                          lam0=2, sigma.det=50, thin.telem=3),
                      simplify=FALSE)
save(sims0.75, file="sims075.gzip")


## rho=0.85

sqrt(600^2 - 600^2*0.85^2)

set.seed(98805)
sims0.85 <- replicate(n=100, sim.scr.move(x=cam.locs, buffer=5000, K=90, 
                                          N=100, rho=0.85, sigma=600, p.cap=0.25,
                                          lam0=2, sigma.det=50, thin.telem=3),
                      simplify=FALSE)
save(sims0.85, file="sims085.gzip")



## rho=0.95

sqrt(600^2 - 600^2*0.95^2)

set.seed(3985)
sims0.95 <- replicate(n=100, sim.scr.move(x=cam.locs, buffer=5000, K=90, 
                                          N=100, rho=0.95, sigma=600, p.cap=0.25,
                                          lam0=2, sigma.det=50, thin.telem=3),
                      simplify=FALSE)
save(sims0.95, file="sims095.gzip")






## Summarize datasets

ssfun <- function(x) {
    y.cap.guys <- x$y.cap>0
    y.det.guys <- rowSums(x$y.det)>0
    out <- cbind(n.cap.guys=sum(y.cap.guys),
                 n.det.guys=sum(y.det.guys),
                 n.both.guys=sum(y.cap.guys & y.det.guys))
    return(out)
}

ssfun(sims0.95[[1]])

summary(t(sapply(sims0.55, ssfun)))
summary(t(sapply(sims0.65, ssfun)))
summary(t(sapply(sims0.75, ssfun)))
summary(t(sapply(sims0.85, ssfun)))
summary(t(sapply(sims0.95, ssfun)))
