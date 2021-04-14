## MCMC algorithm for SCR+telemetry
## No data augmentation. Uses "semi-complete likelihood" approach instead

scrMoveDArsig <- function(ycap,  ## nx1 initial (helicopter-based) capture histories
                          ydet,  ## nxJxK camera trap encounter histories
                          u,     ## nxKx2 array of telemetry coordinates
                          M,     ## Data augmentation dimension
                          x,     ## Trap coordinates
                          oper,  ## Operational status
                          report=10, verbose=FALSE, plotit=TRUE,
                          ## block.rhosigma=FALSE,
                          random.sig=TRUE,
                          block.lam0kappa=FALSE,
                          ## tune order: rho, log.sigma, log.sigma.mu, log.sigma.sd
                          ##             lam0, kappa, p, s, u, cov(lam0,kappa)
                          tune, 
                          n.iters,
                          n.mci,
                          buffer, trim=600, nthreads,
                          su.post.keep=0,
                          pu.prob=0, ## prob of proposing u from prior
                          inits) {
    
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    ## No need to augment ydet with 0's
    ## We'll just use one zero for each augmented individual
    ## We will augment ycap though
    ydet.dims <- dim(ydet)
    n <- ydet.dims[1]
    J <- ydet.dims[2]
    K <- ydet.dims[3]

    if(missing(u)) {
        u <- array(NA_real_, c(n, K, 2))
        u.data <- FALSE
    } else {
        u.in <- u
        u.data <- TRUE
    }
    if(!u.data) u.in <- NULL

    if(nrow(ycap) != nrow(ydet))
        stop("nrow(ycap) != nrow(ydet)")
    if(nrow(ycap) != nrow(u))
        stop("nrow(ycap) != nrow(u)")
    if(J != nrow(x))
        stop("ydet and x should have the same number of traps")
    if(K != dim(u)[2])
        stop("ydet and u should have the same number of occasions")
    if(pu.prob<0 | pu.prob>1)
        stop("pu.prob must be a probability... Pr(proposing u from prior)")

    if(missing(nthreads)) {
        require(parallel)
        availableCores <- detectCores()
        nthreads <- availableCores-1
    }

    n.mci <- as.integer(n.mci)
    nthreads <- as.integer(nthreads)

    u.na <- apply(is.na(u), c(1, 2), any)
    no.u <- sum(u.na)
    has.telem <- !apply(u.na, 1, all)
    has.ydet <- apply(ydet>0, 1, any)
    has.telem.or.ydet <- has.telem | has.ydet

    xlim <- range(x[,1])+c(-buffer,buffer)
    ylim <- range(x[,2])+c(-buffer,buffer)
    Area <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1]) / 1e6 ## km-sq

    if(missing(oper)) {
        oper <- matrix(1L, J, K)
    } else {
        if((nrow(oper) != J) || (ncol(oper) != K))
            stop("oper should be a nSite x nOccasion binary matrix")
    }
    notOper <- !oper

    ## Initial values
    if(missing(inits)) {
        s <- cbind(runif(M, xlim[1], xlim[2]),
                   runif(M, ylim[1], ylim[2]))
        for(i in 1:n) {
            if(has.telem[i]) {
                s[i,] <- colMeans(u[i,,], na.rm=TRUE)
            } else {
                traps.w.dets.i <- rowSums(ydet[i,,])>0
                if(all(!traps.w.dets.i))
                    next
                s[i,] <- colMeans(x[traps.w.dets.i,,drop=FALSE])
            }
        }
        ## beta0 <- runif(1, log(0.3), log(0.4))
        ## beta0 <- runif(1, log(M/2/Area), log(M/1.5/Area)) ## log(E(Density))
        beta1 <- 0
        ## N <- n+ceiling(runif(1, 0, 100))
        psi <- runif(1, 0.2, 0.8)
        z <- c(rep(1L, n), rbinom(M-n, 1, psi))
        N <- sum(z)
        lam0 <- runif(1, 0.5, 1.5)
        if(random.sig) {
            sigma <- runif(M, 600, 700)
        } else {
            sigma <- rep(runif(1, 600, 700), M)
        }
        log.sigma <- log(sigma)
        log.sigma.mu <- mean(log.sigma)
        log.sigma.sd <- sd(log.sigma)
        rho <- runif(1, 0.94, 0.95)
        kappa <- runif(1, 20, 30)
        p <- runif(1)
        for(i in 1:n) {
            if(u.na[i,1]) {
                if(!has.telem[i]) {
                    u[i,1,] <- rnorm(2, s[i,], sigma[i])
                } else {
                    first.u.i <- min(which(!u.na[i,]))
                    u[i,1,] <- rnorm(2, u[i,first.u.i,], 10)
                }
            }
            for(k in 2:K) {
                if(!u.na[i,k])
                    next
                if(has.telem[i])
                    u[i,k,] <- rnorm(2, u[i,k-1,], 10)
                else
                    u[i,k,] <- rnorm(2, s[i,]+(u[i,k-1,]-s[i,])*rho,
                                     sqrt(sigma[i]^2-sigma[i]^2*rho^2))
            }
            for(k in 1:K) {
                dets.ik <- ydet[i,,k]>0
                if(all(!dets.ik) | !u.na[i,k])
                    next
                traps.w.dets.ik <- which(dets.ik)
                if(length(traps.w.dets.ik)>1)
                    warning("Multiple traps during a single occasion")
                u[i,k,] <- x[traps.w.dets.ik[1],]
            }
        }
    } else {
        ## beta0 <- log(inits$theta["ED"])
        psi <- inits$theta["psi"]
        lam0 <- inits$theta["lam0"]
        log.sigma.mu <- inits$theta["log.sigma.mu"]
        log.sigma.sd <- inits$theta["log.sigma.sd"]
        rho <- inits$theta["rho"]
        kappa <- inits$theta["kappa"]
        N <- inits$theta["N"]
        p <- inits$theta["p"]
        z <- inits$z
        sigma <- inits$sigma
        log.sigma <- log(sigma)
        s <- inits$s
        u <- inits$u
        .Random.seed <- inits$RNG.state
    }

    ## Functions to compute log-densities
    ## get.ld.N <- function(N, EN) return(dpois(N, EN, log=TRUE)) ## For Poisson point process
    ## get.ld.n0 <- function(n, N, q.star) {       ## Likelihood of missing n0 guys.
    ##     n0 <- N-n
    ##     if(n0<0) {
    ##         return(-Inf)
    ##     } else if(q.star<0 | q.star>1 | !is.finite(q.star)) {
    ##         warning("q.star should be a probability")
    ##         return(-Inf)
    ##     } else {
    ##         return(lgamma(N+1)-lgamma(n0+1)-lgamma(n+1)+log(q.star)*n0)
    ##     }
    ## }

##    browser()

    ## p(y|u,sigma,lam0,varsigma)
    far <- lambda <- distSq <- array(NA_real_, c(n, J, K))
    for(j in 1:J) {
        for(k in 1:K) {
            distSq[,j,k] <- (u[,k,1]-x[j,1])^2 + (u[,k,2]-x[j,2])^2
            far[,j,k] <- distSq[,j,k]>(trim^2)
            ## lambda[,j,k] <- lam0*exp(-distSq[,j,k]/(2*kappa^2))
            lambda[,j,k] <- ifelse(far[,j,k], 0,
                                   lam0*exp(-distSq[,j,k]/(2*kappa^2))*oper[j,k])
        }
    }

##    browser()
    
    ## p(ydet) for i=1:n
    ld.ydet <- ld.ydet.cand <- get_ld_y_nK(ydet, lambda, nthreads) 
##    ld.ydet <- ld.ydet.cand <- get_ld_y_nK_stats(ydet, lambda, nthreads) ## avoid R::dpois
    ld.ydet.sum <- ld.ydet.cand.sum <- sum(ld.ydet)

    ## p(ydet) for i=(n+1):M
    ld.ydet.aug <- rep(0, M)
    ld.ydet.aug[(n+1):M] <- 0
    for(i in (n+1):M) {
##        if(z[i]==0) 
##            ld.ydet.aug[i] <- 0
##        else
        if(z[i]==1) {
            ld.ydet.aug[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                                              x, K, sigma[i], rho, lam0, kappa,
                                              trim^2, notOper, nthreads)
        }
    }
    ld.ydet.aug.cand <- ld.ydet.aug

    ## Augment ycap
    ycap.aug <- c(ycap, rep(0L, M-n))
    ld.ycap <- dbinom(ycap.aug, 1, p*z, log=TRUE)
    ld.ycap.cand <- ld.ycap
    ld.ycap.sum <- sum(ld.ycap)

    ## p(N|EN(beta0))
    logpixArea <- 0 ## Replace with pixel area if raster covariates are added
    ## logmu <- beta0 + logpixArea ## Could add covariates here
    ##    EN <- sum(exp(logmu))
    ## ED <- exp(logmu)
    ## EN <- ED*Area
    EN <- M*psi
    ED <- EN/Area
    logEN <- log(EN)
    ## ld.N <- get.ld.N(N, EN)

    ld.kappa <- 0 #dgamma(kappa, 1, 1, log=TRUE)

    ## p(n0|N,q.star)
    if(n.mci %% 2 != 0)
        stop("n.mci must be even")
    ## q.star.cam <- get_qstar_cam(n.mci, x, K, sigma, rho, lam0, kappa,
    ##                             buffer, (trim)^2, notOper, nthreads)
    ## q.star.gps <- 1-p
    ## q.star <- q.star.gps*q.star.cam
    ## ld.n0 <- get.ld.n0(n, N, q.star)

    ## p(z|psi)
    ld.z <- ld.z.cand <- dbinom(z, 1, psi, log=TRUE)
    
    ##browser()

    ## p(s|beta)
    ## NOTE: Assuming homogeneous point process for now: s(i)~Unif(s)
    ## ld.s <- ld.s.cand <- logmu - logEN
    ## ld.s.sum <- ld.s.cand.sum <- sum(ld.s)
    ld.s <- ld.s.cand <- 0

    ld.u <- matrix(NA, n, K)
    ## ld.u[,1] <- dnorm(u[,1,1], s[,1], sigma, log=TRUE) +
    ##     dnorm(u[,1,2], s[,2], sigma, log=TRUE)
    ld.u[,1] <- rowSums(dnorm(u[,1,], s[1:n,], sigma[i], log=TRUE))
    for(k in 2:K) {
        ld.u[,k] <-
            rowSums(dnorm(u[,k,], s[1:n,]+(u[,k-1,]-s[1:n,])*rho,
                          sqrt(sigma[1:n]^2-sigma[1:n]^2*rho^2), log=TRUE))
    }
    ld.u.cand <- ld.u
    ld.u.sum <- sum(ld.u)

    u.cand <- u 

    ## deviance.ydet <- -2*ld.ydet.sum 
    deviance.ydet <- -2*(ld.ydet.sum + sum(ld.ydet.aug))
    deviance.ycap <- -2*ld.ycap.sum 
    deviance.u <- -2*sum(ld.u[!u.na])

    ordr <- sample.int(n)
    icols1 <- topo.colors(n)[ordr]
    icols2 <- topo.colors(n, alpha=0.3)[ordr]

    keep.su.post <- su.post.keep>0
    if(keep.su.post) {
        su.keep.seq <- seq(1, n.iters, su.post.keep)
        niter.su.keep <- length(su.keep.seq)
        sigma.post <- matrix(NA_real_, M, niter.su.keep)
        s.post <- array(NA_real_, c(M, 2, niter.su.keep))
        u.post <- array(NA_real_, c(n, K, 2, niter.su.keep))
        su.counter <- 1
    } else {
        sigma.post <- NULL
        s.post <- NULL
        u.post <- NULL
    }

    ## if(block.rhosigma) {
    ##     if(!require(mvtnorm))
    ##         stop("must load 'mvtnorm' package for block updating")
    ##     proposal.vcov.rhosigma <- matrix(c(tune[2]^2, tune[9], tune[9], tune[3]^2), 2)
    ##     ev <- eigen(proposal.vcov.rhosigma, symmetric = TRUE)
    ##     if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
    ##         stop("vcov proposal matrix for rho+sigma is numerically not positive semidefinite. Change tuning parameters.")
    ##     }
    ## }
    if(block.lam0kappa) {
        if(!require(mvtnorm))
            stop("must load 'mvtnorm' package for block updating")
        cov.tune <- tune[length(tune)]
        proposal.vcov.lam0kappa <- matrix(c(tune[5]^2, cov.tune, cov.tune, tune[6]^2), 2)
        ev <- eigen(proposal.vcov.lam0kappa, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            stop("vcov proposal matrix for lam0+kappa is numerically not positive semidefinite. Change tuning parameters.")
        }
    }
    

    out <- matrix(NA, nrow=n.iters, ncol=12)
    colnames(out) <- c("ED", "rho", "log.sigma.mu", "log.sigma.sd", "lam0", "kappa",
                       "psi", "N", "p", "deviance.ycap", "deviance.ydet",
                       "deviance.u")

##    browser()

    if(report>0)
        cat("\nstarting values =",
            round(c(ED, rho, log.sigma.mu, log.sigma.sd, lam0, kappa, psi, N, p,
                    deviance.ycap, deviance.ydet, deviance.u), 3), "\n\n")

    for(iter in 1:n.iters) {

        if(iter %% report == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", round(out[iter-1,], 3), "\n")
            cat("   summary(sigma) =", round(summary(sigma), 1), "\n")
            if(plotit) {
                plot(s, xlim=xlim, ylim=ylim, asp=1, cex=0.3,
                     pch=16,
                     col=icols1)
                points(x, pch="+")
                rect(xlim[1], ylim[1], xlim[2], ylim[2])
                for(i in 1:n) {
                    points(u[i,,], cex=0.2, pch=16,
                           col=icols1[i])
                }
                for(i in which(has.telem)) {
                    suppressWarnings(
                        arrows(u[i,-K,1], u[i,-K,2],
                               u[i,-1,1], u[i,-1,2],
                               length=0.02, col=icols2[i]))
                }
            }
        }

        ## if(verbose)
        ##     cat("   Sampling p(beta0|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ## ## Sample from p(beta0|.)
        ## beta0.cand <- rnorm(1, beta0, tune[1])
        ## logmu.cand <- beta0.cand + logpixArea
        ## ED.cand <- exp(logmu.cand)
        ## EN.cand <- ED.cand*Area
        ## ld.N.cand <- get.ld.N(N, EN.cand)
        ## ld.beta0.cand <- ld.beta0 <- 0
        ## ## if(runif(1) < exp((ld.N.cand + ld.beta0.cand) -
        ## ##                   (ld.N + ld.beta0))) {
        ## mhr.beta0 <- exp((ld.N.cand + ld.beta0.cand) -
        ##                  (ld.N + ld.beta0))
        ## if(!is.finite(mhr.beta0)) stop("beta0")
        ## if(runif(1) < mhr.beta0) {
        ##     ## ld.s <- ld.s.cand
        ##     ## ld.s.sum <- ld.s.cand.sum
        ##     ld.N <- ld.N.cand
        ##     beta0 <- beta0.cand
        ##     ED <- ED.cand
        ##     EN <- EN.cand
        ##     logEN <- log(EN)
        ## }


        ## if(block.rhosigma) {
        ##     if(verbose)
        ##         cat("   Sampling p(rho,sigma|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
            
        ##     ## Sample from p(rho,sigma|.)
        ##     rhosigma.cand <- rmvnorm(1, c(rho, sigma), proposal.vcov.rhosigma)
        ##     rho.cand <- rhosigma.cand[1]
        ##     sigma.cand <- rhosigma.cand[2]
        ##     if(rho.cand>0 & rho.cand<1 & sigma.cand>0) {
        ##         ld.u.cand[,1] <- rowSums(dnorm(u[,1,], s[1:n,], sigma.cand, log=TRUE))
        ##         ld.u.cand[,-1] <-
        ##             dnorm(u[,-1,1], s[1:n,1]+(u[,-K,1]-s[1:n,1])*rho.cand,
        ##                   sqrt(sigma.cand^2-sigma.cand^2*rho.cand^2), log=TRUE) +
        ##             dnorm(u[,-1,2], s[1:n,2]+(u[,-K,2]-s[1:n,2])*rho.cand,
        ##                   sqrt(sigma.cand^2-sigma.cand^2*rho.cand^2), log=TRUE)
        ##         ld.u.cand.sum <- sum(ld.u.cand)
        ##         ld.ydet.aug.cand[(n+1):M] <- 0
        ##         for(i in (n+1):M) {
        ##             if(z[i]==1) {
        ##                 ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
        ##                                                        x, K, sigma.cand, rho.cand, lam0, kappa,
        ##                                                        trim^2, notOper, nthreads)
        ##             }
        ##         }
        ##         ld.rho <- ld.rho.cand <- 0
        ##         ld.sigma <- ld.sigma.cand <- 0
        ##         mhr.rhosigma <-  exp((ld.u.cand.sum + sum(ld.ydet.aug.cand) + ld.rho.cand +
        ##                               ld.sigma.cand) -
        ##                              (ld.u.sum + sum(ld.ydet.aug) + ld.rho + ld.sigma))
        ##         if(!is.finite(mhr.rhosigma)) {
        ##             browser()
        ##             ## stop("mhr.rhosigma = ", mhr.rhosigma)
        ##         }
        ##         if(runif(1) < mhr.rhosigma) {
        ##             ld.u <- ld.u.cand
        ##             ld.u.sum <- ld.u.cand.sum
        ##             ld.ydet.aug <- ld.ydet.aug.cand
        ##             rho <- rho.cand
        ##             sigma <- sigma.cand
        ##         }
        ##     }
        ## } else {
        if(verbose)
            cat("   Sampling p(rho|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
        
        ## Sample from p(rho|.)
        rho.cand <- rnorm(1, rho, tune[1])
        if(rho.cand>0 & rho.cand<1) { 
            ## ld.u.cand <- ld.u  ## Shouldn't be necessary, Paul likes this.
            ld.u.cand[,1] <- ld.u[,1]
            ld.u.cand[,-1] <-
                dnorm(u[,-1,1], s[1:n,1]+(u[,-K,1]-s[1:n,1])*rho.cand,
                      sqrt(sigma[1:n]^2-sigma[1:n]^2*rho.cand^2), log=TRUE) +
                dnorm(u[,-1,2], s[1:n,2]+(u[,-K,2]-s[1:n,2])*rho.cand,
                      sqrt(sigma[1:n]^2-sigma[1:n]^2*rho.cand^2), log=TRUE)
            ld.u.cand.sum <- sum(ld.u.cand)
            ld.ydet.aug.cand[(n+1):M] <- 0
            for(i in (n+1):M) {
                if(z[i]==1) {
                    ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                                                           x, K, sigma[i], rho.cand, lam0, kappa,
                                                           trim^2, notOper, nthreads)
                }
            }
            ld.rho <- ld.rho.cand <- 0
            mhr.rho <-  exp((ld.u.cand.sum + sum(ld.ydet.aug.cand) + ld.rho.cand) -
                            (ld.u.sum + sum(ld.ydet.aug) + ld.rho))
            if(runif(1) < mhr.rho) {
                ld.u <- ld.u.cand
                ld.u.sum <- ld.u.cand.sum
                ld.ydet.aug <- ld.ydet.aug.cand
                rho <- rho.cand
            }
        }
        
        
        if(verbose)
            cat("   Sampling p(sigma[i]|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
        
##        browser()

        ## Sample from p(log(sigma[i])|.)
        if(random.sig) {
            log.sigma.cand <- rnorm(M, log.sigma, tune[2])
        } else {
            log.sigma.cand <- log.sigma+rnorm(1, 0, tune[2])
        }
        sigma.cand <- exp(log.sigma.cand)
        ld.u.cand[,1] <- rowSums(dnorm(u[,1,], s[1:n,], sigma.cand[1:n], log=TRUE))
        ld.u.cand[,-1] <-
            dnorm(u[,-1,1], s[1:n,1]+(u[,-K,1]-s[1:n,1])*rho,
                  sqrt(sigma.cand[1:n]^2-sigma.cand[1:n]^2*rho^2), log=TRUE) +
            dnorm(u[,-1,2], s[1:n,2]+(u[,-K,2]-s[1:n,2])*rho,
                  sqrt(sigma.cand[1:n]^2-sigma.cand[1:n]^2*rho^2), log=TRUE)
        ## ld.u.cand.sum <- sum(ld.u.cand)
        ld.ydet.aug.cand[(n+1):M] <- 0
        for(i in (n+1):M) {
            if(z[i]==1) {
                ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                                                       x, K, sigma.cand[i], rho, lam0, kappa,
                                                       trim^2, notOper, nthreads)
            }
        }
        if(random.sig) {
            ld.sigma <- dnorm(log.sigma, log.sigma.mu, log.sigma.sd, log=TRUE)
            ld.sigma.cand <- dnorm(log.sigma.cand, log.sigma.mu, log.sigma.sd, log=TRUE)
        } else {
            ld.sigma <- ld.sigma.cand <- 0
        }
        ## Tack on the 0's to make dimensions align
        mhr.sigma <-  exp((c(rowSums(ld.u.cand), rep(0, M-n)) + ld.ydet.aug.cand + ld.sigma.cand) -
                          (c(rowSums(ld.u), rep(0, M-n)) + ld.ydet.aug + ld.sigma))
        accept.sigma <- runif(M) < mhr.sigma
##        if(runif(1) < mhr.sigma) {
        ld.u[accept.sigma[1:n],] <- ld.u.cand[accept.sigma[1:n],]
        ld.ydet.aug[accept.sigma] <- ld.ydet.aug.cand[accept.sigma]
        sigma[accept.sigma] <- sigma.cand[accept.sigma]
        log.sigma <- log(sigma)
##        }
        ld.u.sum <- sum(ld.u)


        ## Sample from p(logsig.mu|.) and p(logsig.sd|.)
        ## Could put conjugate inverse gamma prior on variance and sample directly
        if(random.sig) {
            log.sigma.mu.cand <- rnorm(1, log.sigma.mu, tune[3])
            ld.sigma <- dnorm(log.sigma, log.sigma.mu, log.sigma.sd, log=TRUE)
            ld.sigma.cand <- dnorm(log.sigma, log.sigma.mu.cand, log.sigma.sd, log=TRUE)
            if(runif(1) < exp(sum(ld.sigma.cand) - sum(ld.sigma))) {
                log.sigma.mu <- log.sigma.mu.cand
                ld.sigma <- ld.sigma.cand
            }
        
            log.sigma.sd.cand <- rnorm(1, log.sigma.sd, tune[4])
            if(log.sigma.sd.cand>0) {
                ld.sigma.cand <- dnorm(log.sigma, log.sigma.mu, log.sigma.sd.cand, log=TRUE)
                if(runif(1) < exp(sum(ld.sigma.cand) - sum(ld.sigma))) {
                    log.sigma.sd <- log.sigma.sd.cand
                    ld.sigma <- ld.sigma.cand
                }
            }
        }
            
        if(block.lam0kappa) {
            if(verbose)
                cat("   Sampling p(lam0,kappa|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
            
            ## Sample from p(lam0,kappa|.)
            lam0kappa.cand <- rmvnorm(1, c(lam0, kappa), proposal.vcov.lam0kappa)
            lam0.cand <- lam0kappa.cand[1]
            kappa.cand <- lam0kappa.cand[2]
            if(lam0.cand > 0 & kappa.cand>0) {
                lambda.cand <- get_lambda(lam0.cand, distSq, kappa.cand, far, notOper, nthreads) 
                ld.ydet.cand <- get_ld_y_nK(ydet, lambda.cand, nthreads) 
##                ld.ydet.cand <- get_ld_y_nK_stats(ydet, lambda.cand, nthreads) ## avoid R::dpois
                ld.ydet.cand.sum <- sum(ld.ydet.cand)
                ld.ydet.aug.cand[(n+1):M] <- 0
                for(i in (n+1):M) {
                    if(z[i]==1) {
                        ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                                                               x, K, sigma[i], rho,
                                                               lam0.cand, kappa.cand,
                                                               trim^2, notOper, nthreads)
                    }
                }
                ld.lam0 <- ld.lam0.cand <- 0
                ld.kappa <- ld.kappa.cand <- 0
                mhr.lam0kappa  <- exp((ld.ydet.cand.sum + sum(ld.ydet.aug.cand) + ld.lam0.cand +
                                       ld.kappa.cand) -
                                      (ld.ydet.sum + sum(ld.ydet.aug) + ld.lam0 + ld.kappa))
                if(!is.finite(mhr.lam0kappa)) {
##                    browser()
                    stop("mhr.lam0kappa = ", mhr.lam0kappa)
                }
                if(runif(1) < mhr.lam0kappa) {
                    ld.ydet <- ld.ydet.cand
                    ld.ydet.sum <- ld.ydet.cand.sum
                    ld.ydet.aug <- ld.ydet.aug.cand
                    lambda <- lambda.cand
                    lam0 <- lam0.cand
                    kappa <- kappa.cand
                }
            }
        } else {
            if(verbose)
                cat("   Sampling p(lam0|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
            
            ## Sample from p(lam0|.)
            lam0.cand <- rnorm(1, lam0, tune[5])
            if(lam0.cand > 0) {
                lambda.cand <- get_lambda(lam0.cand, distSq, kappa, far, notOper, nthreads) 
                ld.ydet.cand <- get_ld_y_nK(ydet, lambda.cand, nthreads) 
##                ld.ydet.cand <- get_ld_y_nK_stats(ydet, lambda.cand, nthreads) ## avoid R::dpois
                ld.ydet.cand.sum <- sum(ld.ydet.cand)
                ld.ydet.aug.cand[(n+1):M] <- 0
                for(i in (n+1):M) {
                    if(z[i]==1) {
                        ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                                                               x, K, sigma[i], rho,
                                                               lam0.cand, kappa,
                                                               trim^2, notOper, nthreads)
                    }
                }
                ld.lam0 <- ld.lam0.cand <- 0
                mhr.lam0  <- exp((ld.ydet.cand.sum + sum(ld.ydet.aug.cand) + ld.lam0.cand) -
                                 (ld.ydet.sum + sum(ld.ydet.aug) + ld.lam0))
                if(runif(1) < mhr.lam0) {
                    ld.ydet <- ld.ydet.cand
                    ld.ydet.sum <- ld.ydet.cand.sum
                    ld.ydet.aug <- ld.ydet.aug.cand
                    lambda <- lambda.cand
                    lam0 <- lam0.cand
                }
            }

            if(verbose)
                cat("   Sampling p(kappa|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

            ## Sample from p(kappa|.)
            kappa.cand <- rnorm(1, kappa, tune[6]) ##sqrt(1/(2*lam1.cand))
            if(kappa.cand>0) { 
                lambda.cand <- get_lambda(lam0, distSq, kappa.cand, far, notOper, nthreads) 
                ld.ydet.cand <- get_ld_y_nK(ydet, lambda.cand, nthreads) 
                ## ld.ydet.cand <- get_ld_y_nK_stats(ydet, lambda.cand, nthreads) ## avoid R::dpois
                ld.ydet.cand.sum <- sum(ld.ydet.cand)
                ld.ydet.aug.cand[(n+1):M] <- 0
                for(i in (n+1):M) {
                    if(z[i]==1) {
                        ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                                                               x, K, sigma[i], rho,
                                                               lam0, kappa.cand,
                                                               trim^2, notOper, nthreads)
                    }
                }
                ld.kappa.cand <- 0 ##dgamma(kappa.cand, 1, 1, log=TRUE)
                mhr.kappa <- exp((ld.ydet.cand.sum + sum(ld.ydet.aug.cand) + ld.kappa.cand) -
                                 (ld.ydet.sum + sum(ld.ydet.aug) + ld.kappa))
                if(runif(1) < mhr.kappa) {
                    ld.ydet <- ld.ydet.cand
                    ld.ydet.sum <- ld.ydet.cand.sum
                    ld.ydet.aug <- ld.ydet.aug.cand
                    lambda <- lambda.cand
                    kappa <- kappa.cand
                    ld.kappa <- ld.kappa.cand
                }
            }
        }


        if(verbose)
            cat("   Sampling p(N|.): ", format(Sys.time(), "%H:%M:%S"), "\n")


        ## Sample from p(psi|.)
        psi <- rbeta(1, 1+N, 1+M-N)
        EN <- M*psi
        logEN <- log(EN)
        ED <- EN/Area

        
        if(verbose)
            cat("   Sampling p(z|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ## Sample from p(z|.)
        ld.ydet.aug.cand[(n+1):M] <- 0
        for(i in (n+1):M) {
            z.cand.i <- 1-z[i]
            ## Must compute ld.z because psi has changed
            ld.z[i] <- dbinom(z[i], 1, psi, log=TRUE)
            ld.z.cand[i] <- dbinom(z.cand.i, 1, psi, log=TRUE)
            if(z.cand.i==1) {
                ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                                                       x, K, sigma[i], rho,
                                                       lam0, kappa,
                                                       trim^2, notOper, nthreads)
            }
            ld.ycap.cand[i] <- dbinom(ycap.aug[i], 1, p*z.cand.i, log=TRUE)
            mhr.z <- exp((ld.z.cand[i]+ld.ydet.aug.cand[i]+ld.ycap.cand[i])-
                         (ld.z[i]+ld.ydet.aug[i]+ld.ycap[i]))
            if(!is.finite(mhr.z)) {
                ## browser()
                stop("mhr.z = ", mhr.z)
            }
            if(runif(1) < mhr.z) {
                z[i] <- z.cand.i
                ld.ydet.aug[i] <- ld.ydet.aug.cand[i]
                ld.ycap[i] <- ld.ycap.cand[i]
            }
        }
        N <- sum(z)

        ## Sample from p(p|.)
        if(verbose)
            cat("   Sampling p(p|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
        p.cand <- rnorm(1, p, tune[7])
        if(p.cand>0 & p.cand<1) {
            ld.ycap <- dbinom(ycap.aug, 1, p*z, log=TRUE)
            ld.ycap.cand <- dbinom(ycap.aug, 1, p.cand*z, log=TRUE)
            mhr.p <- exp(sum(ld.ycap.cand) - sum(ld.ycap))
            if(!is.finite(mhr.p)) {
                stop("mhr.p = ", mhr.p)
                ## browser()
            }
            if(runif(1) < mhr.p) {
                p <- p.cand
                ld.ycap <- ld.ycap.cand
            }
        }


        ## Sample from p(s|.)
        if(verbose)
            cat("   Sampling p(s|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        sn.ups <- 0
        sM.ups <- 0
        for(i in 1:n) {
            s1.cand <- rnorm(1, s[i,1], tune[8])
            s2.cand <- rnorm(1, s[i,2], tune[8])
            if(s1.cand < xlim[1] | s1.cand > xlim[2] | s2.cand < ylim[1] | s2.cand > ylim[2])
                next
            s.cand <- c(s1.cand, s2.cand)
            ld.u.cand[i,] <- ld.u[i,]
            ld.u.cand[i,1] <- sum(dnorm(u[i,1,], s.cand, sigma[i], log=TRUE))
            ld.u.cand[i,-1] <-
                dnorm(u[i,-1,1], s1.cand+(u[i,-K,1]-s1.cand)*rho,
                      sqrt(sigma[i]^2-sigma[i]^2*rho^2), log=TRUE) +
                dnorm(u[i,-1,2], s2.cand+(u[i,-K,2]-s2.cand)*rho,
                      sqrt(sigma[i]^2-sigma[i]^2*rho^2), log=TRUE)
            ## TODO: Add ld.s after implementing spatial variation in density
            mhr.s <-  exp(sum(ld.u.cand[i,]) - sum(ld.u[i,]))
            if(!is.finite(mhr.s)) {
                stop("mhr.s = ", mhr.s)
                ## browser()
            }
            if(runif(1) < mhr.s) {
                s[i,] <- s.cand
                ld.u[i,] <- ld.u.cand[i,]
                sn.ups <- sn.ups+1
            }
        }
        ld.ydet.aug.cand[(n+1):M] <- 0
        for(i in (n+1):M) {
            s1.cand <- rnorm(1, s[i,1], tune[8]*10)
            s2.cand <- rnorm(1, s[i,2], tune[8]*10)
            if(s1.cand < xlim[1] | s1.cand > xlim[2] | s2.cand < ylim[1] | s2.cand > ylim[2])
                next
            s.cand <- c(s1.cand, s2.cand)
            if(z[i]==1) {
                ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s1.cand, s2.cand, x, K, sigma[i],
                                                       rho, lam0, kappa, trim^2,
                                                       notOper, nthreads)
            }
            ## TODO: Add ld.s after implementing spatial variation in density
            mhr.s <-  exp(ld.ydet.aug.cand[i] - ld.ydet.aug[i])
            if(!is.finite(mhr.s)) {
                stop("mhr.s = ", mhr.s)
                ## browser()
            }
            if(runif(1) < mhr.s) {
                s[i,] <- s.cand
                ld.ydet.aug[i] <- ld.ydet.aug.cand[i]
                sM.ups <- sM.ups+1
            }
        }


        ## Sample from p(u|.)
        if(verbose)
            cat("   Sampling p(uk|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ukups <- 0
        ## distSq.cand.k <- matrix(NA, n, J)
        ## ld.u.cand <- ld.u
        ## Random walk MH does poorly when most locs are unobserved
        propose.u.from.prior <- c(FALSE, runif(K-1)<pu.prob)
        no.prior.props <- all(!propose.u.from.prior)
        u.cand[] <- rnorm(n*K*2, u, tune[9])
        if(no.prior.props) {
            distSq.cand <- get_distSq_ux(u.cand, x, nthreads=nthreads)
            far.cand <- distSq.cand>(trim^2)
            lambda.cand <- get_lambda(lam0, distSq.cand, kappa, far.cand, notOper, nthreads)
            ld.ydet.cand <- get_ld_y_nK(ydet, lambda.cand, nthreads) 
            ## ld.ydet.cand <- get_ld_y_nK_stats(ydet, lambda.cand, nthreads) ## avoid R::dpois
        }
        for(k in 1:K) {
            if(k==1) {
                u.cand.k <- u.cand[,k,] ## Always random walk for k=1
                ld.u.cand[,k] <- rowSums(dnorm(u.cand.k, s[1:n,], sigma[1:n], log=TRUE))
            } else {
                if(propose.u.from.prior[k]) {
                    u.cand.k <- cbind(
                        rnorm(n, s[1:n,1]+(u[,k-1,1]-s[1:n,1])*rho,
                              sqrt(sigma[1:n]^2-sigma[1:n]^2*rho^2)),
                        rnorm(n, s[1:n,2]+(u[,k-1,2]-s[1:n,2])*rho,
                              sqrt(sigma[1:n]^2-sigma[1:n]^2*rho^2)))
                } else {
                    u.cand.k <- u.cand[,k,]
                }
                ld.u.cand[,k] <- rowSums(dnorm(u.cand.k, s[1:n,]+(u[,k-1,]-s[1:n,])*rho,
                                               sqrt(sigma[1:n]^2-sigma[1:n]^2*rho^2), log=TRUE))
            }
            if(k<K) {
                ld.u.cand[,k+1] <- rowSums(dnorm(u[,k+1,], s[1:n,]+(u.cand.k-s[1:n,])*rho,
                                                 sqrt(sigma[1:n]^2-sigma[1:n]^2*rho^2), log=TRUE))
            }
            ## Must do this for all k if there is at least one prior proposal
            if(!no.prior.props) { 
                for(j in 1:J) {
                    distSq.cand[,j,k] <- (u.cand.k[,1]-x[j,1])^2 + (u.cand.k[,2]-x[j,2])^2
                    far.cand[,j,k] <- distSq.cand[,j,k]>(trim^2)
                    lambda.cand[,j,k] <- lam0*exp(-distSq.cand[,j,k]/(2*kappa^2))*oper[j,k] *
                        (!far.cand[,j,k])
                }
                ld.ydet.cand[,k] <- rowSums(dpois(ydet[,,k], lambda.cand[,,k], log=TRUE))
            }
            if(propose.u.from.prior[k]) {
                ## Prior probs cancel out, so only consider p(u[,k+1]) and p(ydet)
                accept.u <- (runif(n) < exp((ld.u.cand[,min(k+1,K)]*(k<K) +
                                             ld.ydet.cand[,k]) -
                                            (ld.u[,min(k+1,K)]*(k<K) +
                                             ld.ydet[,k]))) & u.na[,k] ## Ignore observed u's
            } else {
                ## For RW, must consider p(u[,k]), p(u[,k+1]), and p(ydet)
                accept.u <- (runif(n) < exp((rowSums(ld.u.cand[,k:min(k+1,K),drop=FALSE]) +
                                             ld.ydet.cand[,k]) -
                                            (rowSums(ld.u[,k:min(k+1,K),drop=FALSE]) +
                                             ld.ydet[,k]))) & u.na[,k] ## Ignore observed u's
            }
            if(any(!is.finite(accept.u))) {
                ## browser()
                stop("u update")
            }
            if(any(accept.u)) {
                lambda[accept.u,,k] <- lambda.cand[accept.u,,k]
                distSq[accept.u,,k] <- distSq.cand[accept.u,,k]
                far[accept.u,,k] <- far.cand[accept.u,,k]
                ## lambda[accept.u,,k] <- lambda.cand.k[accept.u,]
                ## distSq[accept.u,,k] <- distSq.cand.k[accept.u,]
                ## far[accept.u,,k] <- !notfar.cand.k[accept.u,]
                ld.ydet[accept.u,k] <- ld.ydet.cand[accept.u,k]
                ld.u[accept.u,k] <- ld.u.cand[accept.u,k]
                if(k<K)
                    ld.u[accept.u,k+1] <- ld.u.cand[accept.u,k+1]
                u[accept.u,k,] <- u.cand.k[accept.u,]
                ukups <- ukups+sum(accept.u)
            }
        }
        

        ld.ydet.sum <- sum(ld.ydet)
        ld.u.sum <- sum(ld.u)

        deviance.ydet <- -2*(ld.ydet.sum + sum(ld.ydet.aug))
        deviance.ycap <- -2*sum(ld.ycap)
        deviance.u <- -2*sum(ld.u[!u.na])
        dev <- deviance.ydet + deviance.ycap + deviance.u

        if(iter %% report == 0) {
            cat("   Acceptance rates\n")
            cat("     sigma =", sum(accept.sigma)/M, "\n")
            cat("     s.n =", sn.ups/n, "\n")
            cat("     s.M =", sM.ups/(M-n), "\n")
            cat("     uk =", ukups/no.u, "\n")
        }

        out[iter,] <- c(ED, rho, log.sigma.mu, log.sigma.sd, lam0, kappa, psi, N,
                        p, deviance.ycap, deviance.ydet, deviance.u)

        if(keep.su.post) {
            if(iter %in% su.keep.seq) {
                sigma.post[,su.counter] <- sigma
                s.post[,,su.counter] <- s
                u.post[,,,su.counter] <- u
                su.counter <- su.counter+1
            }
        }
    }

    return(list(samples=out,
                sigma.post=sigma.post,
                s.post=s.post,
                u.post=u.post,
                settings=list(tune=tune, trim=trim, nthreads=nthreads,
                              buffer=buffer, su.post.keep=su.post.keep,
##                              block.rhosigma=block.rhosigma,
                              block.lam0kappa=block.lam0kappa),
                final.state=list(z=z, s=s, u=u, sigma=sigma,
                                 theta=out[iter,],
                                 RNG.state=.Random.seed),
                data=list(ydet=ydet, ycap=ycap, M=M, x=x, oper=oper,
                          u=u.in),
                seed=seed))
}





