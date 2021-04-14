## MCMC algorithm for SCR data, no telemetry data
## Data augmentation

scrDA <- function(ycap,  ## nx1 initial (helicopter-based) capture histories
                  ydet,  ## nxJxK camera trap encounter histories
                  M,     ## Data augmentation dimension
                  x,     ## Trap coordinates
                  oper,  ## Operational status
                  report=10, verbose=FALSE, plotit=TRUE,
                  ## block.lam0sigma=FALSE,
                  tune, 
                  n.iters,
                  buffer, trim=3000,
                  nthreads,
                  s.post.keep=0,
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

    if(nrow(ycap) != nrow(ydet))
        stop("nrow(ycap) != nrow(ydet)")
    if(J != nrow(x))
        stop("ydet and x should have the same number of traps")

    if(missing(nthreads)) {
        require(parallel)
        availableCores <- detectCores()
        nthreads <- availableCores-1
    }

    ## n.mci <- as.integer(n.mci)
    nthreads <- as.integer(nthreads)

    ## has.ydet <- apply(ydet>0, 1, any)

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
            traps.w.dets.i <- rowSums(ydet[i,,])>0
            if(all(!traps.w.dets.i))
                next
            s[i,] <- colMeans(x[traps.w.dets.i,,drop=FALSE])
        }
        ## beta0 <- runif(1, log(0.3), log(0.4))
        ## beta0 <- runif(1, log(M/2/Area), log(M/1.5/Area)) ## log(E(Density))
        ## beta1 <- 0
        ## N <- n+ceiling(runif(1, 0, 100))
        psi <- runif(1, 0.2, 0.8)
        z <- c(rep(1L, n), rbinom(M-n, 1, psi))
        N <- sum(z)
        lam0 <- runif(1, 0.001, 0.002)
        sigma <- runif(1, 600, 700)
        p <- runif(1)
    } else {
        ## beta0 <- log(inits$theta["ED"])
        psi <- inits$theta["psi"]
        lam0 <- inits$theta["lam0"]
        sigma <- inits$theta["sigma"]
        ## rho <- inits$theta["rho"]
        N <- inits$theta["N"]
        p <- inits$theta["p"]
        z <- inits$z
        s <- inits$s
        .Random.seed <- inits$RNG.state
    }

##    browser()

    ## p(y|u,sigma,lam0,varsigma)
    distSq <- array(NA_real_, c(M, J))
    far <- matrix(FALSE, M, J)
    lambda <- array(NA_real_, c(M, J, K))
    for(j in 1:J) {
        distSq[,j] <- (s[,1]-x[j,1])^2 + (s[,2]-x[j,2])^2
        far[,j] <- distSq[,j]>(trim^2)
        for(k in 1:K) {
            lambda[,j,k] <- lam0*exp(-distSq[,j]/(2*sigma^2))*oper[j,k]*!far[,j]
        }
    }
    lambda.cand <- lambda
    
    ## p(ydet) for i=1:n
##    ld.ydet <- ld.ydet.cand <- get_ld_y_nK(ydet, lambda, nthreads) 
    ##    ld.ydet <- ld.ydet.cand <- get_ld_y_nK_stats(ydet, lambda, nthreads) ## avoid R::dpois
    ld.ydet <- numeric(M)
    ld.ydet[1:n] <- rowSums(dpois(ydet[1:n,,], lambda[1:n,,]*z[1:n], log=TRUE))
    ld.ydet[(n+1):M] <- -rowSums(lambda[(n+1):M,,])*z[(n+1):M]
    ld.ydet.cand <- ld.ydet
    ld.ydet.sum <- ld.ydet.cand.sum <- sum(ld.ydet)

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

    ## p(z|psi)
    ld.z <- ld.z.cand <- dbinom(z, 1, psi, log=TRUE)
    
    ##browser()

    ## p(s|beta)
    ## NOTE: Assuming homogeneous point process for now: s(i)~Unif(s)
    ## ld.s <- ld.s.cand <- logmu - logEN
    ## ld.s.sum <- ld.s.cand.sum <- sum(ld.s)
    ld.s <- ld.s.cand <- 0

    deviance.ydet <- -2*ld.ydet.sum 
##    deviance.ydet <- -2*(ld.ydet.sum + sum(ld.ydet.aug))
    deviance.ycap <- -2*ld.ycap.sum 

    ordr <- sample.int(n)
    icols1 <- topo.colors(n)[ordr]
    icols2 <- topo.colors(n, alpha=0.3)[ordr]

    keep.s.post <- s.post.keep>0
    if(keep.s.post) {
        s.keep.seq <- seq(1, n.iters, s.post.keep)
        niter.s.keep <- length(s.keep.seq)
        s.post <- array(NA_real_, c(M, 2, niter.s.keep))
        ld.ydet.post <- matrix(NA_real_, M, niter.s.keep)
        s.counter <- 1
    } else {
        s.post <- NULL
        ld.ydet.post <- NULL
    }

##    if(block.lam0sigma) {
        if(!require(mvtnorm))
            stop("must load 'mvtnorm' package for block updating")
        cov.tune <- tune[length(tune)]
        proposal.vcov.lam0sigma <- matrix(c(tune[1]^2, cov.tune, cov.tune, tune[2]^2), 2)
        ev <- eigen(proposal.vcov.lam0sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            stop("vcov proposal matrix for lam0+sigma is numerically not positive semidefinite. Change tuning parameters.")
        }
##    }

    out <- matrix(NA, nrow=n.iters, ncol=8)
    colnames(out) <- c("ED", "lam0", "sigma", "psi",
                       "N", "p", "deviance.ycap", "deviance.ydet")

    if(report>0)
        cat("\nstarting values =",
            round(c(ED, lam0, sigma, psi, N, p,
                    deviance.ycap, deviance.ydet), 3), "\n\n")

    for(iter in 1:n.iters) {

        if(iter %% report == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", round(out[iter-1,], 3), "\n")
            if(plotit) {
                plot(s, xlim=xlim, ylim=ylim, asp=1, cex=1,
                     pch=1, col=1)
                points(s[z==1,], pch=16, col="blue")
                points(x, pch="+")
                rect(xlim[1], ylim[1], xlim[2], ylim[2])
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

## if(iter==10)
##     browser()

        if(verbose)
            cat("   Sampling p(lam0,sigma|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
        
        ## Sample from p(lam0,sigma|.)
        lam0sigma.cand <- rmvnorm(1, c(lam0, sigma), proposal.vcov.lam0sigma)
        lam0.cand <- lam0sigma.cand[1]
        sigma.cand <- lam0sigma.cand[2]
        ld.ydet.cand[(n+1):M] <- 0
        if(lam0.cand > 0 & sigma.cand>0) {
            lambda.cand <- get_lambda_nou(lam0.cand, distSq, sigma.cand, far, notOper, nthreads) 
            ## ld.ydet.cand <- get_ld_y_nK(ydet, lambda.cand, nthreads)
            ## for(j in 1:J) {
            ##     for(k in 1:K) {
            ##         lambda.cand[,j,k] <- lam0.cand*exp(-distSq[,j]/(2*sigma.cand^2))*oper[j,k]##*!far[,j]
            ##     }
            ## }

            ld.ydet.cand[1:n] <- rowSums(dpois(ydet[1:n,,], lambda.cand[1:n,,], log=TRUE))
            ## ld.ydet.cand[(n+1):M] <- -rowSums(lambda.cand[(n+1):M,,])*z[(n+1):M]
            which.z1.aug <- which(z[(n+1):M]>0)+n
            ld.ydet.cand[which.z1.aug] <- -rowSums(lambda.cand[which.z1.aug,,])
            ld.ydet.cand.sum <- sum(ld.ydet.cand)
            ld.lam0 <- ld.lam0.cand <- 0
            ld.sigma <- ld.sigma.cand <- 0
            mhr.lam0sigma  <- exp((ld.ydet.cand.sum + ld.lam0.cand + ld.sigma.cand) -
                                  (ld.ydet.sum + ld.lam0 + ld.sigma))
            if(!is.finite(mhr.lam0sigma)) {
                browser()
                ## stop("mhr.lam0kappa = ", mhr.lam0kappa)
            }
            if(runif(1) < mhr.lam0sigma) {
                ld.ydet <- ld.ydet.cand
                ld.ydet.sum <- ld.ydet.cand.sum
                lambda <- lambda.cand
                lam0 <- lam0.cand
                sigma <- sigma.cand
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
        ld.ydet.cand <- ld.ydet
        ld.ydet.cand[(n+1):M] <- 0
        for(i in (n+1):M) {
            z.cand.i <- 1-z[i]
            ## Must compute ld.z because psi has changed
            ld.z[i] <- dbinom(z[i], 1, psi, log=TRUE)
            ld.z.cand[i] <- dbinom(z.cand.i, 1, psi, log=TRUE)
            if(z.cand.i==1) {
                ## ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s[i,1], s[i,2],
                ##                                        x, K, sigma, rho,
                ##                                        lam0, kappa,
                ##                                        trim^2, notOper, nthreads)
                ld.ydet.cand[i] <- -sum(lambda[i,,])  # dpois(0, sum(lambda[i,,]), log=TRUE)
                
            }
            ld.ycap[i] <- dbinom(ycap.aug[i], 1, p*z[i], log=TRUE)
            ld.ycap.cand[i] <- dbinom(ycap.aug[i], 1, p*z.cand.i, log=TRUE)
            mhr.z <- exp((ld.z.cand[i]+ld.ydet.cand[i]+ld.ycap.cand[i]) -
                         (ld.z[i]+ld.ydet[i]+ld.ycap[i]))
            if(!is.finite(mhr.z)) {
                browser()
                ## stop("mhr.z = ", mhr.z)
            }
            if(runif(1) < mhr.z) {
                z[i] <- z.cand.i
                ld.ydet[i] <- ld.ydet.cand[i]
                ld.ycap[i] <- ld.ycap.cand[i]
            }
        }
        N <- sum(z)

        ## Sample from p(p|.)
        if(verbose)
            cat("   Sampling p(p|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ## p.cand <- rnorm(1, p, tune[6])
        ## if(p.cand>0 & p.cand<1) {
        ##     ld.ycap <- dbinom(ycap.aug, 1, p*z, log=TRUE)
        ##     ld.ycap.cand <- dbinom(ycap.aug, 1, p.cand*z, log=TRUE)
        ##     mhr.p <- exp(sum(ld.ycap.cand) - sum(ld.ycap))
        ##     if(!is.finite(mhr.p)) {
        ##         browser()
        ##     }
        ##     if(runif(1) < mhr.p) {
        ##         p <- p.cand
        ##         ld.ycap <- ld.ycap.cand
        ##     }
        ## }

        ## Of N guys (not M), how many were detected. Assumes beta(1,1) prior
        p <- rbeta(1, 1+sum(ycap.aug), 1+N-sum(ycap.aug))
        
        ## Sample from p(s|.)
        if(verbose)
            cat("   Sampling p(s|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        s.ups <- 0
        ld.ydet.cand <- ld.ydet
        ld.ydet.cand[(n+1):M] <- 0
        distSq.cand <- distSq
        lambda.cand <- lambda
        far.cand <- far
        for(i in 1:M) {
            s1.cand <- rnorm(1, s[i,1], tune[3])
            s2.cand <- rnorm(1, s[i,2], tune[3])
            if(s1.cand < xlim[1] | s1.cand > xlim[2] | s2.cand < ylim[1] | s2.cand > ylim[2])
                next
            s.cand <- c(s1.cand, s2.cand)
            distSq.cand[i,] <- (s1.cand-x[,1])^2 + (s2.cand-x[,2])^2
            far.cand[i,] <- distSq.cand[i,]>(trim^2)
            for(k in 1:K) {
                lambda.cand[i,,k] <- lam0*exp(-distSq.cand[i,]/(2*sigma^2))*oper[,k]*(!far.cand[i,])
            }
            if(z[i]==1) {
                ## ld.ydet.aug.cand[i] <- get_ld_ydet_aug(n.mci, s1.cand, s2.cand, x, K, sigma,
                ##                                        rho, lam0, kappa, trim^2,
                ##                                        notOper, nthreads)
                if(i <= n) {
                    ld.ydet.cand[i] <- sum(dpois(ydet[i,,], lambda.cand[i,,], log=TRUE))
                } else {
                    ld.ydet.cand[i] <- -sum(lambda.cand[i,,])
                }
            }
            ## TODO: Add ld.s after implementing spatial variation in density
            mhr.s <-  exp(ld.ydet.cand[i] - ld.ydet[i])
            if(!is.finite(mhr.s))
                browser()
            if(runif(1) < mhr.s) {
                s[i,] <- s.cand
                ld.ydet[i] <- ld.ydet.cand[i]
                lambda[i,,] <- lambda.cand[i,,]
                distSq[i,] <- distSq.cand[i,]
                far[i,] <- far.cand[i,]
                s.ups <- s.ups+z[i]
            }
        }

        ld.ydet.sum <- sum(ld.ydet)

        deviance.ydet <- -2*ld.ydet.sum
        deviance.ycap <- -2*sum(ld.ycap)
        dev <- deviance.ydet + deviance.ycap

        if(iter %% report == 0) {
            cat("   Acceptance rates\n")
            cat("     s =", s.ups/N, "\n")
        }

        out[iter,] <- c(ED, lam0, sigma, psi, N,
                        p, deviance.ycap, deviance.ydet)

        if(keep.s.post) {
            if(iter %in% s.keep.seq) {
                s.post[,,s.counter] <- s
                ld.ydet.post[,s.counter] <- ld.ydet
                s.counter <- s.counter+1
            }
        }
    }

    return(list(samples=out,
                s.post=s.post,
                ld.ydet.post=ld.ydet.post,
                settings=list(tune=tune, trim=trim, ##nthreads=nthreads,
                              buffer=buffer, s.post.keep=s.post.keep),
                              ## block.rhosigma=block.rhosigma,
                              ## block.lam0sigma=block.lam0sigma),
                final.state=list(z=z, s=s, 
                                 theta=out[iter,],
                                 RNG.state=.Random.seed),
                data=list(ydet=ydet, ycap=ycap, M=M, x=x, oper=oper),
                seed=seed))
}





