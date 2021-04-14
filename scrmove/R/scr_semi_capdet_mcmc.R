## MCMC algorithm for SCR0 (no movement or telemetry data)
## No data augmentation. Uses "semi-complete likelihood" appraoch instead



scrSemi <- function(ycap,  ## nx1 initial (helicopter-based) capture histories
                    ydet,  ## nxJxK camera trap encounter histories
                    x,     ## Trap coordinates
                    oper,  
                    aggregate.ydet=FALSE,
                    report=10, verbose=FALSE, plotit=TRUE,
                    tune,
                    n.iters, n.mci, buffer, trim=600, nthreads,
                    s.post.keep=0, do.modsel=FALSE,
                    inits) {

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    ydet.dims <- dim(ydet)
    n <- ydet.dims[1]
    J <- ydet.dims[2]
    K <- ydet.dims[3]

    if(nrow(ycap) != nrow(ydet))
        stop("nrow(ycap) != nrow(ydet)")
    if(J != nrow(x))
        stop("ydet and x should have the same number of traps")

    require(parallel)
    availableCores <- detectCores()
    if(missing(nthreads))
        nthreads <- availableCores-1

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
    trapOper <- rowSums(oper)

    if(aggregate.ydet) {
        ydet.agg <- apply(ydet, c(1, 2), sum)
    }

    ## Initial values
    if(missing(inits)) {
        s <- cbind(runif(n, xlim[1], xlim[2]),
                   runif(n, ylim[1], ylim[2]))
        for(i in 1:n) {
            traps.w.dets.i <- rowSums(ydet[i,,])>0
            if(sum(traps.w.dets.i)>0)
                s[i,] <- colMeans(x[traps.w.dets.i,,drop=FALSE])
        }
        beta0 <- runif(1, log(0.3), log(0.4))
        ## beta0 <- runif(1, log(M/2/Area), log(M/1.5/Area)) ## log(E(Density))
        beta1 <- 0
        N <- n+ceiling(runif(1, 0, 100))
        lam0 <- runif(1, 0.001, 0.005)
        kappa <- runif(1, 800, 1000)
        p <- runif(1)
    } else {
        beta0 <- log(inits$theta["ED"])
        lam0 <- inits$theta["lam0"]
        kappa <- inits$theta["kappa"]
        N <- inits$theta["N"]
        p <- inits$theta["p"]
        s <- inits$s
        .Random.seed <- inits$RNG.state
    }

    ## Note: It would be much faster to make lambda a n x J matrix
    ##       but I'm being lazy and adapting this from scr_move_semi_capdet_mcmc
    
    ## p(y|u,sigma,lam0,varsigma)
    if(aggregate.ydet) {
        far <- lambda <- distSq <- array(NA_real_, c(n, J))
        for(j in 1:J) {
            distSq[,j] <- (s[,1]-x[j,1])^2 + (s[,2]-x[j,2])^2
            far[,j] <- distSq[,j]>(trim^2)
            lambda[,j] <- ifelse(far[,j], 0,
                                 lam0*exp(-distSq[,j]/(2*kappa^2))*trapOper[j])
        }
    } else {
        far <- lambda <- distSq <- array(NA_real_, c(n, J, K))
        for(j in 1:J) {
            for(k in 1:K) {
                distSq[,j,k] <- (s[,1]-x[j,1])^2 + (s[,2]-x[j,2])^2
                far[,j,k] <- distSq[,j,k]>(trim^2)
                lambda[,j,k] <- ifelse(far[,j,k], 0,
                                       lam0*exp(-distSq[,j,k]/(2*kappa^2))*oper[j,k])
            }
        }
    }

    ## browser()

    ## Functions to compute log-densities
    
    get.ld.N <- function(N, EN)
        return(dpois(N, EN, log=TRUE))      ## For Poisson point process
    get.ld.n0 <- function(n, N, q.star) {       ## Likelihood of missing n0 guys.
        n0 <- N-n
        if(n0<0) {
            return(-Inf)
        } else if(q.star<0 | q.star>1 | !is.finite(q.star)) {
            warning("q.star should be a probability")
            return(-Inf)
        } else {
            return(lgamma(N+1)-lgamma(n0+1)-lgamma(n+1)+log(q.star)*n0)
        }
    }

    if(aggregate.ydet) {
        ld.ydet <- ld.ydet.cand <- rowSums(dpois(ydet.agg, lambda, TRUE))
    } else {
        ld.ydet <- ld.ydet.cand <- get_ld_y_noK(ydet, lambda, nthreads)
    }
    ld.ydet.sum <- ld.ydet.cand.sum <- sum(ld.ydet)

    ld.ycap <- dbinom(ycap, 1, p, log=TRUE)
    ld.ycap.sum <- sum(ld.ycap)

    ## p(N|EN(beta0))
    logpixArea <- 0 ## Replace with pixel area if raster covariates are added
    logmu <- beta0 + logpixArea ## Could add covariates here
    ##    EN <- sum(exp(logmu))
    ED <- exp(logmu)
    EN <- ED*Area
    logEN <- log(EN)
    ld.N <- get.ld.N(N, EN)

    ld.kappa <- 0 #dgamma(kappa, 1, 1, log=TRUE)

    ## p(n0|N,q.star)
    if(n.mci %% 2 !=0)
        stop("n.mci must be even")
    q.star.cam <- get_qstar_cam_nomove(n.mci, x, K, lam0, kappa,
                                       buffer, (trim)^2, trapOper, nthreads)
    q.star.gps <- 1-p
    q.star <- q.star.gps*q.star.cam
    ld.n0 <- get.ld.n0(n, N, q.star)

    ##browser()

    ## p(s|beta)
    ## NOTE: Assuming homogeneous point process for now: s(i)~Unif(s)
    ## ld.s <- ld.s.cand <- logmu - logEN
    ## ld.s.sum <- ld.s.cand.sum <- sum(ld.s)
    ld.s <- ld.s.cand <- 0

    deviance.ydet <- -2*ld.ydet.sum 
    deviance.ycap <- -2*ld.ycap.sum 

    ordr <- sample.int(n)
    icols1 <- topo.colors(n)[ordr]
    icols2 <- topo.colors(n, alpha=0.3)[ordr]

    keep.s.post <- s.post.keep>0
    if(keep.s.post) {
        s.keep.seq <- seq(1, n.iters, s.post.keep)
        niter.s.keep <- length(s.keep.seq)
        s.post <- array(NA_real_, c(n, 2, niter.s.keep))
        s.counter <- 1
    } else {
        s.post <- NULL
    }
    
    if(do.modsel) {
        ld.var.ycap <- list(mean=ld.ycap, var=ld.ycap)
        ld.var.ycap$mean[] <- ld.var.ycap$var[] <- 0
        ld.var.ydet <- list(mean=ld.ydet, var=ld.ydet)
        ld.var.ydet$mean[] <- ld.var.ydet$var[] <- 0
        ld.var.N <- list(mean=ld.N, var=ld.N)
        ld.var.N$mean[] <- ld.var.N$var[] <- 0
        ld.var.n0 <- list(mean=ld.n0, var=ld.n0)
        ld.var.n0$mean[] <- ld.var.n0$var[] <- 0
        d.bar.ycap <- ld.ycap
        d.bar.ycap[] <- 0
        d.bar.ydet <- ld.ydet
        d.bar.ydet[] <- 0
        d.bar.N <- 0
        d.bar.n0 <- 0
        running.var <- function(x.mean, x.var, x.new, it, nit) {
            x.mean.old <- x.mean
            x.mean.new <- x.mean+(x.new-x.mean)/it
            x.var.new <- x.var+(x.new-x.mean.old)*(x.new-x.mean.new)/(nit-1)
            return(list(mean=x.mean.new, var=x.var.new))
        }
    } 

    out <- matrix(NA, nrow=n.iters, ncol=8)
    colnames(out) <- c("ED", "lam0", "kappa", "q.star",
                       "N", "p", "deviance.ycap", "deviance.ydet")

    if(report>0)
        cat("\nstarting values =",
            round(c(ED, lam0, kappa, q.star, N, p,
                    deviance.ycap, deviance.ydet), 3), "\n\n")

    for(iter in 1:n.iters) {

        if(iter %% report == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", round(out[iter-1,], 3), "\n")
            if(plotit) {
                plot(s, xlim=xlim, ylim=ylim, asp=1, cex=2,
                     pch=16, col=icols1)
                points(x, pch="+")
                rect(xlim[1], ylim[1], xlim[2], ylim[2])
            }
        }

        if(verbose)
            cat("   Sampling p(beta0|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ## Sample from p(beta0|.)
        beta0.cand <- rnorm(1, beta0, tune[1])
        logmu.cand <- beta0.cand + logpixArea
        ED.cand <- exp(logmu.cand)
        EN.cand <- ED.cand*Area
        ld.N.cand <- get.ld.N(N, EN.cand)
        ld.beta0.cand <- ld.beta0 <- 0
        if(runif(1) < exp((ld.N.cand + ld.beta0.cand) -
                          (ld.N + ld.beta0))) {
            ## ld.s <- ld.s.cand
            ## ld.s.sum <- ld.s.cand.sum
            ld.N <- ld.N.cand
            beta0 <- beta0.cand
            ED <- ED.cand
            EN <- EN.cand
            logEN <- log(EN)
        }


        if(verbose)
            cat("   Sampling p(lam0|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ## Sample from p(lam0|.)
        lam0.cand <- rnorm(1, lam0, tune[2])
        if(lam0.cand > 0) {
            if(aggregate.ydet) {
                lambda.cand <- get_lambda_noK(lam0.cand, distSq, kappa, far, trapOper, nthreads)
                ld.ydet.cand <- rowSums(dpois(ydet.agg, lambda.cand, log=TRUE))
            } else {
                lambda.cand <- get_lambda(lam0.cand, distSq, kappa, far, notOper, nthreads)
                ld.ydet.cand <- get_ld_y_noK(ydet, lambda.cand, nthreads)
            }
            ld.ydet.cand.sum <- sum(ld.ydet.cand)
            q.star.cam.cand <- get_qstar_cam_nomove(n.mci, x, K, 
                                                    lam0.cand, kappa, buffer, (trim)^2,
                                                    trapOper, nthreads)
            q.star.cand <- (1-p)*q.star.cam.cand
            ld.n0.cand <- get.ld.n0(n, N, q.star.cand)
            ld.lam0 <- ld.lam0.cand <- 0
            if(runif(1) < exp((ld.ydet.cand.sum + ld.n0.cand + ld.lam0.cand) -
                              (ld.ydet.sum + ld.n0 + ld.lam0))) {
                ld.ydet <- ld.ydet.cand
                ld.ydet.sum <- ld.ydet.cand.sum
                q.star.cam <- q.star.cam.cand
                q.star <- q.star.cand
                ld.n0 <- ld.n0.cand
                lambda <- lambda.cand
                lam0 <- lam0.cand
            }
        }

        if(verbose)
            cat("   Sampling p(kappa|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ## Sample from p(kappa|.)
        kappa.cand <- rnorm(1, kappa, tune[3]) ##sqrt(1/(2*lam1.cand))
        if(kappa.cand>0) { 
            ## lambda.cand <- get_lambda(lam0, distSq, kappa.cand, far, notOper, nthreads) 
            ## ld.ydet.cand <- get_ld_y_noK(ydet, lambda.cand, nthreads)
            if(aggregate.ydet) {
                lambda.cand <- get_lambda_noK(lam0, distSq, kappa.cand, far, trapOper, nthreads)
                ld.ydet.cand <- rowSums(dpois(ydet.agg, lambda.cand, log=TRUE))
            } else {
                lambda.cand <- get_lambda(lam0, distSq, kappa.cand, far, notOper, nthreads)
                ld.ydet.cand <- get_ld_y_noK(ydet, lambda.cand, nthreads)
            }
            ld.ydet.cand.sum <- sum(ld.ydet.cand)
            q.star.cam.cand <- get_qstar_cam_nomove(n.mci, x, K, lam0,
                                                    kappa.cand, buffer, (trim)^2,
                                                    trapOper, nthreads)
            q.star.cand <- (1-p)*q.star.cam.cand
            ld.n0.cand <- get.ld.n0(n, N, q.star.cand)
            ld.kappa.cand <- 0 ##dgamma(kappa.cand, 1, 1, log=TRUE)
            if(runif(1) < exp((ld.ydet.cand.sum + ld.n0.cand + ld.kappa.cand) -
                              (ld.ydet.sum + ld.n0 + ld.kappa))) {
                ld.ydet <- ld.ydet.cand
                ld.ydet.sum <- ld.ydet.cand.sum
                q.star.cam <- q.star.cam.cand
                q.star <- q.star.cand
                ld.n0 <- ld.n0.cand
                lambda <- lambda.cand
                kappa <- kappa.cand
                ld.kappa <- ld.kappa.cand
            }
        }


        if(verbose)
            cat("   Sampling p(N|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        ## Sample from p(N|.)
        N.cand <- rpois(1, N)
        ld.N.cand <- get.ld.N(N.cand, EN)
        ld.n0.cand <- get.ld.n0(n, N.cand, q.star)
        ld.N.N.cand <- dpois(N, N.cand, log=TRUE)
        ld.N.cand.N <- dpois(N.cand, N, log=TRUE)
        if(runif(1) < exp((ld.n0.cand + ld.N.cand + ld.N.N.cand) -
                          (ld.n0 + ld.N + ld.N.cand.N))) {
            ld.N <- ld.N.cand
            ld.n0 <- ld.n0.cand
            N <- N.cand
        }

        ## Sample from p(p|.)
        if(verbose)
            cat("   Sampling p(p|.): ", format(Sys.time(), "%H:%M:%S"), "\n")
        p.cand <- rnorm(1, p, tune[4])
        if(p.cand>0 & p.cand<1) {
            ld.ycap.cand <- dbinom(ycap, 1, p.cand, log=TRUE)
            q.star.cand <- (1-p.cand)*q.star.cam
            ld.n0.cand <- get.ld.n0(n, N, q.star.cand)
            if(runif(1) < exp((sum(ld.ycap.cand)+ld.n0.cand) -
                              (sum(ld.ycap)+ld.n0))) {
                p <- p.cand
                ld.ycap <- ld.ycap.cand
                q.star <- q.star.cand
                ld.n0 <- ld.n0.cand
            }
        }


        ## Sample from p(s|.)
        if(verbose)
            cat("   Sampling p(s|.): ", format(Sys.time(), "%H:%M:%S"), "\n")

        sups <- 0
        s.cand <- cbind(rnorm(n, s[,1], tune[5]), rnorm(n, s[,2], tune[5]))
        skip <- (s.cand[,1] < xlim[1]) | (s.cand[,1] > xlim[2]) |
                (s.cand[,2] < ylim[1]) | (s.cand[,2] > ylim[2]) 
        ##        lambda.cand <- lambda
        if(aggregate.ydet) {
            far.cand <- lambda.cand <- distSq.cand <- array(NA_real_, c(n, J))
            for(j in 1:J) {
                distSq.cand[,j] <- (s.cand[,1]-x[j,1])^2 + (s.cand[,2]-x[j,2])^2
                far.cand[,j] <- distSq.cand[,j]>(trim^2)
                lambda.cand[,j] <- ifelse(far.cand[,j], 0,
                                          lam0*exp(-distSq.cand[,j]/(2*kappa^2))*trapOper[j])
            }
            ld.ydet.cand <- rowSums(dpois(ydet.agg, lambda.cand, log=TRUE))
        } else {
            far.cand <- lambda.cand <- distSq.cand <- array(NA_real_, c(n, J, K))
            for(j in 1:J) {
                for(k in 1:K) {
                    ## Stupid and slow, being lazy
                    distSq.cand[,j,k] <- (s.cand[,1]-x[j,1])^2 + (s.cand[,2]-x[j,2])^2
                    far.cand[,j,k] <- distSq.cand[,j,k]>(trim^2)
                    lambda.cand[,j,k] <- ifelse(far.cand[,j,k], 0,
                                                lam0*exp(-distSq.cand[,j,k]/(2*kappa^2))*oper[j,k])
                }
            }
            ld.ydet.cand <- get_ld_y_noK(ydet, lambda.cand, nthreads)
        }
        accept.s <- (runif(n) < exp(ld.ydet.cand - ld.ydet)) & (!skip)
        sups <- sum(accept.s)
        if(sups>0) {
            s[accept.s,] <- s.cand[accept.s,]
            ld.ydet[accept.s] <- ld.ydet.cand[accept.s]
            if(aggregate.ydet) {
                distSq[accept.s,] <- distSq.cand[accept.s,]
                far[accept.s,] <- far.cand[accept.s,]
                lambda[accept.s,] <- lambda.cand[accept.s,]
            } else {
                distSq[accept.s,,] <- distSq.cand[accept.s,,]
                far[accept.s,,] <- far.cand[accept.s,,]
                lambda[accept.s,,] <- lambda.cand[accept.s,,]
            }
        }

        ld.ydet.sum <- sum(ld.ydet)

        deviance.ydet <- -2*ld.ydet.sum 
        deviance.ycap <- -2*sum(ld.ycap)

        if(iter %% report == 0) {
            cat("   Acceptance rates\n")
            cat("     s =", sups/n, "\n")
        }
        out[iter,] <- c(ED, lam0, kappa, q.star, N,
                        p, deviance.ycap, deviance.ydet)

        if(keep.s.post) {
            if(iter %in% s.keep.seq) {
                s.post[,,s.counter] <- s
                s.counter <- s.counter+1
            }
        }

        ## Stuff to be used for computing WAIC
        if(do.modsel) {
            d.bar.ycap <- d.bar.ycap+(exp(ld.ycap)-d.bar.ycap)/iter
            ## Does this need to be trap-level?
            d.bar.ydet <- d.bar.ydet+(exp(ld.ydet)-d.bar.ydet)/iter
            d.bar.N <- d.bar.N+(exp(ld.N)-d.bar.N)/iter
            d.bar.n0 <- d.bar.n0+(exp(ld.n0)-d.bar.n0)/iter
            ld.var.ycap <- running.var(ld.var.ycap[[1]], ld.var.ycap[[2]], ld.ycap, iter, n.iters)
            ld.var.ydet <- running.var(ld.var.ydet[[1]], ld.var.ydet[[2]], ld.ydet, iter, n.iters)
            ld.var.N <- running.var(ld.var.N[[1]], ld.var.N[[2]], ld.N, iter, n.iters)
            ld.var.n0 <- running.var(ld.var.n0[[1]], ld.var.n0[[2]], ld.n0, iter, n.iters)
        }

    }

    if(do.modsel) {
        modsel <- list()
        modsel$lppd.ycap <- log(d.bar.ycap)
        modsel$lppd.ydet <- log(d.bar.ydet)
        modsel$lppd.N <- log(d.bar.N)
        modsel$lppd.n0 <- log(d.bar.n0)
        modsel$pen.ycap <- ld.var.ycap
        modsel$pen.ydet <- ld.var.ydet
        modsel$pen.N <- ld.var.N
        modsel$pen.n0 <- ld.var.n0
    } else {
        modsel <- NULL
    }

    return(list(samples=out,
                s.post=s.post,
                settings=list(tune=tune, trim=trim, nthreads=nthreads,
                              buffer=buffer, s.post.keep=s.post.keep),
                final.state=list(s=s, 
                                 theta=out[iter,],
                                 RNG.state=.Random.seed),
                modsel=modsel,
                data=list(ydet=ydet, ycap=ycap, x=x, oper=oper),
                seed=seed))
}




