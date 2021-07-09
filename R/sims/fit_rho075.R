
rho <- 0.75
stub <- "075"

data.file <- paste0("sims", stub, ".gzip")
data.obj <- paste0("sims", rho)
fits.scr.move.name <- paste0("fits.scr.move.", stub)
fits.scr.name <- paste0("fits.scr.", stub)
fits.scr.move.out <- paste0("fits-scr-move-", stub)
fits.scr.out <- paste0("fits-scr-", stub)
out.dir <- paste0("fits", stub, "/")

load(data.file)

n.sims <- length(get(data.obj))

print(.libPaths())
.libPaths("./")      ## Add local library
print(.libPaths())
install.packages("../../scrmove", repos=NULL, type="source", lib=.libPaths()[1])

library(scrmove)


fits.scr.move <- vector(mode="list", length=n.sims)
fits.scr <- vector(mode="list", length=n.sims)



set.seed(43230)

for(i in 1:n.sims) {

    cat("Analyzing simulated dataset", i, ",", format(Sys.time()), "\n")

    tmp <- get(data.obj)[[i]]

    fits.scr.move[[i]] <-
        scrMoveDA(ycap=as.matrix(tmp$y.cap),
                  ydet=tmp$y.det,
                  u=tmp$u,
                  M=225,
                  x=tmp$x,
                  ## oper=oper,
                  plotit=FALSE, 
                  n.iters=12000, n.mci=100L,
                  buffer=5000, trim=1000, nthreads=4,
                  report=1000, verbose=FALSE,
                  block.rhosigma=TRUE,
                  block.lam0kappa=TRUE,
                  ## inits=ini,
                  ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
                  tune=c(0.29, 3e-3, 15, 0.29, 4, 0.1, 300, 300, 3.01e-3, -0.05))
    
    fits.scr[[i]] <-
        scrDA(ycap=as.matrix(tmp$y.cap),
              ydet=tmp$y.det,
              M=225,
              x=tmp$x,
              ## oper=oper,
              plotit=FALSE, 
              n.iters=12000, ##n.mci=100L,
              buffer=5000, trim=Inf, ##nthreads=4,
              report=1000, verbose=FALSE,
              ## block.rhosigma=TRUE,
              ## block.lam0kappa=TRUE,
              ## inits=ini,
              ## tune=c(0.29, 0.0007, 4, 0.15, 1.8, 0.11, 300, 300))
              tune=c(0.01, 20, 300, -0.001))

}



assign(fits.scr.move.name, fits.scr.move)
assign(fits.scr.name, fits.scr)


dir.create(out.dir)


save(list=fits.scr.move.name, file=paste0(out.dir, fits.scr.move.out, ".gzip"))
save(list=fits.scr.name, file=paste0(out.dir, fits.scr.out, ".gzip"))
