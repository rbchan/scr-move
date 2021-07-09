library(coda)

list.files("fits055")

## Load MCMC output for SCR cases
load("fits055/fits-scr-055.gzip")
load("fits065/fits-scr-065.gzip")
load("fits075/fits-scr-075.gzip")
load("fits085/fits-scr-085.gzip")
load("fits095/fits-scr-095.gzip")

## Load MCMC output for SCR+move cases
load("fits055/fits-scr-move-055.gzip")
load("fits065/fits-scr-move-065.gzip")
load("fits075/fits-scr-move-075.gzip")
load("fits085/fits-scr-move-085.gzip")
load("fits095/fits-scr-move-095.gzip")

str(fits.scr.move.055[[1]])

## Summary stats

str(fits.scr.055[[1]]$data)

get.n.cap <- function(x) sum(x$data$ycap)

n.cap.055 <- sapply(fits.scr.move.055, get.n.cap)
n.cap.065 <- sapply(fits.scr.move.065, get.n.cap)
n.cap.075 <- sapply(fits.scr.move.075, get.n.cap)
n.cap.085 <- sapply(fits.scr.move.085, get.n.cap)
n.cap.095 <- sapply(fits.scr.move.095, get.n.cap)

summary(n.cap.055)
summary(n.cap.065)
summary(n.cap.075)
summary(n.cap.085)
summary(n.cap.095)


get.n.capdet <- function(x) sum(x$data$ycap==1 & rowSums(x$data$ydet)>0)

n.capdet.055 <- sapply(fits.scr.move.055, get.n.capdet)
n.capdet.065 <- sapply(fits.scr.move.065, get.n.capdet)
n.capdet.075 <- sapply(fits.scr.move.075, get.n.capdet)
n.capdet.085 <- sapply(fits.scr.move.085, get.n.capdet)
n.capdet.095 <- sapply(fits.scr.move.095, get.n.capdet)

summary(n.capdet.055)
summary(n.capdet.065)
summary(n.capdet.075)
summary(n.capdet.085)
summary(n.capdet.095)



## Extract posterior samples

tmp <- as.mcmc(fits.scr.055[[1]]$samples)
plot(tmp, ask=TRUE)

tmp <- as.mcmc(fits.scr.move.095[[1]]$samples)
plot(tmp, ask=TRUE)


get.samples <- function(x, burn=2000) x$samples[-(1:burn),]

samps.055 <- lapply(fits.scr.055, get.samples) 
samps.065 <- lapply(fits.scr.065, get.samples)
samps.075 <- lapply(fits.scr.075, get.samples)
samps.085 <- lapply(fits.scr.085, get.samples)
samps.095 <- lapply(fits.scr.095, get.samples)

samps.move.055 <- lapply(fits.scr.move.055, get.samples)
samps.move.065 <- lapply(fits.scr.move.065, get.samples)
samps.move.075 <- lapply(fits.scr.move.075, get.samples)
samps.move.085 <- lapply(fits.scr.move.085, get.samples)
samps.move.095 <- lapply(fits.scr.move.095, get.samples)


hist(sapply(samps.055, colMeans)["N",])



## Traceplots of N

samps.N.055 <- as.mcmc(sapply(samps.055, function(x) x[,"N"]))
samps.N.065 <- as.mcmc(sapply(samps.065, function(x) x[,"N"]))
samps.N.075 <- as.mcmc(sapply(samps.075, function(x) x[,"N"]))
samps.N.085 <- as.mcmc(sapply(samps.085, function(x) x[,"N"]))
samps.N.095 <- as.mcmc(sapply(samps.095, function(x) x[,"N"]))

samps.N.055.df <- data.frame(N=as.vector(unclass(samps.N.055)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.065.df <- data.frame(N=as.vector(unclass(samps.N.065)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.075.df <- data.frame(N=as.vector(unclass(samps.N.075)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.085.df <- data.frame(N=as.vector(unclass(samps.N.085)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.095.df <- data.frame(N=as.vector(unclass(samps.N.095)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))

samps.N.move.055 <- as.mcmc(sapply(samps.move.055, function(x) x[,"N"]))
samps.N.move.065 <- as.mcmc(sapply(samps.move.065, function(x) x[,"N"]))
samps.N.move.075 <- as.mcmc(sapply(samps.move.075, function(x) x[,"N"]))
samps.N.move.085 <- as.mcmc(sapply(samps.move.085, function(x) x[,"N"]))
samps.N.move.095 <- as.mcmc(sapply(samps.move.095, function(x) x[,"N"]))

samps.N.move.055.df <- data.frame(N=as.vector(unclass(samps.N.move.055)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.move.065.df <- data.frame(N=as.vector(unclass(samps.N.move.065)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.move.075.df <- data.frame(N=as.vector(unclass(samps.N.move.075)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.move.085.df <- data.frame(N=as.vector(unclass(samps.N.move.085)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))
samps.N.move.095.df <- data.frame(N=as.vector(unclass(samps.N.move.095)),
                             iteration=rep(1:10000, times=100),
                             sim=as.character(rep(1:100, each=10000)))


xyplot(N~iteration|sim, samps.N.055.df, type="l", strip=FALSE,
       layout=c(4,5), as.table=TRUE,
       panel=function(...) {
           panel.xyplot(...)
           lsegments(1, 100, 10000, 100)
       })



## xyplot(samps.N.055[,1:20], layout=c(5,4), as.table=TRUE,
##        ## scales=list(relation="same", alternating=3),
##        strip=FALSE,
##        panel=function(...) {
##            panel.xyplot(...)
##            lsegments(1, 100, 10000, 100)
##        })


png("trace-scr-055.png", width=9, height=20, res=100, units="in")
xyplot(N~iteration|sim, samps.N.055.df, type="l", strip=FALSE,
       layout=c(5,20), as.table=TRUE,
       panel=function(...) {
           panel.xyplot(...)
           lsegments(1, 100, 10000, 100)
       })
dev.off()
system("gopen trace-scr-055.png")


for(i in c("055", "065", "075", "085", "095")) {
    png(paste0("trace-scr-", i, ".png"), width=9, height=20,
        res=50, units="in")
    print(
    xyplot(N~iteration|sim, get(paste0("samps.N.", i, ".df")),
                                type="l", strip=FALSE,
           layout=c(5,20), as.table=TRUE,
           panel=function(...) {
               panel.xyplot(...)
               lsegments(1, 100, 10000, 100)
           })
    )
    dev.off()
}


for(i in c("055", "065", "075", "085", "095")) {
    png(paste0("trace-scr-move-", i, ".png"), width=9, height=20,
        res=50, units="in")
    print(
    xyplot(N~iteration|sim, get(paste0("samps.N.move.", i, ".df")),
                                type="l", strip=FALSE,
           layout=c(5,20), as.table=TRUE,
           panel=function(...) {
               panel.xyplot(...)
               lsegments(1, 100, 10000, 100)
           })
    )
    dev.off()
}





## Posterior means
get.pmean.N <- function(x) mean(x[,"N"])

pmean.N.055 <- sapply(samps.055, get.pmean.N)
pmean.N.065 <- sapply(samps.065, get.pmean.N)
pmean.N.075 <- sapply(samps.075, get.pmean.N)
pmean.N.085 <- sapply(samps.085, get.pmean.N)
pmean.N.095 <- sapply(samps.095, get.pmean.N)

pmean.N.move.055 <- sapply(samps.move.055, get.pmean.N)
pmean.N.move.065 <- sapply(samps.move.065, get.pmean.N)
pmean.N.move.075 <- sapply(samps.move.075, get.pmean.N)
pmean.N.move.085 <- sapply(samps.move.085, get.pmean.N)
pmean.N.move.095 <- sapply(samps.move.095, get.pmean.N)

## Posterior medians
get.pmedian.N <- function(x) median(x[,"N"])

pmedian.N.055 <- sapply(samps.055, get.pmedian.N)
pmedian.N.065 <- sapply(samps.065, get.pmedian.N)
pmedian.N.075 <- sapply(samps.075, get.pmedian.N)
pmedian.N.085 <- sapply(samps.085, get.pmedian.N)
pmedian.N.095 <- sapply(samps.095, get.pmedian.N)

pmedian.N.move.055 <- sapply(samps.move.055, get.pmedian.N)
pmedian.N.move.065 <- sapply(samps.move.065, get.pmedian.N)
pmedian.N.move.075 <- sapply(samps.move.075, get.pmedian.N)
pmedian.N.move.085 <- sapply(samps.move.085, get.pmedian.N)
pmedian.N.move.095 <- sapply(samps.move.095, get.pmedian.N)



## Posterior modes
get.pmode.N <- function(x, breakTies=TRUE) {
    N <- x[,"N"]
    uN <- unique(N)
    if(breakTies) {
        require(nnet)
        out <- uN[which.is.max(tabulate(match(N, uN)))]
    } else {
        out <- uN[which.max(tabulate(match(N, uN)))]
    }
    return(out)
}

pmode.N.055 <- sapply(samps.055, get.pmode.N)
pmode.N.065 <- sapply(samps.065, get.pmode.N)
pmode.N.075 <- sapply(samps.075, get.pmode.N)
pmode.N.085 <- sapply(samps.085, get.pmode.N)
pmode.N.095 <- sapply(samps.095, get.pmode.N)

pmode.N.move.055 <- sapply(samps.move.055, get.pmode.N)
pmode.N.move.065 <- sapply(samps.move.065, get.pmode.N)
pmode.N.move.075 <- sapply(samps.move.075, get.pmode.N)
pmode.N.move.085 <- sapply(samps.move.085, get.pmode.N)
pmode.N.move.095 <- sapply(samps.move.095, get.pmode.N)


## Posterior variance
get.pvar.N <- function(x) var(x[,"N"])

pvar.N.055 <- sapply(samps.055, get.pvar.N)
pvar.N.065 <- sapply(samps.065, get.pvar.N)
pvar.N.075 <- sapply(samps.075, get.pvar.N)
pvar.N.085 <- sapply(samps.085, get.pvar.N)
pvar.N.095 <- sapply(samps.095, get.pvar.N)

pvar.N.move.055 <- sapply(samps.move.055, get.pvar.N)
pvar.N.move.065 <- sapply(samps.move.065, get.pvar.N)
pvar.N.move.075 <- sapply(samps.move.075, get.pvar.N)
pvar.N.move.085 <- sapply(samps.move.085, get.pvar.N)
pvar.N.move.095 <- sapply(samps.move.095, get.pvar.N)

mean(pvar.N.055)
mean(pvar.N.065)
mean(pvar.N.075)
mean(pvar.N.085)
mean(pvar.N.095)

mean(pvar.N.move.055)
mean(pvar.N.move.065)
mean(pvar.N.move.075)
mean(pvar.N.move.085)
mean(pvar.N.move.095)



## CI coverage

cover <- function(x, a=100) {
    N <- x[,"N"]
    lower <- quantile(N, prob=0.025)
    upper <- quantile(N, prob=0.975)
    (lower <= a) & (upper >= a)
}

cover.N <- c(
    '055'=mean(cover.N.055 <- sapply(samps.055, cover)),
    '065'=mean(cover.N.065 <- sapply(samps.065, cover)),
    '075'=mean(cover.N.075 <- sapply(samps.075, cover)),
    '085'=mean(cover.N.085 <- sapply(samps.085, cover)),
    '095'=mean(cover.N.095 <- sapply(samps.095, cover))
)

cover.N.move <- c(
    '055'=mean(cover.N.move.055 <- sapply(samps.move.055, cover)),
    '065'=mean(cover.N.move.065 <- sapply(samps.move.065, cover)),
    '075'=mean(cover.N.move.075 <- sapply(samps.move.075, cover)),
    '085'=mean(cover.N.move.085 <- sapply(samps.move.085, cover)),
    '095'=mean(cover.N.move.095 <- sapply(samps.move.095, cover)))

covr <- cbind(cover.N, cover.N.move)

barplot(covr, beside=TRUE, ylim=c(0,1)); abline(h=0.95)
barplot(t(covr), beside=TRUE, ylim=c(0,1)); abline(h=0.95)

pdf("cover-N.pdf", width=7, height=7)
plot((1:5)-0.0, covr[,1], xlim=c(0.5, 5.5), ylim=c(0.85, 1),
     xaxt="n", ylab="CI coverage",
     xlab=expression(paste("Autocorrelation parameter (", rho, ")")))
axis(1, 1:5, paste(seq(0.55, 0.95, 0.1)))
points((1:5)+0.0, covr[,2], pch=16)
abline(h=0.95, lty=3)
arrows((1:5)-0.0, covr[,1], (1:5)+0.0, covr[,2], length=0.1, col="blue")
legend(0.5, 0.88, c("SCR", "SCR+move"), pch=c(1,16))
dev.off()
system("gopen cover-N.pdf")



## Bias

get.bias <- function(x, a) mean(x)-a

get.bias(pmean.N.055, a=100)
get.bias(pmean.N.065, a=100)
get.bias(pmean.N.075, a=100)
get.bias(pmean.N.085, a=100)
get.bias(pmean.N.095, a=100)

get.bias(pmean.N.move.055, a=100)
get.bias(pmean.N.move.065, a=100)
get.bias(pmean.N.move.075, a=100)
get.bias(pmean.N.move.085, a=100)
get.bias(pmean.N.move.095, a=100)

get.bias(pmedian.N.055, a=100)
get.bias(pmedian.N.065, a=100)
get.bias(pmedian.N.075, a=100)
get.bias(pmedian.N.085, a=100)
get.bias(pmedian.N.095, a=100)

get.bias(pmedian.N.move.055, a=100)
get.bias(pmedian.N.move.065, a=100)
get.bias(pmedian.N.move.075, a=100)
get.bias(pmedian.N.move.085, a=100)
get.bias(pmedian.N.move.095, a=100)

bias.N.mode <- c(
    '055'=get.bias(pmode.N.055, a=100),
    '065'=get.bias(pmode.N.065, a=100),
    '075'=get.bias(pmode.N.075, a=100),
    '085'=get.bias(pmode.N.085, a=100),
    '095'=get.bias(pmode.N.095, a=100))

bias.N.mode.move <- c(
    '055'=get.bias(pmode.N.move.055, a=100),
    '065'=get.bias(pmode.N.move.065, a=100),
    '075'=get.bias(pmode.N.move.075, a=100),
    '085'=get.bias(pmode.N.move.085, a=100),
    '095'=get.bias(pmode.N.move.095, a=100))

pdf("bias-N.pdf", width=7, height=7)
plot((1:5)-0.0, bias.N.mode, xlim=c(0.5, 5.5), ylim=c(-11, 2),
     xaxt="n", ylab="Bias",
     xlab=expression(paste("Autocorrelation parameter (", rho, ")")))
axis(1, 1:5, paste(seq(0.55, 0.95, 0.1)))
points((1:5)+0.0, bias.N.mode.move, pch=16)
abline(h=0, lty=3)
arrows((1:5)-0.0, bias.N.mode, (1:5)+0.0, bias.N.mode.move,
       length=0.1, col="blue")
legend(0.5, -8, c("SCR", "SCR+move"), pch=c(1,16))
dev.off()
system("gopen bias-N.pdf")

## Variance of point estimates

var(pmean.N.055)
var(pmean.N.065)
var(pmean.N.075)
var(pmean.N.085)
var(pmean.N.095)

var(pmean.N.move.055)
var(pmean.N.move.065)
var(pmean.N.move.075)
var(pmean.N.move.085)
var(pmean.N.move.095)

var(pmedian.N.055)
var(pmedian.N.065)
var(pmedian.N.075)
var(pmedian.N.085)
var(pmedian.N.095)

var(pmedian.N.move.055)
var(pmedian.N.move.065)
var(pmedian.N.move.075)
var(pmedian.N.move.085)
var(pmedian.N.move.095)

var.N.mode <- c(
    '055'=var(pmode.N.055),
    '065'=var(pmode.N.065),
    '075'=var(pmode.N.075),
    '085'=var(pmode.N.085),
    '095'=var(pmode.N.095))

var.N.mode.move <- c(
    '055'=var(pmode.N.move.055),
    '065'=var(pmode.N.move.065),
    '075'=var(pmode.N.move.075),
    '085'=var(pmode.N.move.085),
    '095'=var(pmode.N.move.095))

pdf("var-N.pdf", width=7, height=7)
plot((1:5)-0.0, var.N.mode, xlim=c(0.5, 5.5), ylim=c(300, 650),
     xaxt="n", ylab="Sampling variance",
     xlab=expression(paste("Autocorrelation parameter (", rho, ")")))
axis(1, 1:5, paste(seq(0.55, 0.95, 0.1)))
points((1:5)+0.0, var.N.mode.move, pch=16)
abline(h=0, lty=3)
arrows((1:5)-0.0, var.N.mode, (1:5)+0.0, var.N.mode.move,
       length=0.1, col="blue")
legend(0.5, 650, c("SCR", "SCR+move"), pch=c(1,16))
dev.off()
system("gopen var-N.pdf")


## RMSE

get.rmse <- function(x, a) sqrt(mean((x-a)^2))

get.rmse(pmean.N.055, a=100)
get.rmse(pmean.N.065, a=100)
get.rmse(pmean.N.075, a=100)
get.rmse(pmean.N.085, a=100)
get.rmse(pmean.N.095, a=100)

get.rmse(pmean.N.move.055, a=100)
get.rmse(pmean.N.move.065, a=100)
get.rmse(pmean.N.move.075, a=100)
get.rmse(pmean.N.move.085, a=100)
get.rmse(pmean.N.move.095, a=100)

get.rmse(pmedian.N.055, a=100)
get.rmse(pmedian.N.065, a=100)
get.rmse(pmedian.N.075, a=100)
get.rmse(pmedian.N.085, a=100)
get.rmse(pmedian.N.095, a=100)

get.rmse(pmedian.N.move.055, a=100)
get.rmse(pmedian.N.move.065, a=100)
get.rmse(pmedian.N.move.075, a=100)
get.rmse(pmedian.N.move.085, a=100)
get.rmse(pmedian.N.move.095, a=100)

get.rmse(pmode.N.055, a=100)
get.rmse(pmode.N.065, a=100)
get.rmse(pmode.N.075, a=100)
get.rmse(pmode.N.085, a=100)
get.rmse(pmode.N.095, a=100)

get.rmse(pmode.N.move.055, a=100)
get.rmse(pmode.N.move.065, a=100)
get.rmse(pmode.N.move.075, a=100)
get.rmse(pmode.N.move.085, a=100)
get.rmse(pmode.N.move.095, a=100)




## Impact of n-capdet on mean(N) and var(N)

plot(n.capdet.055, pmean.N.move.055)
plot(n.capdet.065, pmean.N.move.065)
plot(n.capdet.075, pmean.N.move.075)
plot(n.capdet.085, pmean.N.move.085)
plot(n.capdet.095, pmean.N.move.095)

plot(n.capdet.055, pvar.N.move.055)
plot(n.capdet.065, pvar.N.move.065)
plot(n.capdet.075, pvar.N.move.075)
plot(n.capdet.085, pvar.N.move.085)
plot(n.capdet.095, pvar.N.move.095)


get.bias.paired <- function(x1, x2, a) mean((x1-x2))





