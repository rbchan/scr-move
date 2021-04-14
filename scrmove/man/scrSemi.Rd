\name{scrSemi}
\alias{scrSemi}

\title{
  MCMC for SCR0 model. Conditional-on-n (semi-likelihood) approach.
}
\description{
  Basic SCR model assuming encounter
  histories are Poisson, and the activity centers follow a homogeneous
  Poisson point process. 
}
\usage{
scrSemi(ycap, ydet, x, oper, aggregate.ydet = FALSE, report = 10,
  verbose = FALSE, plotit = TRUE, tune, n.iters, n.mci, buffer, trim =
  600, nthreads, inits) 
}
\arguments{
  \item{ycap}{
    nIndividuals x 1 binary matrix indicating if each individual was
    captured and collared (before camera data)
}
  \item{ydet}{
    nIndividuals x nTraps x nOccasions array of encounter
    frequencies. Meant to be the camera trap data.
}
  \item{x}{
   nTraps x 2 matrix of camera trap locaions
}
  \item{oper}{
   nTraps x nOccaions binary matrix indicating when cameras were operational
}
  \item{aggregate.ydet}{
    Logical. Should encounter histories be collapsed over time to speed
  up inference?
}
  \item{report}{
    Positive integer for progress report
}
  \item{verbose}{
    Logical. Extra information about speed.
}
  \item{plotit}{
    Visualize the posterior sample of s and u over x in S
}
  \item{tune}{
    MH tuning parameters
}
  \item{n.iters}{
    MCMC iterations
}
  \item{n.mci}{
    nSimulations for Monte Carlo integration
}
  \item{buffer}{
    Used to define rectangular state-space
}
  \item{trim}{
    Used to speed up computation by assuming encounter rate is 0 beyond
    this distance
}
  \item{nthreads}{
    Number of threads for parallel computing with OpenMP
}
  \item{inits}{
    Initial values (from a previous run)
}
}
\details{
  This function uses the conditional-on-n (semi-likelihood) approach to
  capture-recapture. Monte Carlo integration is used to compute the
  probability of not being detected.
}
\value{
A list with the following components
  \item{samples}{MCMC samples for high level parameters}
  \item{s.post}{MCMC samples for latent s variables}
  \item{u.post}{MCMC samples for latent u variables}
  \item{final.state}{State of MCMC chain at final iteration. Can be used
    as \code{inits} to continue sampling} 
  \item{settings}{Arguments}
}
\references{
  King et al. AOAS
  Chandler et al. 2021? Ecology
}
\note{
  Meant to be fitted to camera trap data and telemetry data from a study
  of deer in South Florida. Created for manuscript to appear in special
  issue in Ecology. Work in progress. 
}

\section{Warning }{
  This is not meant to be a general-purpose function. It was designed
  for the dataset \code{\link{deer_scr_telem}}.
}


\seealso{
  \code{\link{scrMoveSemi}}, ~~~
}
\examples{
ls()
}
\keyword{SCR}

