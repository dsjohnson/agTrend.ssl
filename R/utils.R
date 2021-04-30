################################################################################
### Some package imports
################################################################################

#'@importFrom stats coef
#'@importFrom stats plogis
#'@importFrom stats vcov
#'@importFrom stats median
#'@importFrom stats terms
#'@rawNamespace import(mgcv, except = rmvn)
#'@importFrom mvnfast rmvn
#'@importFrom rlang .data
#'@importFrom utils globalVariables
#'@importFrom stringi stri_trans_general
#'@import dplyr purrr tidyr readxl coda

utils::globalVariables(c(".", ".x","SITE","site","year","X","X.mu","X.p","X.phi","mu","p","phi","REGION","obl","RCA"))


################################################################################
### Some utility functions
################################################################################

draw.tw <- function(mu, p, phi, ...){
  nc <- ifelse(is.null(ncol(mu)), 1, ncol(mu))
  nr <- ifelse(is.null(nrow(mu)), length(mu), nrow(mu))
  o <- sapply(1:nr,
              function(i, mu, p, phi){mgcv::rTweedie(mu[i,],p[i],phi[i])},
              mu=mu, p=p, phi=phi
  )
  return(t(o))
}

get.real <- function(N.pred, data, ...){
  surv.times <- which(!is.na(data$count))
  data <- data[surv.times,]
  out <- N.pred
  for(i in 1:nrow(data)){
    if(attr(data,'obl.corr')){
    out[,surv.times[i]] <- data$count[i]*exp(0.03903366*data$obl[i])
    } else{
      out[,surv.times[i]] <- data$count[i]
    }
  }
  return(out)
}

