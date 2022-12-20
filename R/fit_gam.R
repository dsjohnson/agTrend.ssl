################################################################################
### Fit agTrend-like model to data
################################################################################
#' @title Fit penalized GAM to each site for missing data imputation
#' @description This function takes a processed data frame with counts during survey
#' years for each site and fits a multi-site GAM using a Tweedie response distribution.
#' @param data A data frame containing the sites, counts, and years.
#' @param obl.corr Logical. Should oblique photo correction be used?
#' @param alt.mod An alternative to the defaul model. See Details.
#' @param warn Logical. Should fitting warnings be printed. Defaults to \code{FALSE}.
#' @param debug Logical. If set to \code{TRUE} the function will drop into browser
#' mode within the function upon execution. Defaults to \code{FALSE}.
#' @param ... Extra arguments passed to \code{\link[mgcv]{gam}}
#'@details A penalized GAM model is fit to the data with a Tweedie distribution for the response distribution.
#'The mean function for site i in year t is modeled with a smooth year term for all sites plus
#'a factor smooth such that each site has its own smooth year function as well the \code{mgcv} formula is
#'\code{mu.form = count ~ s(year, SITE, bs='fs', k=8, m=1)}. For the power and dispersion formulas,
#'\code{p.form = ~ s(SITE, bs='re')} and \code{phi.form = ~ s(SITE, bs='re')}, so that each
#'site has its own p and phi parameters for the Tweedie specification. Thus, for each site, the
#'mean count is mu(i,t) and the variance is V(i,t) = phi(i)*mu(i,t)^p(i). To specify a different formula for the model, set
#'\code{alt.model = list(mu=mu.form.alt, p=p.form.alt, phi=phi.form.alt)} where the \code{*.alt} signifies
#'the formula for the desired alternative.
#'
#'If \code{obl.corr = TRUE} then a model will be fitted that corrects for the approximately
#'3.8% reduction in the expected count if the data were collected using an oblique
#'photo vs. a medium format vertical photo. Unlike the original \code{agTrend} package
#'the uncertainty in this estimate is not accounted for in the model. A few initial tests
#'revealed that is source of variation seems insignificant when compared to the
#'natural variation of the observed counts and the model was significantly easier to
#'fit and more robust when fixing this quantity.
#'
#'@author Devin S. Johnson
#'@import dplyr
#'@export

fit.gam <- function(data,
                    obl.corr,
                    alt.mod = NULL,
                    warn = FALSE,
                    debug = FALSE,...){

  if(debug) browser()
  data <- dplyr::arrange(data, .data$site, .data$year)
  data$site <- factor(data$site)

  ### Model definition
  if(!is.null(alt.mod)){
    ck.nms <- all(names(alt.mod)%in%c('mu.form','p.form','phi.form'))
    mod <- list(mu.form=NULL, p.form=NULL, phi.form=NULL)
    if(!is.null(alt.mod$mu.form)){
      mod$mu.form <- alt.mod$mu.form
    } else{
      mod$mu.form <- count ~  s(year, site, bs="fs", m=1)
    }
    if(!is.null(alt.mod$p.form)){
      mod$p.form <- alt.mod$p.form
    } else{
      mod$p.form <- ~s(site,bs='re')
    }
    if(!is.null(alt.mod$phi.form)){
      mod$phi.form <- alt.mod$phi.form
    } else{
      mod$phi.form <- ~s(site,bs='re')
    }

  } else{
    mod <- list(
      mu.form = count ~  s(year, site, bs="fs", m=1),
      p.form = ~s(site,bs='re'),
      phi.form = ~s(site,bs='re')
    )
  }

  mod <- list(mod$mu.form, mod$p.form, mod$phi.form)

  if(obl.corr){
    if(! 'obl' %in% colnames(data)) stop("If using oblique photo correction 'obl' must be a column in the data!")
    if(!warn){
      fit.gam <- suppressWarnings(
        gam(mod, offset = -0.03903366*obl, data=data, family=twlss(), select=TRUE, ...)
      )
    } else{
      fit.gam <- gam(mod, offset = -0.03903366*obl, data=data, family=twlss(), select=TRUE, ...)
    }

  } else{
    if(!warn){
      fit.gam <- suppressWarnings(
        gam(mod, data=data, family=twlss(), select=TRUE, ...)
      )
    } else{
      fit.gam <- gam(mod, data=data, family=twlss(), select=TRUE, ...)
    }
  }
  if(obl.corr){
    attr(fit.gam, "obl.corr") <- TRUE
  } else{
    attr(fit.gam, "obl.corr") <- FALSE
  }
  return(fit.gam)
}
