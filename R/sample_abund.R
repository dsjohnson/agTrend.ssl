################################################################################
### Simulate missing counts
################################################################################
#'@title Sample missing counts
#'@description Using a fitted GAM SSL model (see \code{\link{fit.gam}}) and the data used
#'for fitting, missing survey counts are imputed for aggregated trend and abundance
#'analysis.
#'@param fit A fitted model object from a call to \code{\link{fit.gam}}.
#'@param data Data used for model fitting.
#'@param yrs Years used for simulation. Defaults to a sequence of years from the
#'first to the last. But this can be any subset of years in the survey window.
#'@param size The sample size of the random imputation draws. Defaults to 1,000.
#'@param add.site.data An additional data set with a column lableled 'site' matching
#'the sites used in the model fitting data that contain additional groupings of other
#'site-level data.
#'@param keep.par Logical. Should the sample of parameters used be retained? It is not often needed
#'for further analysis.
#'@param debug Enter into the function for debugging.
#'
#'@author Devin S. Johnson
#'@export

sample.abund <- function(fit, data, yrs, size=1000, add.site.data=NULL, keep.par=FALSE, debug=FALSE){
  if(debug) browser()
  ff <- mgcv::interpret.gam(fit$formula)$pred.formula
  vars <- union(c('site','year','count'), attr(terms(ff), "term.labels"))
  if(attr(fit, "obl.corr")) obl.data <- data %>% select(site, year, obl)
  data <- data[,vars]
  data <- arrange(data, site, year)
  min.yr <- min(data$year)
  max.yr <- max(data$year)
  if(missing(yrs)) yrs <- min.yr:max.yr

  ### GAM parameters
  b <- coef(fit)
  V <- vcov(fit, unconditional=TRUE)

  ### Parameter indexes and matrices
  p.idx <- c(grep('\\(Intercept\\).1', names(b)), grep('s.1\\(', names(b)))
  phi.idx <- c(grep('\\(Intercept\\).2', names(b)), grep('s.2\\(', names(b)))
  mu.idx <- c(1:length(b))[!c(1:length(b)) %in% c(p.idx, phi.idx)]
  pred.data <- data %>% select(-year, -count) %>% distinct() %>% rowwise() %>%
    mutate(year=list(yrs)) %>% unnest(cols=year)
  if(attr(fit, "obl.corr")){
    pred.data <- left_join(pred.data, obl.data, by=c('site','year')) %>%
      mutate(obl=ifelse(is.na(obl), 0, obl))
  }
  bbb <- intersect(names(pred.data), names(data))
  pred.data <- pred.data %>% left_join(data,by=bbb) %>%
    group_by(.data$site) %>% nest() %>%
    mutate(
      data = map2(.data$data, .data$site, ~{.x$site=.y; .x}),
      X = map(data, ~{predict(fit, newdata=.x, type='lpmatrix')}),
      X.mu = map(X, ~{.x[,mu.idx]}),
      X.p = map(X, ~{unique(.x[,p.idx]) %>% as.matrix}),
      X.phi = map(X, ~{unique(.x[,phi.idx]) %>% as.matrix})
    ) %>% select(-X) %>% ungroup()

  #check for single p and phi par
  ck.p <- all(map_int(pred.data$X.p, nrow)==1)
  ck.phi <- all(map_int(pred.data$X.phi, nrow)==1)
  if(!(ck.phi&ck.p)) stop("There are multiple p and/or phi parameters for at least some sites. The sampler can't handle that at this time. Please refit the model so there is only 1 p and phi per site.")

  pred.data <- pred.data %>%
    mutate(
      data = map(.data$data, ~{attr(.x,"obl.corr")<-attr(fit,"obl.corr"); .x}),
      surv.times = map(.data$data, ~{unique(.x$year[!is.na(.x$count)])})
    )

  ### Sample parameters
  b.smp <- mvnfast::rmvn(size, b, V)
  b.mu <- b.smp[,mu.idx]
  b.p <- b.smp[,p.idx]
  b.phi <- b.smp[,phi.idx]
  pred.data$mu <- map(pred.data$X.mu, ~{exp(b.mu%*%t(.x))})
  pred.data$p <- map(pred.data$X.p, ~{plogis(b.p%*%t(.x))+1})
  pred.data$phi <- map(pred.data$X.phi, ~{exp(b.phi%*%t(.x))})
  pred.data$N.pred <- pmap(pred.data, draw.tw) %>% map(~{colnames(.x) <- yrs; .x})
  pred.data$N.real <- pmap(pred.data, get.real) %>% map(~{colnames(.x) <- yrs; .x})


  ### Clean up
  pred.data <- select(pred.data, -X.mu, -X.p, -X.phi)
  # pred.data$data <- map(pred.data$data, ~{select(.x, -site)})
  if(!keep.par) pred.data <- select(pred.data, -mu, -p, -phi)
  if(!is.null(add.site.data)){
    by <- names(add.site.data)[names(add.site.data)%in%names(pred.data)]
    pred.data <- pred.data %>% left_join(add.site.data, by=by)
  }
  rm.col <- which(colnames(pred.data$data[[1]])%in%colnames(pred.data))
  if(length(rm.col)>0) pred.data$data <- map(pred.data$data, ~{.x[,-rm.col]})

  attr(pred.data, "site.level") <- TRUE
  return(pred.data)
}
