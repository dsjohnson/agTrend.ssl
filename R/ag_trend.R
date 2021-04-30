################################################################################
### Estimate trend
################################################################################
#'@title Estimate trends of aggregated abundance
#'@description The growth trend of is estimated from the samples obtained from the
#'functions \code{\link{sample.abund}} and \code{\link{ag.abund}}.
#'@param x A abundance sample data frame from \code{\link{sample.abund}} or \code{\link{ag.abund}}
#'@param timeframe a 2-vector giving the start year and end year of the desired
#'trend estimate, e.g., \code{timeframe = c(1989,2019)}.
#'@param ci.prob The probability for the credible interval. Defaults to 0.9.
#'@details The function returns a named list with elements \code{growth}, \code{fitted},
#'and \code{sample}. The \code{growth} element contains a table with the estimated
#'growth of each aggregation in *percent growth* form. The \code{fitted} element contains a
#'table with the fitted trendline values on the count scale. This can be used for plotting.
#'Finally, the \code{sample} element contains the posterior sample. This can be used for further
#'analysis of trends such as comparisons, etc.
#'@export
#'@author Devin S. Johnson
ag.trend <- function(x, timeframe, ci.prob=0.9){
  if(is.null(dim(timeframe))) timeframe <- t(as.matrix(timeframe))
  nms <- pull(x,1)
  cn <- colnames(x)[1]
  yrs <- as.numeric(colnames(x$N.pred[[1]]))
  zeros.pred <- map_lgl(x$N.pred, ~{any(.x==0)})
  zeros.real <- map_lgl(x$N.real, ~{any(.x==0)})
  if(any(zeros.pred) | any(zeros.real)) warning(paste0("There are abundnce values of 0 for: ", nms[zeros.pred | zeros.real], ". Adding '0.5' to allow calculation."))
  ln.Np <- map2(x$N.pred, zeros.pred, ~{log(as.matrix(.x + 0.5*.y))})
  ln.Nr <- map2(x$N.real, zeros.real, ~{log(as.matrix(.x + 0.5*.y))})
  out <- x[,1]
  idx <- yrs>=timeframe[1] & yrs<=timeframe[2]
  H <- cbind(1,1:sum(idx))
  P <-  H %>% {solve(crossprod(.), t(.))}
  out$trend.pred <- map(ln.Np, ~{.x[,idx]}) %>% map(~{t(P%*%t(.x))}) %>% map(coda::mcmc)
  out$fitted.pred <- map(out$trend.pred, ~{exp(t(H%*%t(.x)))}) %>% map(coda::mcmc)
  out$growth.pred <-  map(out$trend.pred, ~{100*(exp(.x[,2])-1)}) %>% map(coda::mcmc)
  out$trend.real <- map(ln.Nr, ~{.x[,idx]}) %>% map(~{t(P%*%t(.x))}) %>% map(coda::mcmc)
  out$fitted.real <- map(out$trend.real, ~{exp(t(H%*%t(.x)))}) %>% map(coda::mcmc)
  out$growth.real <-  map(out$trend.real, ~{100*(exp(.x[,2])-1)}) %>% map(coda::mcmc)
  summ.fitted.pred <- select(out, .data[[cn]], .data$fitted.pred) %>%
    mutate(
      year = rep(list(yrs[idx]),nrow(x)),
      type = 'predicted',
      Est = map(.data$fitted.pred, ~{apply(.x,2,median)}),
      CI = map(.data$fitted.pred, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-.data$fitted.pred) %>% unnest(cols=c('year', 'Est', 'CI'))
  summ.fitted.real <- select(out, .data[[cn]], .data$fitted.real) %>%
    mutate(
      year = rep(list(yrs[idx]),nrow(x)),
      type = 'realized',
      Est = map(.data$fitted.real, ~{apply(.x,2,median)}),
      CI = map(.data$fitted.real, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-.data$fitted.real) %>% unnest(cols=c('year', 'Est', 'CI'))
  summary.fitted <- bind_rows(summ.fitted.pred, summ.fitted.real)
  summ.growth.pred <- select(out, .data[[cn]], .data$growth.pred) %>%
    mutate(
      type='predicted',
      Est = map(.data$growth.pred, median),
      CI = map(.data$growth.pred, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-.data$growth.pred) %>% unnest(cols=c('Est', 'CI'))
  summ.growth.real <- select(out, .data[[cn]], .data$growth.real) %>%
    mutate(
      type='realized',
      Est = map(.data$growth.real, median),
      CI = map(.data$growth.real, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-.data$growth.real) %>% unnest(cols=c('Est', 'CI'))
  summary.growth <- bind_rows(summ.growth.pred, summ.growth.real)

  return(list(growth=summary.growth, fitted=summary.fitted, sample=out))
}
