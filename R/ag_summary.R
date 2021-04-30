################################################################################
### Summarize simulated counts
################################################################################
#'@title Summarize count posterior distributions
#'@description The posterior distribution of the count sample is summarized for
#'each site. The resulting data frame can be used for site-level plots, etc.
#'@param x A data frame created with the sampling function \code{\link{sample.abund}} or
#'\code{\link{ag.abund}}
#'@param ci.prob Probability level for credible intervals. Defaults to 0.9.
#'@author Devin S. Johnson
#'@export
ag.summary <- function(x, ci.prob=0.9){
  site.level <- attr(x, "site.level")
  ### Summarize imputed counts
  df <- select(x, -.data$surv.times, -.data$N.pred, -.data$N.real)
  df$est.pred <- map(x$N.pred, ~{apply(.x, 2, median)})
  df$se.pred <-  map(x$N.pred, ~{apply(.x, 2, sd)})
  df$ci.pred <- map(x$N.pred, ~{
    coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>%
      `colnames<-`(c("ci.pred.lower", "ci.pred.upper"))
  })
  df$est.real <- map(x$N.real, ~{apply(.x, 2, median)})
  df$se.real <-  map(x$N.real, ~{apply(.x, 2, sd)})
  df$ci.real <- map(x$N.real, ~{
    coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>%
      `colnames<-`(c("ci.real.lower", "ci.real.upper"))
  })
  if(site.level){
    df <- unnest(df, cols=c('data','est.pred','se.pred','ci.pred','est.real','se.real','ci.real'))
  } else{
    nc <- as.numeric(colnames(x$N.pred[[1]]))
    df$surveyed <- map(x$surv.times, ~{1.0*(nc%in%.x)})
    df <- unnest(df, cols=c('est.pred','se.pred','ci.pred','est.real','se.real','ci.real','surveyed'))
    df$year <- rep(as.numeric(colnames(x$N.pred[[1]])), nrow(x))
    df <- df[,c(1,ncol(df), 2:(ncol(df)-1))]
  }
  return(df)
}
