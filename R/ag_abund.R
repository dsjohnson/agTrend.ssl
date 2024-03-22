#'@title Aggregate abundance samples
#'@description Takes abundance samples created with the function
#'\code{\link{sample.abund}} and aggregates them based on a factor variable in
#'(or added) to the \code{sample.abund} data set.
#'@param x A sample abundance data frame produced by \code{\link{sample.abund}}
#'@param ag.var A factor (or character) variable in \code{x} that indicates the group membership
#'for aggregation.
#'@author Devin S. Johnson
#'@export
#'
ag.abund <- function(x, ag.var){
  results <- x %>% group_by(.data[[ag.var]]) %>% nest() %>%
    mutate(
      surv.times = map(.data$data, ~{sort(reduce(.x$surv.times, union))}),
      N.pred = map(.data$data, ~{reduce(.x$N.pred, `+`)}),
      N.real = map(.data$data, ~{reduce(.x$N.real, `+`)})
    )
  results <- results %>% ungroup() %>% select(-.data$data)
  attr(results, "site.level") <- FALSE
  return(results)
}
