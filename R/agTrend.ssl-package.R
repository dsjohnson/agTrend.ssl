#' @title Estimating Steller Sea Lion Trends from Aerial Survey Data
#'
#' @description This package contains functions to fit GAM models to Steller sea lion survey data. It uses mgcv to fit Tweedie models and aggreates site estimates to regional trend estimates
#'
#' @note This software package is developed and maintained by scientists at the
#' NOAA Fisheries Alaska Fisheries Science Center and should be
#' considered a fundamental research communication. The recommendations and
#' conclusions presented here are those of the authors and this software should
#' not be construed as official communication by NMFS, NOAA, or the U.S. Dept.
#' of Commerce. In addition, reference to trade names does not imply endorsement
#' by the National Marine Fisheries Service, NOAA. While the best efforts have
#' been made to insure the highest quality, tools such as this are under
#' constant development and are subject to change.
#'
#' @name agTrend.ssl-package
#' @aliases agTrend.ssl-package agTrend.ssl
#' @author Devin S. Johnson <devin.johnson@@noaa.gov> (Maintainer)
"_PACKAGE"



.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
    paste(package, version, paste("(",date, ")", sep=""))
  )
}
