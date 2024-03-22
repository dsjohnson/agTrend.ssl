#' @title Process SSL count data .xlsx file
#' @description The current 'ALLCOUNTS_v(x).xlsx' file is processed into the
#' separate, R friendly, data sets necessary for trend analysis.
#' @param allcounts the file path of the current allcounts excel file.
#' @param age Age needed, must be one of 'pup' or 'nonpup'
#' @param dps DPS needed, must be one of 'wdps' or 'edps'
#' @author Devin S. Johnson
#' @export

proc.data <- function(allcounts, age, dps){
  out <- list()
  age <- tolower(age)
  dsp <- tolower(dps)
  if(dps=='edps'){
    if(age=='pup'){
      out <- read_excel(allcounts, sheet = "edpspup") %>%
        tidyr::gather(key="year", value="count", -SITE, -REGION, convert=T)%>%
        arrange(SITE, year)
      out$SITE <- stringi::stri_trans_general(out$SITE, "latin-ascii")
    } else{
      edps.photo <- read_excel(allcounts, sheet = "edpsnp_corr") %>%
        tidyr::gather(key="year", value="obl", -SITE, -REGION, convert=T) %>%
        arrange(SITE, year) %>% mutate(obl = ifelse(is.na(obl), 0, obl))
      edps.photo$SITE <- stringi::stri_trans_general(edps.photo$SITE, "latin-ascii")

      edpsnp <- read_excel(allcounts, sheet = "edpsnp") %>%
        tidyr::gather(key="year", value="count", -SITE, -REGION, convert=T)%>%
        arrange(SITE, year )
      edpsnp$SITE <- stringi::stri_trans_general(edpsnp$SITE, "latin-ascii")

      out <- full_join(edpsnp, edps.photo, by = c("SITE", "REGION", "year"))
    }
  } else{
    if(age=='pup'){
      out <- read_excel(allcounts, sheet = "wdpspup") %>%
        tidyr::gather(key="year", value="count", -SITE, -REGION, -RCA, convert=T)%>%
        arrange(SITE, year)
    } else{
      out <- read_excel(allcounts, sheet = "wdpsnp") %>%
        tidyr::gather(key="year", value="count", -SITE, -REGION, -RCA, convert=T)%>%
        arrange(SITE, year) %>% mutate(obl = as.integer(year<2004))
    }
  }
  colnames(out) <- tolower(colnames(out))
  out <- droplevels(out)
  out <- out[!is.na(out$count),]
  site.sum <- out %>% group_by(site) %>%
    summarize(
      n.survey = n(),
      n.nonzero = sum(count>0)
    )
  out <- left_join(out, site.sum, by='site')
  return(out)
}

