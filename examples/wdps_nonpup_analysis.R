library(agTrend.ssl)
library(tidyverse)
library(ggplot2)

### read in ALLCOUNTS file and process it. Sites are filtered out that had
### less than 2 surveys or less than 2 positive counts.
wdpsnp <- proc.data(allcounts="ALLCOUNTS_v8.xlsx", age="np", dps="wdps") %>%
  filter(n.survey>1, n.nonzero>1) %>% mutate(obl = as.integer(year<2004))

### Get site information
site.info <- wdpsnp %>% select(site, region, rca) %>% distinct() %>%
  mutate(
    zone = ifelse(region%in%c("W ALEU","C ALEU","E ALEU"), "ALEU", "GULF"),
    wc_aleu = ifelse(region%in%c("W ALEU","C ALEU"), "W/C ALEU", NA)
  )

### Fit Tweedie hierarchical GAM model for used for imputation
fit <- fit.gam(data=wdpsnp, obl.corr=TRUE)

### Sample missing values using the fitted model
set.seed(123) # setting seed to make sure draws are reproducible
N <- sample.abund(fit, wdpsnp, yrs=1989:2019, size=1000, add.site.data=site.info)
N.summ <- ag.summary(N, ci.prob=0.9)

### Run the section below to plot all sites
# for(i in 1:nrow(site.info)){
#   df <- filter(N.summ, site==site.info$site[i])
#   p <- ggplot(data=df) +
#     geom_ribbon(aes(x=year, ymin=ci.pred.lower, ymax=ci.pred.upper), alpha=0.2) +
#     geom_path(aes(x=year, y=est.pred), alpha=1, color='darkred', lwd=1.2) +
#     geom_point(aes(x=year, y=est.real), data=df %>% filter(!is.na(count)))+
#     ggtitle(df$site[[1]]) + xlab('Year') + ylab('Count')
#   print(p)
# }

reg.N <- ag.abund(N, 'region')
reg.summ <- ag.summary(reg.N)
reg.tr <- ag.trend(reg.N, timeframe=c(2000,2019), ci.prob=0.95)

reg.tr$growth

### Plot of fit
ggplot(data=reg.summ) +
  geom_ribbon(aes(x=year, ymin=ci.pred.lower, ymax=ci.pred.upper), alpha=0.1) +
  geom_path(aes(x=year, y=est.pred), alpha=0.2, lwd=1) +
  geom_path(aes(x=year, y=Est), data=reg.tr$fitted %>% filter(type=='predicted'), color='darkred', lwd=1.2) +
  geom_ribbon(aes(x=year, ymin=lower, ymax=upper), reg.tr$fitted %>% filter(type=='predicted'), alpha=0.2, fill='darkred') +
  geom_pointrange(aes(x=year, y=est.real, ymin=ci.real.lower, ymax=ci.real.upper), data=reg.summ %>% filter(surveyed==1)) +
  xlab('Year') + ylab('Count') + facet_wrap(~region, nrow=3, scale='free_y')

wdps.N <- N %>% mutate(total="wdps") %>% ag.abund( 'total')
wdps.summ <- ag.summary(wdps.N)
wdps.tr <- ag.trend(wdps.N, timeframe=c(2000,2019), ci.prob=0.95)

wdps.tr$growth

ggplot(data=wdps.summ) +
  geom_ribbon(aes(x=year, ymin=ci.pred.lower, ymax=ci.pred.upper), alpha=0.1) +
  geom_path(aes(x=year, y=est.pred), alpha=0.2, lwd=1) +
  geom_path(aes(x=year, y=Est), data=wdps.tr$fitted %>% filter(type=='predicted'), color='darkred', lwd=1.2) +
  geom_ribbon(aes(x=year, ymin=lower, ymax=upper), wdps.tr$fitted %>% filter(type=='predicted'), alpha=0.2, fill='darkred') +
  geom_pointrange(aes(x=year, y=est.real, ymin=ci.real.lower, ymax=ci.real.upper), data=wdps.summ %>% filter(surveyed==1)) +
  xlab('Year') + ylab('Count') + ggtitle("wDPS nonpups")
