library(agTrend.ssl)
library(tidyverse)
library(ggplot2)

### read in ALLCOUNTS file and process it. Sites are filtered out that had
### less than 2 surveys or less than 2 positive counts.
edpspup <- proc.data(allcounts="ALLCOUNTS_v8.xlsx", age="pup", dps="edps") %>%
  mutate(region=factor(region)) %>% filter(n.survey>1, n.nonzero>1)

### Get site information
site.info <- edpspup %>% select(site, region) %>% distinct()

### Fit Tweedie hierarchical GAM model for used for imputation
fit <- fit.gam(data=edpspup, obl.corr=FALSE,
               alt.mod=list(
                 mu.form = count ~  s(year, k=8) + s(year, site, bs="fs", k=8),
                 p.form=~s(region, bs='re'), phi.form=~s(region, bs='re')
               ), debug=F)

### Sample missing values using the fitted model
N <- sample.abund(fit, edpspup, yrs=1989:2019, size=5000, add.site.data=site.info)
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
reg.tr <- ag.trend(reg.N, timeframe=c(1989,2019), ci.prob=0.95)

reg.tr$growth

### Plot of fit
ggplot(data=reg.summ) +
  geom_ribbon(aes(x=year, ymin=ci.pred.lower, ymax=ci.pred.upper), alpha=0.1) +
  geom_path(aes(x=year, y=est.pred), alpha=0.2, lwd=1) +
  geom_path(aes(x=year, y=Est), data=reg.tr$fitted %>% filter(type=='predicted'), color='darkred', lwd=1.2) +
  geom_ribbon(aes(x=year, ymin=lower, ymax=upper), reg.tr$fitted %>% filter(type=='predicted'), alpha=0.2, fill='darkred') +
  geom_pointrange(aes(x=year, y=est.real, ymin=ci.real.lower, ymax=ci.real.upper), data=reg.summ %>% filter(surveyed==1)) +
  xlab('Year') + ylab('Count') + facet_wrap(~region, nrow=3, scale='free_y')

edps.N <- N %>% mutate(total="eDPS") %>% ag.abund( 'total')
edps.summ <- ag.summary(edps.N)
edps.tr <- ag.trend(edps.N, timeframe=c(1989,2019), ci.prob=0.95)

edps.tr$growth

ggplot(data=edps.summ) +
  geom_ribbon(aes(x=year, ymin=ci.pred.lower, ymax=ci.pred.upper), alpha=0.1) +
  geom_path(aes(x=year, y=est.pred), alpha=0.2, lwd=1) +
  geom_path(aes(x=year, y=Est), data=edps.tr$fitted %>% filter(type=='predicted'), color='darkred', lwd=1.2) +
  geom_ribbon(aes(x=year, ymin=lower, ymax=upper), edps.tr$fitted %>% filter(type=='predicted'), alpha=0.2, fill='darkred') +
  geom_pointrange(aes(x=year, y=est.real, ymin=ci.real.lower, ymax=ci.real.upper), data=edps.summ %>% filter(surveyed==1)) +
  xlab('Year') + ylab('Count') + ggtitle("eDPS pups")


