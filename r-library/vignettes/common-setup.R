library(tidyverse)
library(patchwork)

here::i_am("vignettes/common-setup.R")

devtools::load_all("~/Git/uk-covid-datatools/")
ukcovidtools::reload("~/Git/uk-covid-datatools/config.yml")

output = ukcovidtools::output("~/Dropbox/sarscov2/r-estimation-methodology")

devtools::load_all("~/Git/standard-print-output/")
standardPrintOutput::setDefaults()

# load jepidemic
# try(detach("package:jepidemic", unload = TRUE),silent = TRUE)
# remove.packages("jepidemic")
# rm(list = ls()) may be required to clear old versions of the library code
# Restarting R maybe also required if there was a running java VM otherwise changes to the jars on the classpath are not picked up.
# install locally compiled R library:
# devtools::install("~/Git/jepidemic/r-library/", upgrade = "never")
J = jepidemic::JavaApi$new()

# TODO: if I am submitting to CRAN I'll need to remove these dependencies
synth = SyntheticDatasetProvider$new(dpc)

source(here::here("vignettes/validation-functions.R"))

## Plots ----

rtPlot = function(out, ylim=c(0,NA)) {
  ggplot(out %>% filter(!is.nan(Rt.Mean)),aes(x=Rt.EndDate, y=Rt.Mean, colour=as.factor(Rt.Window)))+
    geom_point(size=0.5) +
    geom_errorbarh(aes(x=Rt.EndDate, xmin=Rt.StartDate, xmax=Rt.EndDate), alpha=0.3) +
    geom_errorbar(aes(y=Rt.Quantile.0.5, ymin = Rt.Quantile.0.025, ymax = Rt.Quantile.0.975), alpha=0.3)+
    guides(colour="none")+coord_cartesian(ylim=ylim)+
    geom_hline(yintercept=1,colour="grey50",inherit.aes=FALSE)+
    xlab("date")+
    ylab(latex2exp::TeX("$R_t$"))
}

rtPanel = function(out, ylim=c(0,5), colourExpr=expr(as.factor(Rt.Window))) {
  ggplot(out %>% filter(!is.nan(Rt.Mean)), aes(x=Rt.EndDate, y=Rt.Mean, colour=!!colourExpr))+
    geom_point(size=1) +
    geom_errorbarh(aes(x=Rt.EndDate, xmin=Rt.StartDate, xmax=Rt.EndDate), alpha=0.3,colour = "grey") +
    geom_errorbar(aes(y=Rt.Quantile.0.5, ymin = Rt.Quantile.0.025, ymax = Rt.Quantile.0.975), alpha=0.3,colour = "grey")+
    guides(colour="none")+coord_cartesian(ylim=ylim)
}

rtRibbon = function(out, ylim=c(0,NA), mapping = aes(), ...) {
  defaultAes = aes(x=Rt.EndDate, y=Rt.Mean)
  dots = rlang::list2(...)
  combinedAes = modifyList(defaultAes,mapping)
  ggplot(out %>% filter(!is.nan(Rt.Mean)), mapping = combinedAes)+
    geom_ribbon(aes(ymin = Rt.Quantile.0.025, ymax = Rt.Quantile.0.975), alpha=0.2,colour = NA,fill="grey20")+
    geom_ribbon(aes(ymin = Rt.Quantile.0.25, ymax = Rt.Quantile.0.75), alpha=0.3,colour = NA,fill="grey20")+
    geom_line(...) +
    geom_hline(yintercept=1,colour="grey50",inherit.aes=FALSE)+
    guides(colour="none")+coord_cartesian(ylim=ylim)
}

incidencePlot = function(out, ilim = c(NA,NA)) {
  ggplot(out, aes(x=date,y=value,group=subgroup))+
    geom_point(size=0.5)+
    geom_line(alpha=0.1)+
    ylab("cases")
}

plotRollingQuantiles = function(data, quantileVar, orderVar = "date", colours=c("red","blue"), window=14) {
  quantileVar = ensym(quantileVar)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  
  overallLag = data %>% summarise(median = quantile(!!quantileVar,0.5,na.rm = TRUE))
  
  summaryLag = data %>% arrange(!!orderVar) %>% 
    mutate(rollingQuant = slider::slide_index(!!quantileVar, !!orderVar, .before=window,.after=window,.f = ~ enframe(quantile(.x, c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)))
    ) %>%
    select(!!!grps, !!orderVar,rollingQuant) %>%
    unnest(rollingQuant) %>%
    distinct()
  
  list(
    geom_line(data = summaryLag %>% filter(name %in% c("2.5%","97.5%")), mapping=aes(x=!!orderVar,y=value, group=name), linetype="dotted", colour=colours[1] ),
    geom_line(data = summaryLag %>% filter(name %in% c("25%","75%")), mapping=aes(x=!!orderVar,y=value, group=name), linetype="dashed", colour=colours[1] ),
    geom_line(data = summaryLag %>% filter(name %in% c("50%")), mapping=aes(x=!!orderVar,y=value, group=name), linetype="solid", colour=colours[1] ),
    geom_hline(data = overallLag, mapping = aes(yintercept = median), colour=colours[2])
  )
}

plotRollingDeciles = function(data, decileVar, orderVar = "date", colours="grey50", window=14) {
  decileVar = ensym(decileVar)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  
  summaryLag = data %>% arrange(!!orderVar) %>% 
    mutate(rollingQuant = slider::slide_index(!!decileVar, !!orderVar, .before=window,.after=window,.f = ~ enframe(quantile(.x, seq(0.1,0.9,length.out = 5), na.rm=TRUE)))
    ) %>%
    select(!!!grps, !!orderVar,rollingQuant) %>%
    unnest(rollingQuant) %>%
    distinct()
  
  list(
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=value, group=name), colour=colours, alpha=0.5)
  )
}

plotRollingProportion = function(data, binomialExpr, orderVar = "date", colours=c("red","blue"), window=14) {
  binomialExpr = enexpr(binomialExpr)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  
  overallLag = data %>% summarise(p = sum(!!binomialExpr,na.rm = TRUE)/length(na.omit(!!binomialExpr)))
  
  summaryLag = data %>% arrange(!!orderVar) %>% 
    mutate(
      x = slider::slide_index_dbl(!!binomialExpr, !!orderVar, .before=window,.after=window,.f = ~ sum(.x,na.rm = TRUE)),
      n = slider::slide_index_dbl(!!binomialExpr, !!orderVar, .before=window,.after=window,.f = ~ length(na.omit(.x))),
    ) %>%
    select(!!!grps, !!orderVar, x, n) %>%
    distinct() 
  
  # browser()
  
  summaryLag = summaryLag %>%
    filter(!is.na(x) & !is.na(n)) %>%
    mutate(binom::binom.confint(x,n,methods="wilson")) %>%
    rename(lower.0.025 = lower,upper.0.975 = upper) %>% select(-mean) %>%
    mutate(binom::binom.confint(x,n,methods="wilson",conf.level=0.5)) %>%
    rename(lower.0.25 = lower,upper.0.75 = upper)
  
  list(
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=lower.0.025), linetype="dotted", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=upper.0.975), linetype="dotted", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=lower.0.25), linetype="dashed", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=upper.0.75), linetype="dashed", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=mean), linetype="solid", colour=colours[1] ),
    geom_hline(data = overallLag, mapping = aes(yintercept = p), colour=colours[2])
  )
}


plotRollingDensity = function(data, observationVar, n=50, orderVar = "date", colours=c("red","blue"), window=14) {
  observationVar = ensym(observationVar)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  obs = data %>% pull(!!observationVar)
  breaks = seq(min(obs,na.rm = TRUE),max(obs,na.rm = TRUE),length.out = n+1)
  miBreak = na.omit(breaks+lag(breaks))/2
  summaryLag = data %>% select(!!!grps,!!orderVar,!!observationVar) %>% nest(observations = !!observationVar) %>% arrange(!!orderVar) 
  summaryLag = summaryLag %>% mutate(bins = slider::slide(observations, .before=14, .after=14, .f = ~bind_rows(.x) %>% 
                                                       mutate(bin = cut(!!observationVar,breaks = breaks,labels = miBreak)) %>% 
                                                       group_by(bin) %>% 
                                                       summarise(count=n()) %>% 
                                                       mutate(prop = count/sum(count)) %>% 
                                                       ungroup() %>% 
                                                       tidyr::complete(bin, fill=list(count=0,prop=0)) ))
  summaryLag = summaryLag %>% select(!!!grps,!!orderVar,bins) %>% unnest(bins)
  
  list(
    geom_tile(data = summaryLag,  mapping=aes(x=!!orderVar, y=as.numeric(as.character(bin)), fill=prop))
  )
}

printMedianAndCI = function(x) do.call(sprintf, c("%1.2f [95%% CI %1.2f \u2013 %1.2f]",as.list(quantile(x,c(0.5,0.025,0.975),na.rm = TRUE))))
printIQR = function(x) do.call(sprintf, c("%1.2f [IQR %1.2f \u2013 %1.2f]",as.list(quantile(x,c(0.5,0.25,0.75),na.rm = TRUE))))

## EpiEstim data ----
# Flu2009 data
data("Flu2009", package="EpiEstim")


## Synthetic data ----
# generate new
# set.seed(101)
# tmp = synth$getGrowthRateBasedDataset(weekendEffect = 0,smooth=FALSE,seed = 100, bootstraps=10, periodic = TRUE)
# incidencePlot(tmp$ts)
# paste0("c(",paste0(tmp$events$rates,collapse = ", "),")") %>% clipr::write_clip()
# paste0("c(",paste0(tmp$events$breaks,collapse = ", "),")") %>% clipr::write_clip()

synthetic = standardSyntheticDataset()

# set.seed(101)
# options2 = tibble(weekendEffect = c(0,0.03,0.1)) %>% 
#   inner_join(tibble(seed = c(10,1000)), by=character()) %>%
#   inner_join(tibble(smooth = c(TRUE,FALSE)), by=character()) %>%
#   inner_join(tibble(baseline = c(0,-0.005)), by=character())
# 
# periodicAsTibble = function(config) {
#   tmp = synth$getGrowthRateBasedDataset(weekendEffect = config$weekendEffect, smooth= config$smooth, seed = config$seed, baseline=config$baseline,
#                                         periodic=TRUE, bootstraps=10)
#   tmp$serial = NULL
#   enframe(tmp) %>% pivot_wider(names_from = name, values_from=value)
# }
# 
# periodic = options2 %>% group_by_all() %>% group_modify(function(d,g,...) {
#   return(periodicAsTibble(g))
# })

## misc processing functions ----

## Divergence and comparison metrics ----
# KL divergence for 

klGammaVsGamma = function(mean1, sd1, mean2, sd2) {
  if (sd1>mean1 | sd2>mean2) stop("sd must be smaller than mean")
  range = seq(0,10,length=200)
  g1 = dgamma(range, shape=mean1^2/sd1^2,rate = mean1/sd1^2)
  g2 = dgamma(range, shape=mean2^2/sd2^2,rate = mean2/sd2^2)
  kl = flexmix::KLdiv(cbind(g1=g1,g2=g2))[1,2]
  # A matrix of KL divergences where the rows correspond to using the respective distribution as  in the formula above.
  return(kl)
}

klGammaVsUnif = function(value1, precision1, mean2, sd2) {
  if (sd2>mean2) stop("sd must be smaller than mean")
  range = seq(0,10,length=200)
  g1 = dunif(range, min = value1-precision1,max = value1+precision1)
  g2 = dgamma(range, shape=mean2^2/sd2^2,rate = mean2/sd2^2)
  kl = flexmix::KLdiv(cbind(g1=g1,g2=g2))[1,2]
  # A matrix of KL divergences where the rows correspond to using the respective distribution as  in the formula above.
  return(kl)
}


klGammaVsGamma(mean1=2,sd1=1.5,mean2=2.3,sd2=2.1)
klGammaVsUnif(value1=2,precision1=0.0001,mean2=2,sd2=0.002)