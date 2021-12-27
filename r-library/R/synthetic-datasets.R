# Functions to help generate synthetic growth rate and R_t datasets using low-tech non model based approaches

# internal function to generate a random alternating high and low growth rate
alternatingRandom = function(n=6,meanHi=0.005,meanLo=-meanHi,sd=0.005,fn = rnorm) {
  as.vector(t(matrix(c(fn(n/2,mean=meanHi,sd=sd),fn(n/2,mean=meanLo,sd=sd)),nrow = n/2)))
}

#' Generates one or more discrete infectivity profiles based on discretised gamma distributions.
#'
#' @param mean - the mean(s) of the gamma distributions
#' @param sd - the sd(s) of the gamma distribution
#' @param names - the names of the infectivity profiles, if they represent a specific variant for example, if omitted they will be numbered
#' @param range - the number of days to generate - the first time point is always zero, so the resulting length will be one greater than this value.
#'
#' @return a list of class infectivity_profile containing aVector, the end points of the period, and yMatrix the discrete probabilities in the period up to the relevant aVector entry. Each column in yMatrix is a infectivity profile
#' @export
#'
#' @examples
#' discreteGammaInfectivityProfile(c(4,3),c(5,4),10, c("delta","omicron"))
discreteGammaInfectivityProfile = function(mean,sd,names=NULL,range=20) {
  shape = mean^2/sd^2
  rate = mean/sd^2
  a = rep(0:range,length(shape))
  pg = pgamma(a, shape, rate) - pgamma(a-1, shape, rate)
  pg = matrix(pg,ncol=length(shape))
  if (!is.null(names)) colnames(pg) = names
  ip = list(
    aVector = a,
    yMatrix = pg
  )
  class(ip) = "infectivity_profile"
  return(ip)
}

#' Generate a periodic growth timeseries simulation
#'
#' @param name - a meaningful name for this growth rate scenario
#' @param dateAtTime0 - the date at time zero
#' @param length - how many days in the scenario
#' @param period - the period of the periodic function
#' @param baseline - a baseline growth rate around which to oscillate usually 0.
#' @param smooth - smooth (sinusoidal) or triangular wave.
#' @param amplitude - growth rate (or R_t) amplitude
#' @param Gt.mean - mean of the generation time for conversion to R_t
#' @param Gt.sd - sd of the generation time for conversion to R_t
#'
#' @return - a list containing name - the scenario name,  ts - a tibble with the date, generated growth rate and R_t values, events - the time-points of the period, infectivityProfile - an infectivity_profile
#' @export
periodicGrowthRate = function(name, dateAtTime0 = "2020-01-01", length=365, period=91, baseline=0, smooth=TRUE, amplitude=0.1, Gt.mean=5, Gt.sd=4) {
  Gt.alpha = Gt.mean^2/Gt.sd^2
  Gt.beta = Gt.mean/Gt.sd^2
  if(smooth) {
    x = 0:(length-1)
    y = sin(x/period*(2*pi))*amplitude+baseline
    tmp = tibble(
      time = x,
      date = as.Date(dateAtTime0)+time,
      Growth.actual = y
    )
  } else {
    x = seq(0,length,by = period/4)
    y = rep(c(0,1,0,-1),length.out=length(x))
    int = approx(x=x, y=y, n = length, method="linear", rule=2)
    tmp = tibble(
      time = int$x,
      date = as.Date(dateAtTime0)+time,
      #Growth.actual = rates[cut(time,breaks = c(breaks,Inf), include.lowest = TRUE)]
      Growth.actual = int$y*amplitude+baseline
    )
  }
  tmp = tmp %>% mutate(
    Rt.actual = (1+Growth.actual/Gt.beta)^Gt.alpha,
    Gt.mean = Gt.mean,
    Gt.sd = Gt.sd
  )
  events = tibble(
    breaks = seq(0,length,by = period/4),
    rates = rep(c(0,1,0,-1),length.out=length(breaks))*amplitude+baseline,
    `Start date` = as.Date(dateAtTime0)+breaks,
    `End date` = NA, 
    Label = paste0("Growth: ",sprintf("%1.3f", rates))
  )
  #FIXME serial = #TODO SerialIntervalProvider$fixedGamma(controller,mean = Gt.mean, sd=Gt.sd)
  sim = list(
    name = name,
    ts = tmp %>% mutate(source=name), 
    events = events %>% mutate(source=name),
    infectivityProfile = discreteGammaInfectivityProfile(Gt.mean, Gt.sd)
  )
  class(sim) = "growth_simulation"
  return(sim)
}


#' Generate a growth timeseries simulation based on a small set of fixed time points.
#'
#' @param name - a meaningful name for this growth rate scenario
#' @param dateAtTime0 - the date at time zero
#' @param length - how many days in the scenario 
#' @param breaks - time points (between 0 and length) when the growth rate is at a set value.
#' @param rates - the values of the growth rate at the time points up to and including the break points, optionally an additional point for the period to the end of the time-series
#' @param smooth - should the rates be interpolated with a spline. If so it may be a good ideal to specify the growth rate at zero. 
#' @param sawtooth - if the rates are not smooth should they be a step function or a saw-tooth linear function.
#' @param Gt.mean - mean of the generation time for conversion to R_t
#' @param Gt.sd - sd of the generation time for conversion to R_t
#' @param ... - ignored
#'
#' @return - a list containing name - the scenario name,  ts - a tibble with the date, generated growth rate and R_t values, events - the time-points of the period, infectivityProfile - an infectivity_profile
generateGrowthRate = function(name, dateAtTime0 = "2020-01-01", length=365, breaks = sort(c(0,sample(0:364,4))), rates = alternatingRandom(), smooth=TRUE, sawtooth=FALSE, Gt.mean=5, Gt.sd=4, ...) {
  Gt.alpha = Gt.mean^2/Gt.sd^2
  Gt.beta = Gt.mean/Gt.sd^2
  if(length(breaks) == length(rates)-1) breaks=c(breaks,length-1)
  if(length(breaks) != length(rates)) stop("rates must be equal to or one longer than breaks")
  if(smooth) {
    sp = spline(x=breaks, y = rates, xout=0:(length-1))
    tmp = tibble(
      time = sp$x,
      date = as.Date(dateAtTime0)+time,
      Growth.actual = sp$y
    )
  } else {
    int = approx(x=breaks, y=rates, xout = 0:(length-1), method=(if(sawtooth) "linear" else "constant"),ties = "ordered",f=1, rule=2)
    tmp = tibble(
      time = int$x,
      date = as.Date(dateAtTime0)+time,
      #Growth.actual = rates[cut(time,breaks = c(breaks,Inf), include.lowest = TRUE)]
      Growth.actual = int$y
    )
  }
  tmp = tmp %>% mutate(
    Rt.actual = (1+Growth.actual/Gt.beta)^Gt.alpha,
    Gt.mean = Gt.mean,
    Gt.sd = Gt.sd
  )
  events = tibble(
    `Start date` = as.Date(dateAtTime0)+breaks,
    `End date` = NA, 
    Label = paste0("Growth: ",sprintf("%1.3f", rates)), 
    rates=rates, 
    breaks=breaks
    ) %>% 
    arrange(breaks)
  #FIXME serial = SerialIntervalProvider$fixedGamma(controller,mean = Gt.mean, sd=Gt.sd)
  
  sim = list(
    name = name,
    ts = tmp %>% mutate(source=name), 
    events = events %>% mutate(source=name),
    infectivityProfile = discreteGammaInfectivityProfile(Gt.mean, Gt.sd)
  )
  class(sim) = "growth_simulation"
  return(sim)
}

#' Generate an R_t timeseries simulation based on a small set of fixed time points.
#'
#' @param name - a meaningful name for this growth rate scenario
#' @param dateAtTime0 - the date at time zero
#' @param length - how many days in the scenario 
#' @param breaks - time points (between 0 and length) when the growth rate is at a set value.
#' @param R_t - the values of the reproduction at the time points up to and including the break points, optionally an additional point for the period to the end of the time-series
#' @param smooth - should the rates be interpolated with a spline. If so it may be a good ideal to specify the growth rate at zero. 
#' @param sawtooth - if the rates are not smooth should they be a step function or a saw-tooth linear function.
#' @param Gt.mean - mean of the generation time for conversion to R_t
#' @param Gt.sd - sd of the generation time for conversion to R_t
#' @param ... - ignored
#'
#' @return - a list containing name - the scenario name,  ts - a tibble with the date, generated growth rate and R_t values, events - the time-points of the period, infectivityProfile - an infectivity_profile
#' @export
generateRt = function(name, dateAtTime0 = "2020-01-01", length=365, breaks = sort(c(0,sample(0:364,5))), R_t = alternatingRandom(fn = rlnorm), smooth=TRUE, sawtooth=FALSE, Gt.mean=5, Gt.sd=4, ...) {
  Gt.alpha = Gt.mean^2/Gt.sd^2
  Gt.beta = Gt.mean/Gt.sd^2
  if(length(breaks) == length(R_t)-1) breaks=c(breaks,length-1)
  if(length(breaks) != length(R_t)) stop("R_t must be equal to or one longer than breaks")
  if(smooth) {
    sp = spline(x=breaks , y = log(R_t),xout= 0:(length-1))
    tmp = tibble(
      time = sp$x,
      date = as.Date(dateAtTime0)+time,
      Rt.actual = exp(sp$y)
    )
  } else {
    int = approx(x=breaks, y=R_t, xout = 0:(length-1), method=(if(sawtooth) "linear" else "constant"),f=1, rule=2)
    tmp = tibble(
      time = int$x,
      date = as.Date(dateAtTime0)+time,
      Rt.actual = int$y
    )
  }
  tmp = tmp %>% mutate(
    Growth.actual = (Rt.actual^(-Gt.alpha)-1)*Gt.beta,
    Gt.mean = Gt.mean,
    Gt.sd = Gt.sd
  )
  events = tibble(
    `Start date` = as.Date(dateAtTime0)+breaks,
    `End date` = NA, 
    Label = paste0("Rt: ",sprintf("%1.3f", rt)), 
    rt=rt, 
    breaks=breaks) %>% arrange(breaks)
  sim = list(
    name = name,
    ts = tmp %>% mutate(source=name), 
    events = events %>% mutate(source=name),
    infectivityProfile = discreteGammaInfectivityProfile(Gt.mean, Gt.sd)
  )
  class(sim) = "growth_simulation"
  return(sim)
}

#' Add importations to the simulation. Without importations there will be no growth so some cases must be imported.
#'
#' @param growthSimulation - a timeseries with an integer time column.
#' @param importDf - an optional dataframe with 2 columns: "time" and "import" counts.
#' @param rate - if no importDf is specified the rate parameter will define the number of cases on the zero day.
#'
#' @return a growth_simulation with an additional import column
#' @export
addImportations = function(growthSimulation, importDf = tibble(time = 0, import = rate), rate=100) {
  growthSimulation$ts = growthSimulation$ts %>% left_join(importDf, by="time") %>% mutate(import = ifelse(is.na(import),0,import))
  return(growthSimulation)
}

#' Calculate a simulated set of expected incidences based on imports and growth rate 
#'
#' @param growthSimulation - a timeseries with an Growth.actual column and a import column.
#'
#' @return the timeseries with an expected poisson rate column Est.actual.
#' @export
addPoissonRate = function(growthSimulation) {
  if(!("import" %in% colnames(growthSimulation$ts))) growthSimulation = growthSimulation %>% addImportations()
  x = growthSimulation$ts$import[[1]]
  for(i in 2:nrow(growthSimulation$ts)) {
    last = x[length(x)]
    x = c(x, last*exp(growthSimulation$ts$Growth.actual[[i-1]])+growthSimulation$ts$import[[i]])
  }
  growthSimulation$ts = growthSimulation$ts %>% mutate(Est.actual = x)
  return(growthSimulation)
}

#' Calculate a simulated set of expected proportions based on imports, growth rates and a baseline total
#'
#' @param growthSimulation - a timeseries with an Growth.actual column and a import column.
#' @param baselineExpr - an expression to derive the baseline total - either a single number or a way of calculating that number from growthSimulation
#'
#' @return the timeseries with an expected binomial proportion in the column Proportion.actual.
#' @export
addBinomialRate = function(growthSimulation,baselineExpr = rpois(n(),Est.actual)) {
  baselineExpr = enexpr(baselineExpr)
  if(!("Est.actual" %in% colnames(growthSimulation$ts))) growthSimulation = growthSimulation %>% addPoissonRate()
  growthSimulation$ts = growthSimulation$ts %>% mutate(total = !!baselineExpr) %>% mutate(Proportion.actual = Est.actual/total)
  # x = growthSimulation$ts$import[[1]]/growthSimulation$ts$total[[1]]
  # for(i in 2:nrow(growthSimulation$ts)) {
  #   last_p_t = x[length(x)]
  #   r_t = growthSimulation$ts$Growth.actual[[i-1]]
  #   x = c(x, 
  #         # p_{t+1} &= \frac{1}{1 + \frac{(1-p_t)}{p_t}e^{-r_t}}
  #         1/(1+(1-p_t)/p_t*exp(-r_t))+
  #           growthSimulation$ts$import[[i]]/growthSimulation$ts$total[[i]])
  # }
  # growthSimulation$ts = growthSimulation$ts %>% mutate(Proportion.actual = x)
  return(growthSimulation)
}

#' Add in a simulated rate of observed cases based on the underlying growth rate and imports with a weekly variation
#'
#' @param growthSimulation - the growth rate simulation
#' @param observedFractionExpr - what proportion of the cases are observed, either as a number, a vector of same length as the simulation, or a way of calculating that number from the growthSimulation. Typically under 1.
#' @param weekendEffect - what fraction of cases can we expect to not be observed on a saturday or sunday.
#'
#' @return the timeseries with an expected observed rate in the column Est.observed
#' @export
addObservedRate = function(growthSimulation, observedFractionExpr = 1, weekendEffect = 0.1) {
  observedFractionExpr = enexpr(observedFractionExpr)
  if(!("Est.actual" %in% colnames(growthSimulation$ts))) growthSimulation = growthSimulation %>% addPoissonRate()
  weekendProb = case_when(
    format(growthSimulation$ts$date,"%A") == "Saturday" ~ -1,
    format(growthSimulation$ts$date,"%A") == "Sunday" ~ -1,
    format(growthSimulation$ts$date,"%A") == "Monday" ~ 2,
    TRUE ~ 0
  )*weekendEffect+1
  growthSimulation$ts = growthSimulation$ts %>% 
    ungroup() %>%
    mutate(
      Fraction.observed = !!observedFractionExpr,
      weekday = format(date,"%A"),
      Fraction.weekend = case_when(
          weekday == "Saturday" ~ -1,
          weekday == "Sunday" ~ -1,
          weekday == "Monday" ~ 2,
          TRUE ~ 0)*weekendEffect+1
    ) %>% 
    mutate(Est.observed = Est.actual*Fraction.weekend*Fraction.observed)
  return(growthSimulation)
}

#' Generate a random set of observations from a poisson rate
#'
#' @param growthSimulation - the growth_simulation
#' @param bootstraps - the number of generated bootstraps
#' @param delayMean 
#' @param delaySd 
#' @param lastObservation 
#'
#' @return a growth_simulation with a time series including bootstrapped observations (as "value" column, bootstrap id in "subgroup")
#' @export
addBootstrappedObservations = function(growthSimulation, bootstraps = 100, delayMean = 0, delaySd = 1, lastObservation = nrow(growthSimulation$ts)-1) {
  if(!("Est.observed" %in% colnames(growthSimulation$ts))) growthSimulation = growthSimulation %>% addObservedRate()
  delayShape = delayMean^2/delaySd^2
  delayRate = delayMean/delaySd^2
  datedProb = growthSimulation$ts %>% filter(time %in% lastObservation) %>% select(endTime = time,observationDate = date) %>% distinct()
  if(delayMean == 0) {
    datedProb = datedProb %>% group_by_all() %>% summarise(
      time = 0:endTime,
      Fraction.available = 1
    )
  } else {
    datedProb = datedProb %>% group_by_all() %>% summarise(
      time = 0:endTime,
      Fraction.available = rev(pgamma((0:endTime+1), delayShape, delayRate))
    )
  }
  output = tibble(subgroup = 1:bootstraps) %>% inner_join(datedProb, by=character())
  
  growthSimulation$ts = growthSimulation$ts %>% 
    inner_join(output, by="time") %>%
    mutate(
      Est.delayed = Est.observed*Fraction.available,
      value=rpois(n(),Est.delayed)
    ) %>%
    mutate(statistic = "case",type = "incidence", code="XYZ", name="Test",codeType = "TEST",gender=NA_character_,ageCat=NA_character_)
  return(growthSimulation)
}

#' Convenience method to generate a random growth rate simulation
#'
#' @param growthSimulation - optional growth rate timeseries, if missing a random one is generated
#' @param bootstraps - ho many replicates
#' @param seed - what initial case loads
#' @param weekendEffect - any weekend effect?
#' @param Gt.mean - generation time mean
#' @param Gt.sd - generation time sd
#' @param periodic - a random periodic or
#' @param name - a name for the simulation
#' @param ... 
#'
#' @return a single strain random simulation with sensible parameters and no observation delay
#' @export
getGrowthRateBasedDataset = function(
  growthSimulation = NULL,
  bootstraps = 100, 
  seed=100, 
  weekendEffect = 0.1, 
  Gt.mean=5, Gt.sd=4, 
  periodic=FALSE, 
  name = "synthetic",
  ...) {
  
  if (is.null(growthSimulation)) {
    if(periodic) {
      growthSimulation = periodicGrowthRate(..., name = name, Gt.mean = Gt.mean, Gt.sd = Gt.sd) 
    } else {
      growthSimulation = generateGrowthRate(..., name = name, Gt.mean = Gt.mean, Gt.sd = Gt.sd)
    }
  }
  
  out = growthSimulation %>%
    addImportations(importDf = tibble(time = 0, import = seed)) %>%
    addPoissonRate() %>%
    addObservedRate(weekendEffect = weekendEffect) %>%
    addBootstrappedObservations(bootstraps = bootstraps,delayMean = 0, delaySd = 1)
  
  return(out)
}

#' Convenience wrapper to generate a 2 strain type comparison. strain 1 will have suffix .s1 and strain 2 will have .s2
#'
#' @param scenario1 a growthSimulation scenario for the s1 strain
#' @param scenario2 a growthSimulation scenario for the s2 strain 
#' @param delayMean1 a delay distribution parameter for the s1 strain
#' @param delayMean2 a delay distribution parameter for the s2 strain
#' @param delaySd1 a delay distribution parameter for the s1 strain
#' @param delaySd2 a delay distribution parameter for the s2 strain
#' @param timepoints times at which to generate observed time series.
#' @param ... 
#'
#' @return a 2 strain dataset in wide format.
getTwoAlternativesDataset = function(
  scenario1,
  scenario2,
  delayMean1, 
  delayMean2,
  delaySd1 = delayMean1/4, 
  delaySd2 = delayMean2/4,
  timepoints = nrow(scenario1$ts), #c(35,45,55,100)
  ...
) {
  
  scenario1Ts = scenario1 %>%
        addBootstrappedObservations(delayMean = delayMean1,delaySd = delaySd1, lastObservation = timepoints)

  scenario2Ts = scenario2 %>%
         addBootstrappedObservations(delayMean = delayMean1,delaySd = delaySd1, lastObservation = timepoints)
  
  combinedTs = scenario1Ts$ts %>% inner_join(scenario2Ts$ts, by=c("time","date","observationDate","subgroup","statistic","type","code","name","codeType","gender","ageCat"), suffix=c(".s1",".s2")) 
  combinedTs = combinedTs %>%
    mutate(
      Growth.relative = Growth.actual.s1 - Growth.actual.s2,
      Rt.advantage = Rt.actual.s1/Rt.actual.s2,
      Proportion.actual = Est.actual.s1 / (Est.actual.s1+Est.actual.s2),
      Total.actual = Est.actual.s1+Est.actual.s2,
      proportion = value.s1/(value.s1+value.s2),
      total = value.s1+value.s2
    )
    
  combinedEvents = bind_rows(scenario1$events,scenario2$events)
  ip = list() 
  ip[[scenario1$name]]=scenario1$infectivityProfile
  ip[[scenario2$name]]=scenario2$infectivityProfile
  
  sim = list(
    name = paste0(scenario1$name," vs ",scenario2$name),
    ts = combinedTs,
    events = combinedEvents,
    infectivityProfile = ip
  )
  
  class(sim) = "growth_simulation"
  return(sim)
}


