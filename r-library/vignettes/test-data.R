library(tidyverse)

devtools::load_all("~/Git/uk-covid-datatools/")
dpc = DataProviderController$setup("~/Data/maps")
dpc$loadSpimSources("~/S3/encrypted/")
tsp = dpc$timeseriesProcessor()

devtools::install("~/Git/r6-generator-maven-plugin-test/r-library/")

J = testRapi::JavaApi$new()

tmp = dpc$datasets$getPHEApiNations()
tmp = tmp %>% tsp$estimateRt()
tmpSi = SerialIntervalProvider$default(dpc)
tmpConfig = tmpSi$getBasicConfig(quick=FALSE, priorR0=1, priorR0Sd=2)

J$Serialiser$serialiseDataframe(dataframe = tmp,filename = "/home/terminological/Git/jepidemic/src/main/resources/pheApi.ser")
J$Serialiser$serialiseNamedList(dataframe = tmpConfig,filename = "/home/terminological/Git/jepidemic/src/main/resources/epiestimCovidConfig.ser")

# incidence = tmp %>% filter(name=="England" && statistic=="case") %>% ungroup() %>% select(dates=date, I=RollMean.value)
# siDist = EpiEstim::discr_si(0:14,4.1,2)
# stdConfig = EpiEstim::make_config(method = "non_parametric_si",si_distr = siDist, seed = 100, mean_prior = 1, std_prior = 2)
# EpiEstim::estimate_R(incidence, config = stdConfig)

library(EpiEstim)
## load data
data(Flu2009)
## incidence:
# head(Flu2009$incidence)
## serial interval (SI) distribution:
# Flu2009$si_distr

J$Serialiser$serialiseDataframe(dataframe = Flu2009$incidence, filename = "/home/terminological/Git/jepidemic/src/main/resources/flu2009.ser")

res <- estimate_R(Flu2009$incidence, config = make_config(method = "non_parametric_si",si_distr = Flu2009$si_distr, mean_prior=5, std_prior=5))

siFlu = tibble(
  day = 0:(length(Flu2009$si_distr)-1),
  prob = Flu2009$si_distr
)
J$Serialiser$serialiseDataframe(dataframe = siFlu, filename = "/home/terminological/Git/jepidemic/src/main/resources/flu2009SerialInterval.ser")

out = res$R %>% mutate(start_date = res$dates[t_start],end_date = res$dates[t_end])
J$Serialiser$serialiseDataframe(dataframe = out, filename = "/home/terminological/Git/jepidemic/src/main/resources/flu2009EpiEstim.ser")
