library(tidyverse)
here::i_am("data-raw/synthetic-data-builder.R")
J = testRapi::JavaApi$get()

devtools::load_all("~/Git/uk-covid-datatools/")
ukcovidtools::reload("~/Git/uk-covid-datatools/config.yml")
synth = SyntheticDatasetProvider$new(dpc)
# generate new

# tmp = synth$getGrowthRateBasedDataset(weekendEffect = 0,smooth=FALSE,seed = 1000)
# tsp$plotIncidenceQuantiles(tmp$ts, events = tmp$events)
# 
# paste0("c(",paste0(tmp$events$rates,collapse = ", "),")") %>% clipr::write_clip()
# paste0("c(",paste0(tmp$events$breaks,collapse = ", "),")") %>% clipr::write_clip()

# kept

syntheticTimeseries = sapply(c(0,0.03,0.1), USE.NAMES = TRUE, simplify=FALSE, FUN=function(weekend) {
  tmp2 = synth$getGrowthRateBasedDataset(
    rates = c(0.0213732805102556, -0.0261293028360427, -0.00472741863318508, 0.0421778375461848, -0.042958304850032, -0.0693039179026735),
    breaks = c(0, 45, 86, 197, 254, 351),
    weekendEffect = weekend
  )
  
  # J$Serialiser$serialiseDataframe(dataframe = tmp2$ts, filename = paste0("/home/terminological/Git/jepidemic/src/main/resources/synthetic_timeseries_",weekend,".ser"))
  tmp2$serial = NULL
  return(tmp2)
})
names(syntheticTimeseries) = c(0,0.03,0.1)
usethis::use_data(syntheticTimeseries)

tsp$plotIncidenceQuantiles(tmp2$ts, events = tmp2$events)
J$Serialiser$serialiseNamedList(list(inf = syntheticTimeseries$`0`$infectivityProfile, events = syntheticTimeseries$`0`$events), filename = "/home/terminological/Git/jepidemic/src/main/resources/synthetic_cfg.ser")
