library(tidyverse)
library(patchwork)

here::i_am("vignettes/common-setup.R")
# load jepidemic
# try(detach("package:jepidemic", unload = TRUE),silent = TRUE)
# remove.packages("jepidemic")
# rm(list = ls()) may be required to clear old versions of the library code
# Restarting R maybe also required if there was a running java VM otherwise changes to the jars on the classpath are not picked up.
# install locally compiled R library:
# devtools::install("~/Git/jepidemic/r-library/", upgrade = "never")
devtools::load_all()
source(here::here("vignettes/validation-functions.R"))
source(here::here("vignettes/plotting-functions.R"))

output = outputter("~/Dropbox/sarscov2/r-estimation-methodology")

devtools::load_all("~/Git/standard-print-output/")
standardPrintOutput::setDefaults()

J = jepidemic::JavaApi$new()

## EpiEstim data ----
data("Flu2009", package="EpiEstim")




