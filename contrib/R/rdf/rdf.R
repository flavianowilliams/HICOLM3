rm(list = ls())

graphics.off()

setwd(getwd())

library(tidyverse)

source("import_and_clean_data.R",local = knitr::knit_global())
source("model.R",local = knitr::knit_global())

Sys.getenv(c("rdf_a1","rdf_a2"))

a1="rdf_a1"
a2="rdf_a2"
