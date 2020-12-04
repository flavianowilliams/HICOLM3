#obtendo dataframe

setwd("/home/flaviano/Documentos/GitHub/HICOLM/contrib/R/report")

head=readr::read_table(file = "HICOLM.thermodynamics",col_names = FALSE,n_max = 1,skip = 3)

thermo=readr::read_table2(file="HICOLM.thermodynamics",col_names = TRUE,skip = 5)

thermo = thermo %>% gather(VOLUME,TEMPERATURE,PRESSURE,ENERGY,DENSITY,key="data",value="Value")

atoms=readr::read_table2(file="HICOLM.atoms",col_names = TRUE,skip = 5)
