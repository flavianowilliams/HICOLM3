#obtendo dataframe

df=readr::read_table(file = "HICOLM.thermodynamics",col_names = FALSE,n_max = 1)

df2=readr::read_table2(file="HICOLM.thermodynamics",col_names = TRUE,skip = 2)

df2 = df2 %>% gather(VOLUME,TEMPERATURE,PRESSURE,ENERGY,DENSITY,key="data",value="Value")
