#obtendo dataframe

df=readr::read_table(file = "HICOLM.thermodynamics",col_names = c("V1","V2"),n_max = 6)

df2=readr::read_table2(file="HICOLM.thermodynamics",col_names = TRUE,skip = 7)

df2 = df2 %>% gather(VOLUME,TEMPERATURE,PRESSURE,ENERGY,key="data",value="Value")
