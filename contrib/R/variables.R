#obtendo informacoes da MD

library(dplyr)

df=readr::read_table2(file = "HICOLM.out",col_names = c("V1","V2","V3","V4"))

var1=slice(df,9L)

natoms=as.integer(var1[4])

var1=slice(df,27L)

ensemble=as.character(paste(var1[2],var1[3]))

var1=slice(df,3L)

var_host=as.character(var1[2])

var1=slice(df,4L)

var_date=as.integer(var1[2])

var1=slice(df,34L)

timestep=as.double(var1[2])

#var1=slice(df,27L)

#thermostat=as.numeric(var1[2])

#var1=slice(df,28L)

#temperature=as.numeric(var1[2])

#var1=slice(df,29L)

#pressure=as.numeric(var1[2])

# carregando dataframe

#dff=readr::read_table2(file="HICOLM.md",col_names =
#                         c("V1","V2","V3","V4","V5","V6","V7"),
#                       skip_empty_rows = TRUE)

#nxx=as.integer(dff[1,4])

#dfn=as.integer(3*nxx+9)

#df=dff[seq(4,nrow(dff),dfn),]
#df[,7]=dff[seq(5,nrow(dff),dfn),1]

df2=readr::read_table2(file="HICOLM.df",col_names = TRUE)
