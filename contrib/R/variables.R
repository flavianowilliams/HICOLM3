# carregando dataframe

dff=readr::read_table2(file="HICOLM.md",col_names =
                         c("V1","V2","V3","V4","V5","V6"),
                       skip_empty_rows = TRUE)

nxx=as.integer(dff[1,4])

dfn=as.integer(3*nxx+7)

df=dff[seq(4,nrow(dff),dfn),]
df[,7]=dff[seq(5,nrow(dff),dfn),1]

names(df)[7]="V7"

#obtendo informacoes da MD

library(dplyr)

df2=readr::read_table2(file = "HICOLM.out")

var1=slice(df2,26L)

ensemble=as.character(paste(var1[2],var1[3]))

var1=slice(df2,27L)

#thermostat=as.numeric(var1[2])

#var1=slice(df2,28L)

#temperature=as.numeric(var1[2])

#var1=slice(df2,29L)

#pressure=as.numeric(var1[2])
