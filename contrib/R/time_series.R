#install.packages("ggplot2","plotly")

#library(ggplot2)
#library(plotly)
library(tibble)
library(dplyr)

x11()

Sys.sleep(3)

meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
Optimization=ts(meusdados,frequency = 1)
i0=dim.data.frame(meusdados)
Sys.sleep(3)
meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
Optimization=ts(meusdados,frequency = 1)
i=dim.data.frame(meusdados)
plot(Optimization)

while(i[1]!=i0[1]) {
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  Optimization=ts(meusdados,frequency = 1)
  i0=dim.data.frame(meusdados)
  Sys.sleep(3)
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  Optimization=ts(meusdados,frequency = 1)
  i=dim.data.frame(meusdados)
  plot(Optimization)
}

dev.off()

cat("\n","Ending of simulation!")
cat("\n")

png("Optimization.png")
plot(Optimization)
summary(Optimization)

cat("\n")
