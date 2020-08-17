suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(readr))

Sys.sleep(3)

nt=readr::read_table2(file = "HICOLM.out") %>% slice(27)
nt2=as.numeric(nt[2])

x11()

i=c(0,0)
while (i[1]<1) {
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  Optimization=ts(meusdados,frequency = 1)
  i=dim.data.frame(meusdados)
}

plot(Optimization,main="Evolution of variable at each step",xlab="Step",col="red")
Sys.sleep(3)

while(i[1]<nt2[1]) {
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  Optimization=ts(meusdados,frequency = 1)
  i=dim.data.frame(meusdados)
  plot(Optimization,main="Evolution of variables at each step",xlab="Step",col="red")
  Sys.sleep(3)
}

cat("\n")

png("Optimization.png")
summary(Optimization)

cat("\n")
