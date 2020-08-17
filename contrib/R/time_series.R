suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(readr))

Sys.sleep(3)

nt2=c(0)
while (nt2[1]==0) {
  nt=readr::read_table2(file = "HICOLM.out") %>% slice(27)
  nt2=as.numeric(nt[2])
}

i0=c(0,0)

while(i0[1]<1) {
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  i0=dim.data.frame(meusdados)
}

x11()

Sys.sleep(3)
meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
Optimization=ts(meusdados,frequency = 1)
i=dim.data.frame(meusdados)
plot(Optimization,main="Evolution of variable at each step",xlab="Step",col="red")

while(i[1]<nt2[1]) {
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  i0=dim.data.frame(meusdados)
  Sys.sleep(3)
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  Optimization=ts(meusdados,frequency = 1)
  i=dim.data.frame(meusdados)
  plot(Optimization,main="Evolution of variables at each step",xlab="Step",col="red")
}

cat("\n")

png("Optimization.png")
summary(Optimization)

cat("\n")
