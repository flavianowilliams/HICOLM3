library("magrittr")
library("tibble")
library("readr")

suppressMessages(library("dplyr"))

Sys.sleep(5)

df=readr::read_table2(file = "HICOLM.out",col_names = c("V1","V2"))

mtd=df %>% slice(25L)

opcao=paste(mtd[1],mtd[2])

if (opcao=="Steepest descent") {
  nt=df %>% slice(27L)
} else if (opcao=="Molecular dynamics") {
  nt=df %>% slice(31L)
}

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
  plot(Optimization,main="Evolution of variable at each step",xlab="Step",col="red")
  Sys.sleep(3)
}

cat("\n")

png("Optimization.png")
summary(Optimization)

cat("\n")
