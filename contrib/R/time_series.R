
x11()

Sys.sleep(3)

meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
Optimization=ts(meusdados,frequency = 1)
i0=dim.data.frame(meusdados)
Sys.sleep(3)
meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
Optimization=ts(meusdados,frequency = 1)
i=dim.data.frame(meusdados)
plot(Optimization,main="Evolution of variables at each step",xlab="Step",col="red")

while(i[1]!=i0[1]) {
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  Optimization=ts(meusdados,frequency = 1)
  i0=dim.data.frame(meusdados)
  Sys.sleep(3)
  meusdados=read.table("HICOLM.df",sep = "",header = TRUE)
  Optimization=ts(meusdados,frequency = 1)
  i=dim.data.frame(meusdados)
  plot(Optimization,main="Evolution of variables at each step",xlab="Step",col="red")
}

dev.off()

cat("\n","Ending of simulation!")
cat("\n")

png("Optimization.png")
plot(Optimization,main="Evolution of variables at each step",xlab="Step",col="red")
summary(Optimization)

cat("\n")
