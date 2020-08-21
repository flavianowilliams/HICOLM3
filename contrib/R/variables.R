#carregando pacotes

library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)
library(ggplot2)

#setwd("/home/flaviano/Dropbox/Projeto/simulacoes/H2O_MD")
#setwd("/home/flaviano/Documentos/maria_eduarda")

# carregando dataframe

nx=readr::read_table2(file="HICOLM.md",col_names = FALSE,n_max=1)

dff=data.frame(c("A1","X2","X3","X4","X5","X6","X7"))

dff=readr::read_table2(file="HICOLM.md",col_names = FALSE,skip = 2,skip_empty_rows = TRUE)
dff$X7=c(seq(1,dim(dff)[1],1))

n1=as.integer(nx[1,"X1"])
n2=as.numeric(nx[1,"X2"])
n3=as.integer(nx[1,"X3"])
n4=as.integer(nx[1,"X4"])
n5=as.integer(nx[1,"X5"])

nxx=c(n1,n2,n3,n4,n5)

dfn=as.integer(3*nxx[4]+7)

dt=dff %>% filter(mod(dff$X7,dfn)==2) %>% select("X1") %>% rename(TIME=X1)
dv=dff %>% filter(mod(dff$X7,dfn)==2) %>% select("X2") %>% rename(VOLUME=X2)
dtemp=dff %>% filter(mod(dff$X7,dfn)==2) %>% select("X3") %>% rename(TEMPERATURE=X3)
dpress=dff %>% filter(mod(dff$X7,dfn)==2) %>% select("X4") %>% rename(PRESSURE=X4)
de=dff %>% filter(mod(dff$X7,dfn)==2) %>% select("X6") %>% rename(ENERGY=X6)
drho=dff %>% filter(mod(dff$X7,dfn)==3) %>% select("X1") %>% rename(DENSITY=X1)

df=data.frame(dt,dv,dtemp,dpress,de,drho)

#remove(nxx,nx,n1,n2,n3,n4,n5,dt,dv,dtemp,dpress,de,drho,dff,dfn)

# seção volume versus tempo

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("Linear regression of Volume versus time","\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("\n")

ajuste=lm(VOLUME ~ TIME,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=TIME,y=VOLUME))+
  geom_line(color="black")+
  geom_smooth(method = "lm",color="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="Time (ps)",y=expression(paste("Volume ",(ring(A)^3))))+
  ggsave("volume.png",width = 6,height = 4)

summary(ajuste)

# iniciando nova seção

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("Linear regression of Temperature versus time","\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("\n")

ajuste=lm(TEMPERATURE ~ TIME,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=TIME,y=TEMPERATURE))+
  geom_line(color="black")+
  geom_smooth(method = "lm",color="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="Time (ps)",y=expression(paste("Temperature (K)")))+
  ggsave("temperature.png",width = 6,height = 4)

summary(ajuste)

# iniciando nova seção

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("Linear regression of Pressure versus time","\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("\n")

ajuste=lm(PRESSURE ~ TIME,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=TIME,y=as.double(PRESSURE)))+
  geom_line(color="black")+
  geom_smooth(method = "lm",color="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="Time (ps)",y=expression(paste("Pressure (atm)")))+
  ggsave("pressure.png",width = 6,height = 4)

summary(ajuste)

# iniciando nova seção

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("Linear regression of Energy versus time","\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("\n")

ajuste=lm(ENERGY ~ TIME,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=TIME,y=as.double(ENERGY)))+
  geom_line(color="black")+
  geom_smooth(method = "lm",color="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="Time (ps)",y=expression(paste("Energy (kcal/mol)")))+
  ggsave("energy.png",width = 6,height = 4)

summary(ajuste)

# iniciando nova seção

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("Linear regression of Density versus time","\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("\n")

ajuste=lm(DENSITY ~ TIME,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=TIME,y=as.double(DENSITY)))+
  geom_line(color="black")+
  geom_smooth(method = "lm",color="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="Time (ps)",y=expression(paste("Density ",(kg/m^3))))+
  ggsave("density.png",width = 6,height = 4)

summary(ajuste)

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n")
