#carregando pacotes

library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)
library(ggplot2)

setwd("/home/flaviano/Dropbox/Projeto/simulacoes/H2O_MD")
#setwd("/home/flaviano/Documentos/maria_eduarda")

# carregando dataframe

dff=readr::read_table2(file="HICOLM.md",col_names = c("X1","X2","X3","X4","X5","X6"),skip_empty_rows = TRUE)
#dff=dff[rowSums(is.na(dff))!=ncol(dff),]

nx=list(dff[1,1:5])
nxx=as.integer(dff[1,4])

dfn=as.integer(3*nxx+7)

df1=dff[seq(4,nrow(dff),dfn),]
df2=dff[seq(5,nrow(dff),dfn),]

dt=df1 %>% select("X1") %>% rename(TIME=X1)
dv=df1 %>% select("X2") %>% rename(VOLUME=X2)
dtemp=df1 %>% select("X3") %>% rename(TEMPERATURE=X3)
dpress=df1 %>% select("X4") %>% rename(PRESSURE=X4)
de=df1 %>% select("X6") %>% rename(ENERGY=X6)
drho=df2 %>% select("X1") %>% rename(DENSITY=X1)

df=data.frame(dt,dv,dtemp,dpress,de,drho)

remove(nxx,nx,dt,dv,dtemp,dpress,de,drho,dfn,df1,df2,dff)

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
