# clean up workspace

rm(list = ls())

# close all figure windows created with x11()

graphics.off()

# load packages

library(tibble)
library(ggplot2)

setwd("/home/flaviano/MEGA/maria_eduarda/simulacoes/imidazol+H2O/simulacoes/MD/berendsen")

# carregando dataframe

dff=readr::read_table2(file="HICOLM.md",col_names =
                         c("V1","V2","V3","V4","V5","V6"),
                       skip_empty_rows = TRUE)

nxx=as.integer(dff[1,4])

dfn=as.integer(3*nxx+9)

df=dff[seq(4,nrow(dff),dfn),]
df[,7]=dff[seq(5,nrow(dff),dfn),1]

names(df)[7]="V7"

# seção volume versus tempo

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("Linear regression of Volume versus time","\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","\n")
cat("\n")

ajuste=lm(V2 ~ V1,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=V1,y=V2))+
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

ajuste=lm(V3 ~ V1,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=V1,y=V3))+
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

ajuste=lm(V4 ~ V1,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=V1,y=as.double(V4)))+
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

ajuste=lm(V6 ~ V1,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=V1,y=as.double(V6)))+
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

ajuste=lm(V7 ~ V1,data=df)

# imprimindo variaveis termodinamicas

  df %>% ggplot(aes(x=V1,y=as.double(V7)))+
  geom_line(color="black")+
  geom_smooth(method = "lm",color="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="Time (ps)",y=expression(paste("Density ",(kg/m^3))))+
  ggsave("density.png",width = 6,height = 4)

summary(ajuste)

cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n")
cat("\n")

head(df)
cat("\n")
