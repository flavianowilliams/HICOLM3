# importanto e transformando dataframe de dados termodinamicos

#setwd("/home/flaviano/Documentos/GitHub/HICOLM/examples/H2O")
#library(tidyverse)

head=readr::read_table2(file = "HICOLM.out",col_names = c("X1","X2","X3"),n_max = 10,skip = 42)

head$X1=NULL
head$X3[4]=0
head$X3[5]=0
head$X3[9]=0
head=t(head)
head=as_tibble(head, .name_repair = ~ c("ensemble","thermostat","barostat","temperature",
                                        "pressure","nstep","nrelax","reuse","timestep","r_cutoff"))

head$ensemble=as.character(head$ensemble)
head$thermostat=as.numeric(head$thermostat)
head$barostat=as.numeric(head$barostat)
head$temperature=as.numeric(head$temperature)
head$pressure=as.numeric(head$pressure)
head$nstep=as.integer(head$nstep)
head$nrelax=as.integer(head$nrelax)
head$reuse=as.integer(head$reuse)
head$timestep=as.numeric(head$timestep)
head$r_cutoff=as.numeric(head$r_cutoff)

system=readr::read_table2(file = "HICOLM.out",col_names = c("X1","X2","X3"),n_max = 1,skip = 67)

system$X1=NULL

thermo=readr::read_csv(file="thermodynamics.csv",col_names = TRUE)

thermo = thermo %>% gather(volume,temperature,pressure,energy,density,key="data",value="Value")

# importanto e transformando dataframe de dados atomicos

atoms=readr::read_csv(file="atoms.csv",col_names = TRUE)

atoms$step=as.integer(atoms$step)
atoms$site=as.integer(atoms$site)
atoms$molecule=as.integer(atoms$molecule)
atoms$Z=as.integer(atoms$Z)
