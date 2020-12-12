# importanto e transformando dataframe de dados termodinamicos

head=readr::read_table(file = "HICOLM.thermodynamics",col_names = FALSE,n_max = 1,skip = 3)

head$X1=NULL
head$X4=as.integer(head$X4)

thermo=readr::read_table2(file="HICOLM.thermodynamics",col_names = TRUE,skip = 5)

thermo = thermo %>% gather(volume,temperature,pressure,energy,density,key="data",value="Value")

# importanto e transformando dataframe de dados atomicos

atoms=readr::read_table2(file="HICOLM.atoms",col_names = TRUE,skip = 5)

atoms$step=as.integer(atoms$step)
atoms$site=as.integer(atoms$site)
atoms$molecule=as.integer(atoms$molecule)
atoms$Z=as.integer(atoms$Z)