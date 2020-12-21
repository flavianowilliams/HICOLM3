# calculo da forca maxima

value_max = atoms %>%
  mutate(f=sqrt(fx**2+fy**2+fz**2)) %>%
  group_by(time) %>%
  summarize(fmax=max(f))

# Estatistica das variaveis termodinamicas

data_mean = thermo %>%
  group_by(data) %>%
  summarise(mean(value))
data_sd =  thermo %>%
  group_by(data) %>%
  summarise(sd(value))
data_min =  thermo %>%
  group_by(data) %>%
  summarise(min(value))
data_max =  thermo %>%
  group_by(data) %>%
  summarise(max(value))
data_mean = data_mean %>%
  inner_join(data_sd, by="data") %>%
  inner_join(data_min, by="data") %>%
  inner_join(data_max, by="data") %>%
  rename("Mean"="mean(value)","Standard deviation"="sd(value)",
         "Minimum value"="min(value)","Maximum value"="max(value)")

# limpando memoria

rm(data_sd,data_min,data_max)
