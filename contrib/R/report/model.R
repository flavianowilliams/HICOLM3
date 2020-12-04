# calculo da forca maxima

value_max = atoms %>%
  mutate(f=sqrt(fx**2+fy**2+fz**2)) %>%
  group_by(time) %>%
  summarize(fmax=max(f))
