library("tidyverse")
library("broom")

## Rhizosphere model predictions
rhizo.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/rhizosphere_model_2006.csv")  

rhizo.f %>%
  nest(data = -sample) %>%
  mutate(
    fit = map(data, ~ lm(measurement ~ relabundance, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance)
  ) %>%
  unnest(glanced)

ggplot(data = rhizo.f, aes(x=relabundance,y=measurement))+
  geom_smooth(aes(color = factor(sample)), method = "lm", se = F)


## Avena benchmark
avena.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/avena/TOC_data_Hopland_Avena.csv")
avena.f$diff <- avena.f$TOCsim - avena.f$TOCmM
model <- lm(TOCmM ~ TOCsim, data = avena.f)
summary(model)
modeldiff <- lm(diff ~ TOCsim, data = avena.f)
summary(modeldiff)
