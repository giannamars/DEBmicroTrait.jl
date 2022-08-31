library("tidyverse")
library("lme4")
library("MuMIn")
library("FSA")
library("PMCMRplus")
library("randomForest")
library("gbm")
set.seed(1)

## Batch model BGE predictions
BGE.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_BGE.csv") %>% filter(!is.nan(BGE))

##########################################################################################
## Taxonomic and Resource Variance Partitioning

# Model 1: Isolate Identity
fm1 <- lmer(BGE ~ as.factor(isolate) + (1|monomer/ontology), data = BGE.f)
summary(fm1)$AIC
r.squaredGLMM(fm1)
qqnorm(resid(fm1))

# Model 2: Class Order
fm2 <- lmer(BGE ~ as.factor(class) + (1|monomer/ontology), data = BGE.f)
summary(fm2)$AIC
r.squaredGLMM(fm2)
qqnorm(resid(fm2))

# Model 3: Phylum Order
fm3 <- lmer(BGE ~ as.factor(phylum) + (1|monomer), data = BGE.f)
summary(fm3)$AIC
r.squaredGLMM(fm3)
qqnorm(resid(fm3))

# Model 4: Monomer Identity
fm4 <- lmer(BGE ~ as.factor(monomer) + (1|isolate/class/phylum), data = BGE.f)
summary(fm4)$AIC
r.squaredGLMM(fm4)
qqnorm(resid(fm4))

# Model 5: Monomer Class
fm5 <- lmer(BGE ~ as.factor(ontology) + (1|isolate/class/phylum), data = BGE.f)
summary(fm5)$AIC
r.squaredGLMM(fm5)
qqnorm(resid(fm5))

##########################################################################################
## Model selection

BGE.f.high <- BGE.f[BGE.f$rgrowth > 0.041, ]
BGE.f.low <- BGE.f[BGE.f$rgrowth < 0.0407, ]

# Growth rate
# low
lm.f = lm(rgrowth ~ rrn + rrn:genomesize, data = BGE.f.low)
summary(lm.f)
AIC(lm.f)

lm.f = lm(rgrowth ~ genomesize + rrn:genomesize, data = BGE.f.low)
summary(lm.f)
AIC(lm.f)

lm.f = lm(rgrowth ~ rrn*genomesize, data = BGE.f.low)
summary(lm.f)
AIC(lm.f)

# high
lm.f = lm(rgrowth ~ rrn + rrn:genomesize, data = BGE.f.high)
summary(lm.f)
AIC(lm.f)

lm.f = lm(rgrowth ~ genomesize + rrn:genomesize, data = BGE.f.high)
summary(lm.f)
AIC(lm.f)

lm.f = lm(rgrowth ~ rrn*genomesize, data = BGE.f.high)
summary(lm.f)
AIC(lm.f)

# CUE
# low
lm.f = lm(BGE ~ rrn + rrn:genomesize, data = BGE.f.low)
summary(lm.f)
AIC(lm.f)

lm.f = lm(BGE ~ genomesize + rrn:genomesize, data = BGE.f.low)
summary(lm.f)
AIC(lm.f)

lm.f = lm(BGE ~ rrn*genomesize, data = BGE.f.low)
summary(lm.f)
AIC(lm.f)

# high
lm.f = lm(BGE ~ rrn + rrn:genomesize, data = BGE.f.high)
summary(lm.f)
AIC(lm.f)

lm.f = lm(BGE ~ genomesize + rrn:genomesize, data = BGE.f.high)
summary(lm.f)
AIC(lm.f)

lm.f = lm(BGE ~ rrn*genomesize, data = BGE.f.high)
summary(lm.f)
AIC(lm.f)

##########################################################################################
## BGE-growth Regressions
BGE.f.grouped <- BGE.f %>% group_by(response, ontology) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth), BP_med = median(BP), BR_med = median(BR))
# low growth regime
growthLOW <- lm(rgrowth_med ~ BGE_med, data = BGE.f.grouped[BGE.f.grouped$rgrowth_med < 0.0407, ])
summary(growthLOW)
# high growth regine
growthHIGH <- lm(rgrowth_med ~ BGE_med, data = BGE.f.grouped[BGE.f.grouped$rgrowth_med > 0.041, ])
summary(growthHIGH)

# BP-BR scaling
ttest <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  pstat <- 2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
  return(list = c(tstat, reg$df.residual, pstat))
}

BGE.f$Group <- NA
for (i in 1:dim(BGE.f)[1])
  if (BGE.f$rgrowth[i] < 0.0407){
    BGE.f$Group[i] <- "A"
  } else {
    BGE.f$Group[i] <- "B"
  }

modhigh <- lm(log(BP) ~ log(BR), data =  BGE.f[BGE.f$Group == "B",])
summary(modhigh)
ttest(modhigh,2,1)

modlow <- lm(log(BP) ~ log(BR), data =  BGE.f[BGE.f$Group == "A",])
summary(modlow)
ttest(modlow,2,1)

#########################################################################################
## Random forest analysis
# CUE
df_train <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_train.csv")
df_train$yield <- 1- df_train$yield
df_test <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_test.csv")
df_test$yield <- 1- df_test$yield
# log transform data
logdf_train <- log(df_train) %>% filter(!is.nan(BGE))
logdf_test<- log(df_test) %>% filter(!is.nan(BGE))
# boost
boost.bge = gbm(BGE~., data=logdf_train, distribution="gaussian", n.trees=10000, interaction.depth = 8, shrinkage = 0.001)
summary(boost.bge)

#########################################################################################
## Flux variance
# high
pchigh.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_fluxes_high.csv")
kruskal.test(pc1 ~ response, data = pchigh.f)
dunnTest(pc1 ~ response, data = pchigh.f, method = "bh")
pchigh.f$response = as.factor(pchigh.f$response)
kwAllPairsConoverTest(pc1 ~ response, data = pchigh.f)

kruskal.test(pc2 ~ response, data = pchigh.f)
dunnTest(pc2 ~ response, data = pchigh.f, method = "bh")
pchigh.f$response = as.factor(pchigh.f$response)
kwAllPairsConoverTest(pc2 ~ response, data = pchigh.f)

# low
pclow.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_fluxes_low.csv")
kruskal.test(pc1 ~ response, data = pclow.f)
dunnTest(pc1 ~ response, data = pclow.f, method = "bh")
pclow.f$response = as.factor(pclow.f$response)
kwAllPairsConoverTest(pc1 ~ response, data = pclow.f)

kruskal.test(pc2 ~ response, data = pclow.f)
dunnTest(pc2 ~ response, data = pclow.f, method = "bh")
pclow.f$response = as.factor(pclow.f$response)
kwAllPairsConoverTest(pc2 ~ response, data = pclow.f)


pcflux.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_fluxes_low.csv")
pcflux.f.num <- pcflux.f[,2:15]
pcflux.f.num["rMco2_frac"] = (pcflux.f.num$rMco2)/(pcflux.f.num$rMco2+ pcflux.f.num$rGco2 + pcflux.f.num$rXco2)
pcflux.f.num["turnover"] = pcflux.f.num$jE + pcflux.f.num$jV
pcflux.f.pn <- pcflux.f.num %>% filter(response=="positive" | response =="negative")
pcflux.pn.stats <- pcflux.f.pn %>% group_by(response) %>% summarise_each(funs(median, sd))
pval <- sapply(2:16, function(x) kruskal.test(pcflux.f.pn[,x], pcflux.f.pn$response)$p.value)

#########################################################################################
# Levin's index
levins <- function(p_xi = ""){
  p = 0
  for (i in p_xi){
    p = p + i^2
  }
  nb = 1 / (length(p_xi) * p)
  return(nb)
}

# Niche breadth
Levin.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_levin.csv") %>%
  replace(is.na(.), 0.0)
Levin.f.std <- Levin.f[,1:83]/(apply(Levin.f[,1:83], 1, sum))

Levin.f.grouped <- BGE.f %>% group_by(isolate) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth))
Levin.f.grouped$levins <- levins(Levin.f.std)
Levin.f.grouped$response <- Levin.f$response
plot(Levin.f.grouped$levins, Levin.f.grouped$BGE_med)
plot(density(Levin.f.grouped$levins))
