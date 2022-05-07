
# MODEL COMPARASION FOR GROWTH

library("brms")
library(data.table)
rm(list=ls(all=TRUE))

###########################################################################################################
# First get the data
getwd()
setwd("~/Dropbox/Jaime M/Projects_JM/FSU/Pool_manipulation/KG_git/")


center = 18
# Guppy data cleaning -----------------------------------------------------

Gdata <- read.csv("data/GuppyIPM.csv")
Kdata <- read.csv("data/KillifishIPM.csv")

Gdata$z <- Gdata$SL1_mm - center
Gdata$z1 <- Gdata$SL2_mm 
Gdata <- Gdata[-which(Gdata$Sex2=="M"),]

# Make sure the non reproductive have zeros
Gdata$Recr[which(Gdata$Repr==0 & Gdata$surv ==1)] = 0



# Killifish data cleaning -------------------------------------------------

Kdata$z <- Kdata$SL1_mm - center
Kdata$z1 <- Kdata$SL2_mm 


Kdata <- Kdata[-which(Kdata$Sex2=="M"),]

# Make sure the non reproductive have zeros
Kdata$Recr[which(Kdata$Repr==0 & Kdata$surv ==1)] = 0



# Normalize the covariates (killifish biomass and canopy) -----------------

Gdata$FishBiom <- Gdata$BiomassG1 + Gdata$BiomassK1
Kdata$FishBiom <- Kdata$BiomassG1 + Kdata$BiomassK1
meanC = mean(unique(c(Gdata$FishBiom, Kdata$FishBiom ) ))
sdC = sd(unique(c(Gdata$FishBiom, Kdata$FishBiom ) ))
Gdata$FishBiom <- (Gdata$FishBiom - meanC) /  sdC
Kdata$FishBiom <- (Kdata$FishBiom - meanC) /  sdC


# fish density biomass/Area

Gdata$FishDen <- Gdata$FishBiom  / Gdata$area
Kdata$FishDen <- Kdata$FishBiom  / Kdata$area

meanC = mean(unique(c(Gdata$FishDen, Kdata$FishDen ) ))
sdC = sd(unique(c(Gdata$FishDen, Kdata$FishDen ) ))
Gdata$FishDen <- (Gdata$FishDen - meanC) /  sdC
Kdata$FishDen <- (Kdata$FishDen - meanC) /  sdC

#


meanC = mean(unique(c(Gdata$canopy, Kdata$canopy) ))
sdC = sd(unique(c(Gdata$canopy, Kdata$canopy) ))
Gdata$canopy <- (Gdata$canopy - meanC) /  sdC
Kdata$canopy <- (Kdata$canopy - meanC) /  sdC




# Re-code the random effects (drainage id) --------------------------------


Kdata$stream <- factor(Kdata$Location)
levels(Kdata$stream) <- 1:length(levels(Kdata$stream))
Kdata$stream <- as.numeric(Kdata$stream)

Gdata$stream <- factor(Gdata$Location)
levels(Gdata$stream) <- 1:length(levels(Gdata$stream))
Gdata$stream <- as.numeric(Gdata$stream)



# Compare growh models ----------------------------------------------------

surv = bf(surv ~ z*NK + FishBiom + FishDen + canopy, family = bernoulli())
growth = bf(z1 ~ z*NK + FishBiom + FishDen + canopy, family = gaussian())


ModG1 = brm(z1 ~ z*NK + FishBiom + FishDen + canopy, Gdata)
ModG1 = add_criterion(ModG1, c("waic", "loo"))

ModG2 = brm(z1 ~ z*NK +  FishBiom + canopy, Gdata)
ModG2 = add_criterion(ModG2, c("waic", "loo"))

ModG3 = brm(z1 ~ z*NK +  FishDen + canopy, Gdata)
ModG3 = add_criterion(ModG3, c("waic", "loo"))


ModG4 = brm(z1 ~ z*NK +  canopy, Gdata)
ModG4 = add_criterion(ModG4, c("waic", "loo"))

loo_compare(ModG1, ModG2, ModG3, ModG4, criterion = "loo")



loo_compare( ModG2, ModG3, criterion = "loo")

loo_model_weights(ModG1, ModG2, ModG3, ModG4, criterion = "loo")

summary(ModG1)
summary(ModG2)

# pp_check(ModG1, nsamples = 100)
# pp_check(ModG2, nsamples = 100)


ModK1 = brm(z1 ~ z*NG + FishBiom + FishDen + canopy, Kdata)
ModK1 = add_criterion(ModK1, c("waic", "loo"))


ModK2 = brm(z1 ~ z*NG +  FishBiom + canopy, Kdata)
ModK2 = add_criterion(ModK2, c("waic", "loo"))

ModK3 = brm(z1 ~ z*NG +  FishDen + canopy, Kdata)
ModK3 = add_criterion(ModK3, c("waic", "loo"))


ModK4 = brm(z1 ~ z*NG +  canopy, Kdata)
ModK4 = add_criterion(ModK4, c("waic", "loo"))

loo_compare(ModK1, ModK2, ModK3, ModK4, criterion = "loo")
