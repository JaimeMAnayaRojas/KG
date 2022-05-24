####  
library(data.table)
library(brms)
rm(list=ls(all=TRUE))



###########################################################################################################
# First get the data
getwd()
#setwd("~/Dropbox/Jaime M/Projects_JM/FSU/Pool_manipulation/KG_git/")
source("R/Functions.R")

center = 18
# Guppy data cleaning -----------------------------------------------------

Gdata <- read.csv("data/GuppyIPM.csv")
Kdata <- read.csv("data/KillifishIPM.csv")

Gdata$z <- Gdata$SL1_mm - center
Gdata$z2 <- Gdata$SL1_mm^2 - center^2
Gdata$z1 <- Gdata$SL2_mm 

# Make sure the non reproductive have zeros
Gdata$Recr[which(Gdata$Repr==0 & Gdata$surv ==1)] = 0



# Killifish data cleaning -------------------------------------------------

Kdata$z <- Kdata$SL1_mm - center
Kdata$z1 <- Kdata$SL2_mm 
Kdata$z2 <- Kdata$SL1_mm^2 - center^2


# Make sure the non reproductive have zeros
Kdata$Recr[which(Kdata$Repr==0 & Kdata$surv ==1)] = 0



# Normalize the covariates within streams -----------------

Gdata$FishBiom <- Gdata$BiomassG1 + Gdata$BiomassK1
Kdata$FishBiom <- Kdata$BiomassG1 + Kdata$BiomassK1

Gdata$Density <-( Gdata$BiomassG1 + Gdata$BiomassK1) / Gdata$area
Kdata$Density <- (Kdata$BiomassG1 + Kdata$BiomassK1) / Kdata$area


library("plyr")
names(Gdata)
vars = c("site", "Location", "Pool_1", "area", "FishBiom", "canopy", "Density")
df = rbind(Gdata[,vars], Kdata[,vars])
df = ddply(df,c('site','Location','Pool_1'),summarise, area=mean(area), Biomass=mean(FishBiom), canopy=mean(canopy), Density = mean(Density))

means = ddply(df,c('Location'),summarise, mean_area=mean(area), mean_Biomass=mean(Biomass), mean_canopy=mean(canopy), mean_Density= mean(Density))
sd = ddply(df,c('Location'),summarise, sd_area=sd(area), sd_Biomass=sd(Biomass), sd_canopy=sd(canopy), sd_Density= sd(Density))

Gdata = merge(Gdata, means, by.x = "Location")
Gdata = merge(Gdata, sd, by.x = "Location")

Kdata = merge(Kdata, means, by.x = "Location")
Kdata = merge(Kdata, sd, by.x = "Location")




Gdata$FishBiom <- (Gdata$FishBiom - Gdata$mean_Biomass) /  Gdata$sd_Biomass
Kdata$FishBiom <- (Kdata$FishBiom - Kdata$mean_Biomass) /  Kdata$sd_Biomass

Gdata$area <- (Gdata$area - Gdata$mean_area) /  Gdata$sd_area
Kdata$area <- (Kdata$area - Kdata$mean_area) /  Kdata$sd_area


Gdata$canopy <- (Gdata$canopy - Gdata$mean_canopy) /  Gdata$sd_canopy
Kdata$canopy <- (Kdata$canopy - Kdata$mean_canopy) /  Kdata$sd_canopy


Gdata$Density <- (Gdata$Density - Gdata$mean_Density) /  Gdata$sd_Density
Kdata$Density <- (Kdata$Density - Kdata$mean_Density) /  Kdata$sd_Density



# Re-code the random effects (drainage id) --------------------------------


Kdata$stream <- factor(Kdata$Location)
levels(Kdata$stream) <- 1:length(levels(Kdata$stream))
Kdata$stream <- as.numeric(Kdata$stream)

Gdata$stream <- factor(Gdata$Location)
levels(Gdata$stream) <- 1:length(levels(Gdata$stream))
Gdata$stream <- as.numeric(Gdata$stream)

Gdata$growth = log(Gdata$SL2_mm/Gdata$SL1_mm)
Kdata$growth = log(Kdata$SL2_mm/Kdata$SL1_mm)



# Model selection for Guppies ---------------------------------------------

# Survival ----------------------------------------------------------------
 Gdata <- Gdata[-which(Gdata$Sex2=="M"),]
# 
 Kdata <- Kdata[-which(Kdata$Sex2=="M"),]


Gdata$NK = factor(Gdata$NK)  


# Survival model ----------------------------------------------------------

S1 <- brm(surv ~ z*NK + FishBiom + Density + canopy + (1|stream), family = bernoulli(), Gdata,
         iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13), chains = 1)
summary(S1)

Survival_G = Model_selection(S1, name = "Survival", species = "Guppy")
summary(Survival_G)


post = data.frame(posterior_samples(Survival_G))
write.csv(post, "outputs/Post_survivalG.csv")

# Growth model ------------------------------------------------------------

G1 <- brm(z1 ~ z*NK + FishBiom + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

summary(G1)

Growth_G = Model_selection(G1)
Growth_G$Selection
Growth_G$Best_model


# Reproduction model ------------------------------------------------------------

R1 <- brm(Recr ~ z*NK + FishBiom + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

Repro_G = Model_selection(R1)


## predicted responses of patient 1 for new data

post = posterior_samples(Repro_G)

p_linkG = function(post, NK=0, size, center = center){
  z = size - center
  pp = post$b_Intercept + post$b_z *z  + post$b_NK1 * NK + post$`b_z:NK1`* (z*NK)
  return(pp)
  
}

size <-  seq(from= 5, to=30, length.out = 200)


# Predictions for Growth
post_G = posterior_samples(Growth_G)
p.KG_G= (sapply(1:length(size), function(i) p_linkG(post = post_G, NK = 0, size = size, center = 18)))
p.KG_G.mean <- apply(p.NK_G, 2, mean)
p.KG_G.PI <- apply(p.NK_G, 2, HPDI,prob = .95)

p.NK_G= exp(sapply(1:length(size), function(i) p_linkG(post = post, NK = 1, size = size, center = 18)))
p.NK_G.mean <- apply(p.NK_G, 2, mean)
p.NK_G.PI <- apply(p.NK_G, 2, HPDI,prob = .95)


# Predictions for Reproduction
p.KG_G= exp(sapply(1:length(size), function(i) p_linkG(post = post, NK = 0, size = size, center = 18)))
p.KG_G.mean <- apply(p.NK_G, 2, mean)
p.KG_G.PI <- apply(p.NK_G, 2, HPDI,prob = .95)

p.NK_G= exp(sapply(1:length(size), function(i) p_linkG(post = post, NK = 1, size = size, center = 18)))
p.NK_G.mean <- apply(p.NK_G, 2, mean)
p.NK_G.PI <- apply(p.NK_G, 2, HPDI,prob = .95)



LOS(p.NK_G- p.KG_G)


svg("plots/Recr_IntSlope.svg", width = 7, height = 3.5)
op<-par(mfrow=c(1,2), mar = c(4, 4 , 2, 1), oma = c(0.5, 1, 1, 0.5)) # c(bottom, left, top, right)
dens((post$b_Intercept + post$b_NK1), show.HPDI = T, main = paste("Intercept ", "(", LOS(post$b_NK1), "%", ")", sep = "" ), col= "red", ylim = c(0,1.8))
dens((post$b_Intercept), add = T, show.HPDI = T)
dens((post$b_z + post$`b_z:NK1`), show.HPDI = T,  col= "red", ylim = c(0,15), main = paste("Slope ", "(", LOS(post$`b_z:NK1`), "%", ")", sep = "" ))
dens((post$b_z), add = T, show.HPDI = T)
graphics.off()



svg("plots/Recr_IntSlope.svg", width = 7, height = 3.5)
op<-par(mfrow=c(1,2), mar = c(4, 4 , 2, 1), oma = c(0.5, 1, 1, 0.5)) # c(bottom, left, top, right)
dens((post$b_Intercept + post$b_NK1), show.HPDI = T, main = paste("Intercept ", "(", LOS(post$b_NK1), "%", ")", sep = "" ), col= "red", ylim = c(0,1.8))
dens((post$b_Intercept), add = T, show.HPDI = T)

dens((post$b_z + post$`b_z:NK1`), show.HPDI = T,  col= "red", ylim = c(0,15), main = paste("Slope ", "(", LOS(post$`b_z:NK1`), "%", ")", sep = "" ))
dens((post$b_z), add = T, show.HPDI = T)
graphics.off()
