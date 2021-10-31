####  
library(rethinking)
library(data.table)
rm(list=ls(all=TRUE))
###########################################################################################################
# First get the data
getwd()
setwd("~/Dropbox/Projects_JM/FSU/Pool_manipulation/KG_git/")

Gdata <- read.csv("data/GuppyIMP.csv")
Kdata <- read.csv("data/KillifishIPM.csv")

Gdata$z <- log(Gdata$SL1_mm) - log(18)
Gdata$z1 <- log(Gdata$SL2_mm) 
Gdata <- Gdata[-which(Gdata$Sex2=="M"),]
Gdata$Recr[which(Gdata$Repr==0 & Gdata$surv ==1)] = 0
Gdata$area <- Gdata$area - mean(unique(Gdata$area))
Gdata$canopy <- Gdata$canopy - mean(unique(Gdata$canopy))


Kdata$z <- log(Kdata$SL1_mm) - log(18)
Kdata$z1 <- log(Kdata$SL2_mm) 



Kdata <- Kdata[-which(Kdata$Sex2=="M"),]
Kdata$Recr[which(Kdata$Repr==0 & Kdata$surv ==1)] = 0
Kdata$area <- Kdata$area - mean(unique(Kdata$area))
Kdata$canopy <- Kdata$canopy - mean(unique(Kdata$canopy))




Kdata$stream <- factor(Kdata$Location)
levels(Kdata$stream) <- 1:length(levels(Kdata$stream))
Kdata$stream <- as.numeric(Kdata$stream)

Gdata$stream <- factor(Gdata$Location)
levels(Gdata$stream) <- 1:length(levels(Gdata$stream))
Gdata$stream <- as.numeric(Gdata$stream)



data_stan = list(
  # variables size 
  N_survG = length(Gdata$surv),
  N_growG = length(which(Gdata$surv ==1)),
  N_repG = length(which(Gdata$surv ==1)),
  N_recrG = length(which(Gdata$Repr ==1)),
  
  N_survK = length(Kdata$surv),
  N_growK = length(which(Kdata$surv ==1)),
  N_repK = length(which(Kdata$surv ==1)),
  N_recrK = length(which(Kdata$Repr ==1)),
  
  N_stream =  length(unique(Gdata$stream)),
  
  #Response variables
  Surv_G = Gdata$surv,
  z_survG = Gdata$z,
  NK_survG = Gdata$NK,
  stream_survG = Gdata$stream,
  area_survG = Gdata$area,
  canopy_survG = Gdata$canopy,
  
  z1_G   = Gdata$z1[which(Gdata$surv ==1)],
  z_growG = Gdata$z[which(Gdata$surv ==1)],
  NK_growG = Gdata$NK[which(Gdata$surv ==1)],
  stream_growG = Gdata$stream[which(Gdata$surv ==1)],
  area_growG = Gdata$area[which(Gdata$surv ==1)],
  canopy_growG = Gdata$canopy[which(Gdata$surv ==1)],
  
  
  Repr_G = Gdata$Repr[which(Gdata$surv ==1)],
  z_repG = Gdata$z[which(Gdata$surv ==1)],
  NK_repG = Gdata$NK[which(Gdata$surv ==1)],
  stream_repG = Gdata$stream[which(Gdata$surv ==1)],
  area_repG = Gdata$area[which(Gdata$surv ==1)],
  canopy_repG = Gdata$canopy[which(Gdata$surv ==1)],
  
  
  Recr_G = Gdata$Recr[which(Gdata$Repr ==1)],
  z_recrG = Gdata$z[which(Gdata$Repr ==1)],
  NK_recrG = Gdata$NK[which(Gdata$Repr ==1)],
  stream_recrG = Gdata$stream[which(Gdata$Repr ==1)],
  area_recrG = Gdata$area[which(Gdata$Repr ==1)],
  canopy_recrG = Gdata$canopy[which(Gdata$Repr ==1)],
  
  
  Surv_K = Kdata$surv,
  z_survK = Kdata$z,
  NG_survK = Kdata$NG,
  stream_survK = Kdata$stream,
  area_survK = Kdata$area,
  canopy_survK = Kdata$canopy,
  
  
  z1_K   = Kdata$z1[which(Kdata$surv ==1)],
  z_growK = Kdata$z[which(Kdata$surv ==1)],
  NG_growK = Kdata$NG[which(Kdata$surv ==1)],
  stream_growK = Kdata$stream[which(Kdata$surv ==1)],
  area_growK = Kdata$area[which(Kdata$surv ==1)],
  canopy_growK = Kdata$canopy[which(Kdata$surv ==1)],
  z_gK = log(Kdata$z[which(Kdata$surv ==1)] + 18) - log(18),
  z1_gK = log(Kdata$z1[which(Kdata$surv ==1)]) - log(18),
  
  
  Repr_K = Kdata$Repr[which(Kdata$surv ==1)],
  z_repK = Kdata$z[which(Kdata$surv ==1)],
  NG_repK = Kdata$NG[which(Kdata$surv ==1)],
  stream_repK = Kdata$stream[which(Kdata$surv ==1)],
  area_repK = Kdata$area[which(Kdata$surv ==1)],
  canopy_repK = Kdata$canopy[which(Kdata$surv ==1)],
  
  
  
  Recr_K = Kdata$Recr[which(Kdata$Repr ==1)],
  z_recrK = Kdata$z[which(Kdata$Repr ==1)],
  NG_recrK = Kdata$NG[which(Kdata$Repr ==1)],
  stream_recrK = Kdata$stream[which(Kdata$Repr ==1)],
  area_recrK = Kdata$area[which(Kdata$Repr ==1)],
  canopy_recrK = Kdata$canopy[which(Kdata$Repr ==1)] 

  
  
)
# 

# modG = stan("R/models/pool_mod.stan", data = data_stan, cores = 4, chains = 4, 
#                 iter = 6000, warmup = 4500, control = list(adapt_delta = 0.92, max_treedepth = 12))

 
 modG = stan("R/models/LOG_mod.stan", data = data_stan, cores = 4, chains = 4, 
             iter = 6000, warmup = 4500, control = list(adapt_delta = 0.92, max_treedepth = 12))


# saveRDS(modG, "Model_priors_Bassar_2017.RDS")

# modG <- readRDS("Model.RDS")
#modG <- readRDS("Model_priors_Bassar_2017.RDS")
# 
 
 LOS <- function(x=NULL){
   
   out =   (length(which(x > 0)) / length(x))*100
   
   return(round(out, 3))
   
 }
 
 
 
precis(modG, digits = 5, prob = .95, depth = 2)




#traceplot_ulam(modG)

#

pos <- extract(modG)
write.csv(as.data.frame(pos), "Posteriors.csv")
write.csv(as.data.frame(pos), "PosteriorsLOG_K.csv")

PP = as.vector(apply(as.data.frame(pos), 2, LOS))

model.summary <- as.data.frame(precis(modG, prob = .95, digits = 3, depth = 2))

model.summary$LOS_l = PP[-65]
model.summary$LOS_u = 100- PP[-65]
model.summary[model.summary$`2.5%` > 0,]
model.summary[model.summary$`97.5%` < 0,] 
model.summary[model.summary$`2.5%` > 0,] 

model.summary

model.summary[model.summary$LOS_l > 90,] 




#rm(list=ls(all=TRUE))
# 
# head(as.data.frame(pos))
# 
# 
# model.summary
# 
# 100-LOS(pos$b_NK_survG)
# 
# LOS(pos$b_zNK_survG)
# 
# LOS(pos$b_z_survG + pos$b_zNK_survG )
# # Plots for each vital rate
# 
# # Make survival plot
# 
# surv_link <- function(post, size, G = 1, NK, NG, center){
#   
#   z = size - center
#   
#   if(G == 1){
#     p =  inv_logit(with(post, Intercept_survG + b_NK_survG * NK + b_z_survG * z + b_zNK_survG* NK*z))    
#   }else{
#     p =  inv_logit(with(post,Intercept_survK + b_NG_survK * NG + b_z_survK * z + b_zNG_survK* NG*z))
#   }
#   
#   return(p)
# }
# 
# 
# 
# grow_link <-function(post, size, G = 1, NK, NG, center){
#   
#   z = size - center
#   
#   if(G == 1){
#     mu = with(post,Intercept_growG + b_NK_growG * NK + b_z_growG * z + b_zNK_growG * NK * z ) 
#   }else{
#     mu = with(post,Intercept_growK + b_NG_growK * NG + b_z_growK * z + b_zNG_growK * NG * z)
#   }
#   
#   
#   return(mu-size)
# }
# 
# 
# 
# 
# recr_link <-function(post, size, G = 1, NK, NG, center){
#   
#   z = size - center
#   
#   if(G == 1){
#     lambda = with(post,Intercept_recrG + b_NK_recrG * NK + b_z_recrG * z + b_zNK_recrG * NK * z )
#   }else{
#     lambda = with(post, Intercept_recrK + b_NG_recrK * NG + b_z_recrK * z + b_zNG_recrK * NG * z )
#   }
#   return(exp(lambda))
# }
# 


# ##
# 
# #------------------------ IPMs
#  source("R/IPM_guppy.R") # Guppies
# 
# G_lamda.est <- read.csv("G_lamda.est.csv")[,-1]
# Gsurv.mat <- read.csv("Gsurv.mat.csv")[,-1]
# Ggrow.mat <- read.csv("Ggrow.mat.csv")[,-1]
# Gfec.mat <- read.csv( "Gfec.matc.csv")[,-1]
# Grcz.mat <- read.csv( "Grcz.mat.csv")[,-1]
# 
# names(G_lamda.est)
# 
# mins = abs(apply(Gsurv.mat + Ggrow.mat + Gfec.mat +  Grcz.mat, 1, min))
# mats = list(surv = Gsurv.mat, grow=  Ggrow.mat, fec= Gfec.mat, recz=   Grcz.mat)
# mats$surv_p= G_lamda.est$surv.p
# mats$grow_p= G_lamda.est$grow.p
# mats$fec_p= G_lamda.est$fec.p
# mats$rcz_p= G_lamda.est$rcz.p
# 
# i=1
# for(i in 1:dim(Gsurv.mat)[1]){
#   d  = sum(abs(mats$surv[i,]) + abs(mats$grow[i,]) + abs(mats$fec[i,]) + abs(mats$recz[i,]))
#   mats$surv_p[i]= sum(mats$surv[i,]) / d
#   mats$grow_p[i]= sum(mats$grow[i,]) /d
#   mats$fec_p[i]= sum(mats$fec[i,]) /d
#   mats$rcz_p[i]= sum(mats$recz[i,]) /d
#   
#   
# }
# 
# mean(mats$surv_p)
# HPDI(mats$surv_p, prob = .95)
# 
# mean(mats$grow_p)
# HPDI(mats$grow_p, prob = .95)
# 
# mean(mats$fec_p)
# HPDI(mats$fec_p, prob = .95)
# 
# 
# G_lamda.est[1,]
# 
# mats$lam_fit =  mats$surv + mats$grow + mats$fec + mats$recz
# 
# apply(mats$lam_fit, 1, sum)
# 
# 
# 
# 
# Gsum = as.data.frame(precis(G_lamda.est, prob = .95, digits = 3))
# Gsum$ppL <- apply(G_lamda.est, 2, LOS)
# Gsum$ppU <- 100- apply(G_lamda.est, 2, LOS)
# Gsum
# 
# #source("R/IPM_killifish.R") # Killifish
# 
# K_lamda.est <- read.csv("K_lamda.est.csv")[,-1]
# Ksurv.mat <- read.csv("Ksurv.mat.csv")[,-1]
# Kgrow.mat <- read.csv("Kgrow.mat.csv")[,-1]
# Kfec.mat <- read.csv( "Kfec.matc.csv")[,-1]
# Krcz.mat <- read.csv( "Krcz.mat.csv")[,-1]
# 
# 
# 
# Ksum = as.data.frame(precis(K_lamda.est, prob = .95, digits = 3))
# Ksum$ppL <- apply(K_lamda.est, 2, LOS)
# Ksum$ppU <- 100- apply(K_lamda.est, 2, LOS)
# Ksum
# 
# ###
# 
# graphics.off()
# plot(G_lamda.est$lamda.dif ~ G_lamda.est$sum.comp.lamda)
# abline(0,1)
# plot(K_lamda.est$lamda.dif ~ K_lamda.est$sum.comp.lamda)
# abline(0,1)
# 
# 
# toPlot<- Gsum[c("Sum.lam.effect", "lam.surv", "lam.growth", "lam.fec", "lam.rcz"),]
# toPlot$x <- 1:5
# 
# 
# # size specific LTRE
# 
# nBigMatrix <- 100
# # Make the meshpoints
# min.size <- (2) #min(IPMdata$z1, na.rm = T) - sd(IPMdata$z1, na.rm = T)*0.5
# max.size <- (40) #max(IPMdata$z1, na.rm = T) + sd(IPMdata$z1, na.rm = T)*0.5
# 
# U=max.size
# L=min.size
# m = nBigMatrix
# h <- (U - L)/m
# meshpts <- z1 <- z <- L + ((1:m) - 1/2) * h
# size.cen <- (18)
# mean <- apply(Gsurv.mat,2, mean)
# ci <- apply(Gsurv.mat,2, HPDI, prob=.95)
# meshpts[t(ci)[,1] >0]
# meshpts[t(ci)[,2] < 0]
# 
# library(latex2exp)
# plot(mean ~ meshpts, ylab = TeX("Fitness effects ($\\Delta \\lambda$)"), pch="", ylim = c(-0.15,0.05), xlim = c(5,30))
# abline(h=0, lty=2)
# lines(mean ~ meshpts, lwd=1.5, col = "red")
# shade(ci, meshpts, col = col.alpha('red', 0.2))
# 
# mean <- apply(Ggrow.mat,2, mean)
# ci <- apply(Ggrow.mat,2, HPDI, prob=.95)
# range(meshpts[which(apply(Ggrow.mat, 2, LOS) > 95)])
# 
# lines(mean ~ meshpts, lwd=1.5, col = "blue")
# shade(ci, meshpts, col = col.alpha('blue', 0.2))
# 
# 
# plot(mean ~ meshpts, ylab = TeX("Fitness effects ($\\Delta \\lambda$)"), pch="", ylim = c(-0.1,0.02), xlim = c(5,30))
# abline(h=0, lty=2)
# 
# ## Killifish
# 
# nBigMatrix <- 100
# 
# # Make the meshpoints
# min.size <- (2) 
# max.size <- (110) 
# 
# U=max.size
# L=min.size
# m = nBigMatrix
# h <- (U - L)/m
# meshpts <- z1 <- z <- L + ((1:m) - 1/2) * h
# size.cen <- (18)
# 
# mean <- apply(Ksurv.mat,2, mean)
# ci <- apply(Ksurv.mat,2, HPDI, prob=.95)
# meshpts[t(ci)[,1] >0]
# meshpts[t(ci)[,2] < 0]
# 
# meshpts <- c(0,meshpts)
# 
# plot(mean ~ meshpts, ylab = TeX("Fitness effects ($\\Delta \\lambda$)"), pch="", ylim = c(-0.1,0.1), xlim = c(10,100))
# abline(h=0, lty=2)
# lines(mean ~ meshpts, lwd=1.5, col = "red")
# shade(ci, meshpts, col = col.alpha('red', 0.2))
# 
# 100-apply(Ksurv.mat, 2, LOS)
# 
# a = range(meshpts[which(apply(Ksurv.mat, 2, LOS) > 95)])
# abline(v=a[1], lty=3, col='red' )
# abline(v=a[2], lty=3, col='red' )
# 
# plot(mean ~ meshpts, ylab = TeX("Fitness effects ($\\Delta \\lambda$)"), pch="", ylim = c(-0.001,0.005), xlim = c(5,100))
# abline(h=0, lty=2)
# mean <- apply(Kgrow.mat,2, mean)
# ci <- apply(Kgrow.mat,2, HPDI, prob=.95)
# meshpts[t(ci)[,1] >0]
# meshpts[t(ci)[,2] < 0]
# 
# meshpts[apply(Kgrow.mat, 2, LOS) < 5]
# 
# a = range(meshpts[which(apply(Kgrow.mat, 2, LOS) < 5)])
# abline(v=a[1], lty=3, col='red' )
# abline(v=a[2], lty=3, col='red' )
# 
# lines(mean ~ meshpts, lwd=1.5, col = "blue")
# shade(ci, meshpts, col = col.alpha('blue', 0.2))
# 
# 
# plot(mean ~ meshpts, ylab = TeX("Fitness effects ($\\Delta \\lambda$)"), pch="", ylim = c(-0.1,0.02), xlim = c(5,100))
# abline(h=0, lty=2)
# mean <- apply(Kfec.mat,2, mean)
# ci <- apply(Kfec.mat,2, HPDI, prob=.95)
# meshpts[t(ci)[,1] >0]
# meshpts[t(ci)[,2] < 0]
# lines(mean ~ meshpts, lwd=1.5, col = "black")
# shade(ci, meshpts, col = col.alpha('black', 0.2))
# apply(Kfec.mat, 2, LOS)
# 
# 
# mean <- apply(Gfec.mat,2, mean)
# ci <- apply(Gfec.mat,2, HPDI, prob=.95)
# meshpts[t(ci)[,1] >0]
# meshpts[t(ci)[,2] < 0]
# lines(mean ~ meshpts, lwd=1.5, col = "black")
# shade(ci, meshpts, col = col.alpha('black', 0.2))
# apply(Gfec.mat, 2, LOS)
# 
