rm(list=ls(all=TRUE))
library("psych") 
library(rethinking)
library("plyr")
library("dplyr")
library(lme4)
library(lmerTest)
library("matlib")
library("MASS")
library("stringr")

data = read.csv("data/data for Pub.csv")

data$count <- 1

for (i in 1:dim(data)[1]){
  data[i,'stream'] <- str_split(data$Full_name[i],'-')[[1]][1]
}
unique(data$stream)

data[which(data$stream=='El_Cedro'),'drainage'] <- 'El_Cedro'

ddply(data,c('drainage','stream','Habitat'),summarise,N=sum(count))

ddply(data,c('drainage','stream','Habitat'),summarise,N=mean(d15Nb))
ddply(data,c('drainage','stream','Habitat'),summarise,N=mean(d13Cb))


center.size <- 10

names(data)
data$Length.cen <- (data$Length - center.size)

# sample sizes
ddply(data,c('SP',"drainage","stream",'Habitat'),summarise,N=sum(count))

##################################################################################################################################
##################################################################################################################################
#15N Analyses
##################################################################################################################################
##################################################################################################################################
data$diff.d15 <- (data$d15N - data$d15Nb)

data$TP <- data$diff.d15 


data$Length.cen <- (data$Length - center.size)


mod <- lmer(TP ~ 1 + FG*Habitat + FG*Habitat*Length.cen + (1|drainage), data=data )
summary(mod)
av <- anova(mod)

mod <- lmer(TP ~ 0 + FG:Habitat + FG:Habitat:Length.cen + (1|drainage), data=data )
summary(mod)


plot(mod)

rand <- ranef(mod)$drainage
rand$drainage <- dimnames(rand)[[1]]
dimnames(rand)[[2]][1] <- 'rand'

write.csv(summary(mod)$coef,"./Results/TP Parameters.csv")
write.csv(as.matrix(summary(mod)$vcov),"./Results/d15N vcov.csv")
write.csv(av,"./Results/TP Anova.csv")
#=================================================================================================================================
#15N Figures
#=================================================================================================================================
source('./R/2 d15N Figures New.R')




#=================================================================================================================================
#Contrasts
#=================================================================================================================================
coef <- summary(mod)$coef[,'Estimate']
vcov <- as.matrix(summary(mod)$vcov)
df <- summary(mod)$coef[,'df']




#Killifish natural sites slope comparison
L <- array(0,c(length(coef),2))
L[which(names(coef)=='FGKillifish:HabitatKGP:Length.cen'),] <- 1
L[which(names(coef)=='FGKillifish:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coef)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKGP:Length.cen')] <- 1/3
L.df[which(names(df)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 1/3
L.df[which(names(df)=='FGKillifish:HabitatKO:Length.cen')] <- 1/3


F.stat <- (coef %*% L %*% ginv(t(L)%*%vcov%*%L) %*% t(L) %*% coef) / 2
d <- L.df%*%df

F.stat
d

pf(F.stat, 2, d, lower.tail = FALSE, log.p = FALSE)




#Guppy natural site slope comparison
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 1
L[which(names(coef)=='FGGuppy:HabitatKGP:Length.cen')] <- -1


L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 0.5
L.df[which(names(df)=='FGGuppy:HabitatKGP:Length.cen')] <- 0.5


diff <- L%*%coef
se <- sqrt(L %*% vcov %*% L)

t <- diff / se
d <- L.df%*%df
t
d
2*pt(abs(t),df=d,lower.tail=FALSE)


 
#Killifish versus guppies in natural KG
#intercept-KG
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGKillifish:HabitatKG.nat')] <- 1
L[which(names(coef)=='FGGuppy:HabitatKG.nat')] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKG.nat')] <- 0.5
L.df[which(names(df)=='FGGuppy:HabitatKG.nat')] <- 0.5


diff <- L%*%coef
se <- sqrt(L %*% vcov %*% L)

t <- diff / se
d <- L.df%*%df
t
d
2*pt(abs(t),df=d,lower.tail=FALSE)

#Slopes
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 1
L[which(names(coef)=='FGGuppy:HabitatKG.nat:Length.cen')] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 0.5
L.df[which(names(df)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 0.5


diff <- L%*%coef
se <- sqrt(L %*% vcov %*% L)

t <- diff / se
d <- L.df%*%df
t
d
2*pt(abs(t),df=d,lower.tail=FALSE)













#Slopes All
L <- array(0,c(length(coef),3))
L[which(names(coef)=='FGGuppy:HabitatKGP:Length.cen'),] <- 1
L[which(names(coef)=='FGGuppy:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coef)=='FGGuppy:HabitatKG.exp:Length.cen'),2] <- -1
L[which(names(coef)=='FGGuppy:HabitatKG.old:Length.cen'),3] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGGuppy:HabitatKGP:Length.cen')] <- 1/4
L.df[which(names(df)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 1/4
L.df[which(names(coef)=='FGGuppy:HabitatKG.exp:Length.cen')] <- 1/4
L.df[which(names(coef)=='FGGuppy:HabitatKG.old:Length.cen')] <- 1/4


F.stat <- (coef %*% L %*% ginv(t(L)%*%vcov%*%L) %*% t(L) %*% coef) / 4
d <- L.df%*%df

F.stat
d

pf(F.stat, 3, d, lower.tail = FALSE, log.p = FALSE)




 
#Slopes All
L <- array(0,c(length(coef),4))
L[which(names(coef)=='FGKillifish:HabitatKGP:Length.cen'),] <- 1
L[which(names(coef)=='FGKillifish:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coef)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1
L[which(names(coef)=='FGKillifish:HabitatKG.exp:Length.cen'),3] <- -1
L[which(names(coef)=='FGKillifish:HabitatKG.old:Length.cen'),4] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKGP:Length.cen')] <- 1/5
L.df[which(names(df)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 1/5
L.df[which(names(df)=='FGKillifish:HabitatKO:Length.cen')] <- 1/5
L.df[which(names(coef)=='FGKillifish:HabitatKG.exp:Length.cen')] <- 1/5
L.df[which(names(coef)=='FGKillifish:HabitatKG.old:Length.cen')] <- 1/5


F.stat <- (coef %*% L %*% ginv(t(L)%*%vcov%*%L) %*% t(L) %*% coef) / 4
d <- L.df%*%df

F.stat
d

pf(F.stat, 4, d, lower.tail = FALSE, log.p = FALSE)








#=================================================================================================================================
#Similarity Contrasts for within site comparisons
#=================================================================================================================================
# Two kinds. First is the pairwise difference. Second is the squared pairwise difference or similarity.
sizes <- seq(center.size,30,length.out=25)#c(center.size:35)

#Function to calculate for N and C
dS.func.comm <- function(sizes,means,vcov,its){
  S.KG.nat <- S.KG.old <- S.KG.exp <-  S.KGP  <- array(0,c(its))

  dS.nat <- array(0,c(its))
  dS.old <- array(0,c(its))
  dS.exp <- array(0,c(its))

  D.KG.nat <- D.KG.old <- D.KG.exp <-  D.KGP  <- array(0,c(its))

  dD.nat <- array(0,c(its))
  dD.old <- array(0,c(its))
  dD.exp <- array(0,c(its))

  D.K.KGP.KG <- array(0,c(its))
  D.K.KO.KG <- array(0,c(its))

  S.K.KGP.KG <- array(0,c(its))
  S.K.KO.KG <- array(0,c(its))

  S.G.KGP.KG <- array(0,c(its))
  D.G.KGP.KG <- array(0,c(its))


  for (ii in 1:its){
    beta <- mvrnorm(1,means,vcov)

    #Nat
    b.G.KG.nat <- beta['FGGuppy:HabitatKG.nat']#post_TP$N_alpha_KG_TG[ii]
    b.z.G.KG.nat <- beta['FGGuppy:HabitatKG.nat:Length.cen']#post_TP$N_slope_KG_TG[ii]

    b.K.KG.nat <- beta['FGKillifish:HabitatKG.nat']#post_TP$N_alpha_KG_TK[ii]
    b.z.K.KG.nat <- beta['FGKillifish:HabitatKG.nat:Length.cen']#post_TP$N_slope_KG_TK[ii]

    b.G.KGP <- beta['FGGuppy:HabitatKGP']#post_TP$N_alpha_KGP_TG[ii]
    b.z.G.KGP <- beta['FGGuppy:HabitatKGP:Length.cen']#post_TP$N_slope_KGP_TG[ii]

    b.K.KGP <- beta['FGKillifish:HabitatKGP']#post_TP$N_alpha_KGP_TG[ii]
    b.z.K.KGP <- beta['FGKillifish:HabitatKGP:Length.cen']#post_TP$N_slope_KGP_TG[ii]

    b.K.KO <- beta['FGKillifish:HabitatKO']#post_TP$N_alpha_KO_TK[ii]
    b.z.K.KO <- beta['FGKillifish:HabitatKO:Length.cen']#post_TP$N_slope_KO_TK[ii]


    #old Exp
    b.G.KG.old <- beta['FGGuppy:HabitatKG.old']#post_TP$N_alpha_KG_TG[ii]
    b.z.G.KG.old <- beta['FGGuppy:HabitatKG.old:Length.cen']#post_TP$N_slope_KG_TG[ii]

    b.K.KG.old <- beta['FGKillifish:HabitatKG.old']#post_TP$N_alpha_KG_TK[ii]
    b.z.K.KG.old <- beta['FGKillifish:HabitatKG.old:Length.cen']#post_TP$N_slope_KG_TK[ii]


    #New Exp
    b.G.KG.exp <- beta['FGGuppy:HabitatKG.exp']#post_TP$N_alpha_KG_TG[ii]
    b.z.G.KG.exp <- beta['FGGuppy:HabitatKG.exp:Length.cen']#post_TP$N_slope_KG_TG[ii]

    b.K.KG.exp <- beta['FGKillifish:HabitatKG.exp']#post_TP$N_alpha_KG_TK[ii]
    b.z.K.KG.exp <- beta['FGKillifish:HabitatKG.exp:Length.cen']#post_TP$N_slope_KG_TK[ii]



    for (i in 1:length(sizes)){
      for (j in 1:length(sizes)){

        #Similarity: Square pairwise differences
        S.KG.nat[ii] <- S.KG.nat[ii] + (-( ( b.G.KG.nat + b.z.G.KG.nat*(sizes[i] - center.size) ) - ( b.K.KG.nat + b.z.K.KG.nat*(sizes[j] - center.size) ) )^2 ) / length(sizes)^2
        S.KG.old[ii] <- S.KG.old[ii] + (-( ( b.G.KG.old + b.z.G.KG.old*(sizes[i] - center.size) ) - ( b.K.KG.old + b.z.K.KG.old*(sizes[j] - center.size) ) )^2 ) / length(sizes)^2
        S.KG.exp[ii] <- S.KG.exp[ii] + (-( ( b.G.KG.exp + b.z.G.KG.exp*(sizes[i] - center.size) ) - ( b.K.KG.exp + b.z.K.KG.exp*(sizes[j] - center.size) ) )^2 ) / length(sizes)^2

        S.KGP[ii] <- S.KGP[ii] + (-( ( b.G.KGP + b.z.G.KGP*(sizes[i] - center.size) ) - ( b.K.KGP + b.z.K.KGP*(sizes[j] - center.size) ) )^2 ) / length(sizes)^2



        #Directional distance
        D.KG.nat[ii] <- D.KG.nat[ii] + (( ( b.G.KG.nat + b.z.G.KG.nat*(sizes[i] - center.size) ) - ( b.K.KG.nat + b.z.K.KG.nat*(sizes[j] - center.size) ) ) ) / length(sizes)^2
        D.KG.old[ii] <- D.KG.old[ii] + (( ( b.G.KG.old + b.z.G.KG.old*(sizes[i] - center.size) ) - ( b.K.KG.old + b.z.K.KG.old*(sizes[j] - center.size) ) ) ) / length(sizes)^2
        D.KG.exp[ii] <- D.KG.exp[ii] + (( ( b.G.KG.exp + b.z.G.KG.exp*(sizes[i] - center.size) ) - ( b.K.KG.exp + b.z.K.KG.exp*(sizes[j] - center.size) ) ) ) / length(sizes)^2

        D.KGP[ii] <- D.KGP[ii] + (( ( b.G.KGP + b.z.G.KGP*(sizes[i] - center.size) ) - ( b.K.KGP + b.z.K.KGP*(sizes[j] - center.size) ) ) ) / length(sizes)^2


        ###Killifish Differences
        S.K.KGP.KG[ii] <- S.K.KGP.KG[ii] + (-( ( b.K.KGP + b.z.K.KGP*(sizes[i] - center.size) ) - ( b.K.KG.nat + b.z.K.KG.nat*(sizes[j] - center.size) ) )^2 ) / length(sizes)^2
        S.K.KO.KG[ii] <- S.K.KO.KG[ii] + (-( ( b.K.KO + b.z.K.KO*(sizes[i] - center.size) ) - ( b.K.KG.nat + b.z.K.KG.nat*(sizes[j] - center.size) ) )^2 ) / length(sizes)^2


        D.K.KGP.KG[ii] <- D.K.KGP.KG[ii] + (( ( b.K.KGP + b.z.K.KGP*(sizes[i] - center.size) ) - ( b.K.KG.nat + b.z.K.KG.nat*(sizes[j] - center.size) ) ) ) / length(sizes)^2
        D.K.KO.KG[ii] <- D.K.KO.KG[ii] + (( ( b.K.KO + b.z.K.KO*(sizes[i] - center.size) ) - ( b.K.KG.nat + b.z.K.KG.nat*(sizes[j] - center.size) ) ) ) / length(sizes)^2

        ###Guppy Differences
        S.G.KGP.KG[ii] <- S.G.KGP.KG[ii] + (-( ( b.G.KGP + b.z.G.KGP*(sizes[i] - center.size) ) - ( b.G.KG.nat + b.z.G.KG.nat*(sizes[j] - center.size) ) )^2 ) / length(sizes)^2
        D.G.KGP.KG[ii] <- D.G.KGP.KG[ii] + (( ( b.G.KGP + b.z.G.KGP*(sizes[i] - center.size) ) - ( b.G.KG.nat + b.z.G.KG.nat*(sizes[j] - center.size) ) ) ) / length(sizes)^2


      }
    }



    dS.nat[ii] <- S.KG.nat[ii] - S.KGP[ii]
    dS.old[ii] <- S.KG.old[ii] - S.KGP[ii]
    dS.exp[ii] <- S.KG.exp[ii] - S.KGP[ii]

    dD.nat[ii] <- D.KG.nat[ii] - D.KGP[ii]
    dD.old[ii] <- D.KG.old[ii] - D.KGP[ii]
    dD.exp[ii] <- D.KG.exp[ii] - D.KGP[ii]

  }

  results <- list(dS.nat,dS.old,dS.exp,S.KG.nat,S.KG.exp,S.KGP,S.KG.old,dD.nat,dD.old,dD.exp,D.KG.nat,D.KG.exp,D.KGP,D.KG.old,D.K.KGP.KG,D.K.KO.KG,S.K.KGP.KG,S.K.KO.KG,S.G.KGP.KG,D.G.KGP.KG)
  names(results) <- c('dS.nat','dS.old','dS.exp','S.KG.nat','S.KG.exp','S.KGP','S.KG.old','dD.nat','dD.old','dD.exp','D.KG.nat','D.KG.exp','D.KGP','D.KG.old','D.K.KGP.KG','D.K.KO.KG','S.K.KGP.KG','S.K.KO.KG','S.G.KGP.KG','D.G.KGP.KG')
  return(results);

}

#Run Function
d <- dS.func.comm(sizes=sizes,means=summary(mod)$coef[,c('Estimate')],vcov=summary(mod)$vcov,its=10000)


write.csv(as.data.frame(precis(d, prob = .95)), "Results/N15_DS.csv")

precis(d, prob = .95)






##################################################################################################################################
##################################################################################################################################
#delta 13C
##################################################################################################################################
##################################################################################################################################

D <- 7.018
I <- 0.048

data$L <- 93 / ( 1 + (0.246*data$CNratio - 0.775)^(-1)) #This relies on CN data, but there is not CN data for all individuals?
data$d13Cc <- data$d13C + D * (I + (3.9/(1+287/data$L))) #Jaime had a negative after d13C

plot(data$Length,data$CNratio)
cor.test(data$Length,data$CNratio)

plot(data$Length,data$L)
cor.test(data$Length,data$L)


plot(data$d13C ~ data$CNratio)
summary(lm(data$d13C ~ data$CNratio)) # as expected there is a correlation between the CN rations, so we should normilize the data


plot(data$d13Cc ~ data$CNratio)
summary(lm(data$d13Cc ~ data$CNratio)) # as expected there is a correlation between the CN rations, so we should normilize the data


#Biplots of d13C and d15N
toplot <- data
plot(toplot[which(toplot$Habitat=='KGP' & toplot$FG=='Killifish'),'d13Cc'],toplot[which(toplot$Habitat=='KGP' & toplot$FG=='Killifish'),'TP'],
     xlim=c(-35,-22),pch=16,col='orange')
points(toplot[which(toplot$Habitat=='KGP' & toplot$FG=='Guppy'),'d13Cc'],toplot[which(toplot$Habitat=='KGP' & toplot$FG=='Guppy'),'TP'],pch=16,col='grey')


plot(toplot[which(toplot$Habitat=='KG.nat' & toplot$FG=='Killifish'),'d13Cc'],toplot[which(toplot$Habitat=='KG.nat' & toplot$FG=='Killifish'),'TP'],
     xlim=c(-35,-22),pch=16,col='orange')
points(toplot[which(toplot$Habitat=='KG.nat' & toplot$FG=='Guppy'),'d13Cc'],toplot[which(toplot$Habitat=='KG.nat' & toplot$FG=='Guppy'),'TP'],pch=16,col='grey')




data$diff.d13Cc <- data$d13Cc - data$d13Cb


#=================================================================================================================================
#Stats
#=================================================================================================================================

mod <- lmer(diff.d13Cc ~ 1 + FG*Habitat + FG*Habitat*Length.cen + (1|drainage), data=data )
summary(mod)
av <- anova(mod)

mod <- lmer(diff.d13Cc ~ 0 + FG:Habitat + FG:Habitat:Length.cen + (1|drainage), data=data )
summary(mod)


plot(mod)



rand <- ranef(mod)$drainage
rand$drainage <- dimnames(rand)[[1]]
dimnames(rand)[[2]][1] <- 'randC'


#=================================================================================================================================
#plots
#=================================================================================================================================
source('./R/2 d13C Figures New.R')


#=================================================================================================================================
#Contrasts
#=================================================================================================================================
coef <- summary(mod)$coef[,'Estimate']
vcov <- as.matrix(summary(mod)$vcov)
df <- summary(mod)$coef[,'df']




#Killifish natural sites slope comparison
L <- array(0,c(length(coef),2))
L[which(names(coef)=='FGKillifish:HabitatKGP:Length.cen'),] <- 1
L[which(names(coef)=='FGKillifish:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coef)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKGP:Length.cen')] <- 1/3
L.df[which(names(df)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 1/3
L.df[which(names(df)=='FGKillifish:HabitatKO:Length.cen')] <- 1/3


F.stat <- (coef %*% L %*% ginv(t(L)%*%vcov%*%L) %*% t(L) %*% coef) / 2
d <- L.df%*%df

F.stat
d

pf(F.stat, 2, d, lower.tail = FALSE, log.p = FALSE)




#Guppy natural site slope comparison
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 1
L[which(names(coef)=='FGGuppy:HabitatKGP:Length.cen')] <- -1


L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 0.5
L.df[which(names(df)=='FGGuppy:HabitatKGP:Length.cen')] <- 0.5


diff <- L%*%coef
se <- sqrt(L %*% vcov %*% L)

t <- diff / se
d <- L.df%*%df
t
d
2*pt(abs(t),df=d,lower.tail=FALSE)



#Killifish versus guppies in natural KG
#intercept-KG
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGKillifish:HabitatKG.nat')] <- 1
L[which(names(coef)=='FGGuppy:HabitatKG.nat')] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKG.nat')] <- 0.5
L.df[which(names(df)=='FGGuppy:HabitatKG.nat')] <- 0.5


diff <- L%*%coef
se <- sqrt(L %*% vcov %*% L)

t <- diff / se
d <- L.df%*%df
t
d
2*pt(abs(t),df=d,lower.tail=FALSE)

#Slopes
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 1
L[which(names(coef)=='FGGuppy:HabitatKG.nat:Length.cen')] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 0.5
L.df[which(names(df)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 0.5


diff <- L%*%coef
se <- sqrt(L %*% vcov %*% L)

t <- diff / se
d <- L.df%*%df
t
d
2*pt(abs(t),df=d,lower.tail=FALSE)













#Slopes All
L <- array(0,c(length(coef),3))
L[which(names(coef)=='FGGuppy:HabitatKGP:Length.cen'),] <- 1
L[which(names(coef)=='FGGuppy:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coef)=='FGGuppy:HabitatKG.exp:Length.cen'),2] <- -1
L[which(names(coef)=='FGGuppy:HabitatKG.old:Length.cen'),3] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGGuppy:HabitatKGP:Length.cen')] <- 1/4
L.df[which(names(df)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 1/4
L.df[which(names(coef)=='FGGuppy:HabitatKG.exp:Length.cen')] <- 1/4
L.df[which(names(coef)=='FGGuppy:HabitatKG.old:Length.cen')] <- 1/4


F.stat <- (coef %*% L %*% ginv(t(L)%*%vcov%*%L) %*% t(L) %*% coef) / 4
d <- L.df%*%df

F.stat
d

pf(F.stat, 3, d, lower.tail = FALSE, log.p = FALSE)





#Slopes All
L <- array(0,c(length(coef),4))
L[which(names(coef)=='FGKillifish:HabitatKGP:Length.cen'),] <- 1
L[which(names(coef)=='FGKillifish:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coef)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1
L[which(names(coef)=='FGKillifish:HabitatKG.exp:Length.cen'),3] <- -1
L[which(names(coef)=='FGKillifish:HabitatKG.old:Length.cen'),4] <- -1

L.df <- array(0,c(length(coef)))
L.df[which(names(df)=='FGKillifish:HabitatKGP:Length.cen')] <- 1/5
L.df[which(names(df)=='FGKillifish:HabitatKG.nat:Length.cen')] <- 1/5
L.df[which(names(df)=='FGKillifish:HabitatKO:Length.cen')] <- 1/5
L.df[which(names(coef)=='FGKillifish:HabitatKG.exp:Length.cen')] <- 1/5
L.df[which(names(coef)=='FGKillifish:HabitatKG.old:Length.cen')] <- 1/5


F.stat <- (coef %*% L %*% ginv(t(L)%*%vcov%*%L) %*% t(L) %*% coef) / 4
d <- L.df%*%df

F.stat
d

pf(F.stat, 4, d, lower.tail = FALSE, log.p = FALSE)






write.csv(summary(mod)$coef,"./Results/d13C Parameters.csv")
write.csv(as.matrix(summary(mod)$vcov),"./Results/d13C vcov.csv")
write.csv(av,"./Results/d13C Anova.csv")

#=================================================================================================================================
#Similarity and Differences
#=================================================================================================================================
# Two kinds. First is the pairwise differences. Second is the squared pairwise difference or similarity.
sizes <- seq(center.size,30,length.out=25)#c(center.size:35)

d <- dS.func.comm(sizes=sizes,means=summary(mod)$coef[,c('Estimate')],vcov=summary(mod)$vcov,its=10000)


write.csv(as.data.frame(precis(d, prob = .95)), "Results/C13_DS.csv")

precis(d, prob = .95)


#################################################################################################################
#Figure 3
################################################################################################################
#Center on 20mm
data$Length.cen20 <- data$Length - 20
nitro.mod <- lmer(TP ~ 0 + FG:Habitat + FG:Habitat:Length.cen20 + (1|drainage), data=data)
summary(nitro.mod)


carbo.mod <- lmer(diff.d13Cc ~ 0 + FG:Habitat + FG:Habitat:Length.cen20 + (1|drainage), data=data )
summary(carbo.mod)


nitro.coef <- summary(nitro.mod)$coef[,'Estimate']
nitro.vcov <- as.matrix(summary(nitro.mod)$vcov)
nitro.df <- summary(nitro.mod)$coef[,'df']


carbo.coef <- summary(carbo.mod)$coef[,'Estimate']
carbo.vcov <- as.matrix(summary(carbo.mod)$vcov)
carbo.df <- summary(carbo.mod)$coef[,'df']



L.N <- array(0,c(length(nitro.coef),5))

L.N[which(names(nitro.coef)=='FGGuppy:HabitatKGP'),1] <- 1
L.N[which(names(nitro.coef)=='FGGuppy:HabitatKG.nat'),2] <- 1
L.N[which(names(nitro.coef)=='FGKillifish:HabitatKG.nat'),3] <- 1
L.N[which(names(nitro.coef)=='FGKillifish:HabitatKO'),4] <- 1
L.N[which(names(nitro.coef)=='FGKillifish:HabitatKGP'),5] <- 1


nitro.means <- nitro.coef%*%L.N
nitro.errors <- sqrt(diag(t(L.N)%*%nitro.vcov%*%L.N))


L.C <- array(0,c(length(carbo.coef),5))
L.C[which(names(carbo.coef)=='FGGuppy:HabitatKGP'),1] <- 1
L.C[which(names(carbo.coef)=='FGGuppy:HabitatKG.nat'),2] <- 1
L.C[which(names(carbo.coef)=='FGKillifish:HabitatKG.nat'),3] <- 1
L.C[which(names(carbo.coef)=='FGKillifish:HabitatKO'),4] <- 1
L.C[which(names(carbo.coef)=='FGKillifish:HabitatKGP'),5] <- 1


carbo.means <- carbo.coef%*%L.C
carbo.errors <- sqrt(diag(t(L.C)%*%carbo.vcov%*%L.C))


means <- cbind(t(carbo.means),t(nitro.means))
errors <- cbind((carbo.errors),(nitro.errors))


jpeg(file = "./Figures/Figure-4.jpeg",width = 5.0, height = 5.0, units='in', res=600)

par(mfrow=c(1,1), mar = c(4, 5 ,2, 1), oma = c(0.5, 1, 1, 0.5))
scale.min <- 0.75
scale.max <- 1.25
plot(means[,1],means[,2],xlim=c(0.1,6),ylim=c(3,4.5),
     pch=16,cex=2,bty="L",col=c('grey','grey','orange','orange','orange'),
     xlab=expression(paste(Delta^{13}, "C (\u2030)")),ylab=expression(paste(Delta^{15}, "N (\u2030)")))


# Vertical arrow
arrows(y0=means[,2]-errors[,2], x0=means[,1], y1=means[,2]+errors[,2], x1=means[,1], code=3, angle=90, length=0.05, col=c('grey','grey','orange','orange','orange'), lwd=2)
# Horizontal arrow
arrows(x0=means[,1]-errors[,1], y0=means[,2], x1=means[,1]+errors[,1], y1=means[,2], code=3, angle=90, length=0.05, col=c('grey','grey','orange','orange','orange'), lwd=2)
points(means[,1],means[,2],pch=1,cex=2)

legend('topleft',c('Guppy','Killifish'),pch=c(16,16),col=c('grey','orange'),bty='n',cex=1.25)


text(x=0.87, y = 3.4, labels = 'KG-Nat')
text(x=4.7, y = 3.95, labels = 'KGP')

text(x=2.4, y = 3.85, labels = 'KO')
text(x=2.4, y = 3.35, labels = 'KG-Nat')
text(x=4.2, y = 3.3, labels = 'KGP')


graphics.off()






