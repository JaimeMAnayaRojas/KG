
library(data.table)
library(reshape)

rm(list=ls(all=TRUE))
###########################################################################################################
# First get the data

setwd("~/Dropbox/Projects_JM/FSU/Pool_manipulation/Pool_R/")


### Check the data for mistakes
g.data <- fread("data/Guppy.csv")
g.data$site = paste(g.data$Location, g.data$Pool_1, sep='-') 

k.data <- fread("data/Killifish.csv")
k.data$site = paste(k.data$Location, k.data$Pool_1, sep='-') 


env.data <- fread("data/EnvData.csv")
env.data$site = paste(env.data$Population, env.data$Pool_1, sep='-')

Rep.data <- fread("data/Rep.data.csv")
G_Rep <- Rep.data[sp== 1,,]
K_Rep <- Rep.data[sp== 0,,]

G_R = merge(g.data, G_Rep[, sp:= NULL], by= "ID", all = T)
K_R = merge(k.data, K_Rep[, sp:= NULL], by= "ID", all = T)

G_RE = merge(G_R, env.data, by = "site", all.x = T)
K_RE = merge(K_R, env.data, by = "site", all.x = T)


write.csv(G_RE,"data/G_data_to_check.csv", sep = ',')
write.csv(K_RE,"data/K_data_to_check.csv", sep = ',')

