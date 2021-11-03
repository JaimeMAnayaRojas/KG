using Distributions
using CSV
using DataFrames
using LinearAlgebra
using StatsPlots
using JLD2
using RCall

cd("$(homedir())")
cd("Dropbox/Projects_JM/FSU/Pool_manipulation/KG_git/")
pwd()


# R"""
#     source("R/MainScript.R")
# """
#@rget post
post = CSV.read("Posteriors.csv", DataFrame)# Statistical analyses
# I am running the models in stan via R, and getting the posteriors for the guppy and killifish IPM# 



include("Julia/IPM_G.jl") # run from IED
## include("IPM_G.jl") # run from terminal
#Guppy_IPM(DataFrame(post[1:5,:]); nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0)




include("Julia/IPM_K.jl")
## include("IPM_K.jl")
# Killifish_IPM(DataFrame(post[1,:]); nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0)


# IPMs = [Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0),
# Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0)]



@save "IPMs.jld2" IPMs  
@load "IPMs.jld2" IPMs  



# size(G_sim)[1]

# G_sim[1][1]

#summ_tab = DataFrame()
df = IPMs[1][1]
ci = (mapcols(x -> HDI(x, credible_mass=0.95), df))
# p_val = (mapcols(x -> boot_p(x), Δ13C_net))
G = DataFrame(parameter = names(df), mean = round.(mean.(eachcol(df)), digits=3),
lc = round.(Vector(ci[1, :]), digits =3), up = round.(Vector(ci[2, :]), digits =3) );

#append!(summ_tab, summ_tab1)

df = IPMs[2][1]
ci = (mapcols(x -> HDI(x, credible_mass=0.95), df))
# p_val = (mapcols(x -> boot_p(x), Δ13C_net))
K = DataFrame(parameter = names(df), mean = round.(mean.(eachcol(df)), digits=3),
lc = round.(Vector(ci[1, :]), digits =3), up = round.(Vector(ci[2, :]), digits =3) );


sTab = [G, K]

sTab[1]
sTab[2]


## Make plots
a = 1:4
med_b = sTab[1][1:2,:mean]
med_b = append!(sTab[2][1:2,:mean], med_b)
max_b = sTab[1][1:2,:up]
max_b =append!(sTab[2][1:2,:up], max_b)
min_b = sTab[1][1:2,:lc]
min_b =append!(sTab[2][1:2,:lc], min_b)

x = [0, 4.5 ]
y = [1,1]
plot(x,y, xlims = (0.5, 4.5), c = :black, line = (:dash, 1))
x = [1, 2 , NaN, 3, 4]
y = [med_b[1], med_b[2],  NaN, med_b[3], med_b[4]]
plot!(x,y, c = :black)
scatter!((a, med_b),	yerror=(med_b-min_b, max_b-med_b), markersize = 10, 
    xlims = (0.5, 4.5), legend = false, colour = [:cyan3, :cyan3, :coral2, :coral2] )
ylabel!("Fitness (λ)")
xlabel!("Pool treatment")
xticks!([1,2,3,4], ["KG", "NK", "KG", "NG"])
savefig("Figure 1_julia.png")

#

# LTRE

## Make plots


a = 1:5
med_b = sTab[1][[9,8,5,6,7],:mean]
max_b = sTab[1][[9,8,5,6,7],:up]
min_b = sTab[1][[9,8,5,6,7],:lc]

bar((a, med_b),	yerror=(med_b-min_b, max_b-med_b), markersize = 10, 
    xlims = (0.5, 4.5), lab ="Guppy"  )

med_b = sTab[2][[9,8,5,6,7],:mean]
max_b = sTab[2][[9,8,5,6,7],:up]
min_b = sTab[2][[9,8,5,6,7],:lc]

bar((a, med_b),	yerror=(med_b-min_b, max_b-med_b), markersize = 10, xlims = (0.5, 4.5), lab ="killifish"  )


ylabel!("Fitness (λ)")
xlabel!("Pool treatment")
xticks!([1,2,3,4], ["KG", "NK", "KG", "NG"])
savefig("Figure 1_julia.png")
