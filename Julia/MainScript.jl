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

post = CSV.read("Posteriors.csv", DataFrame)# Statistical analyses
# I am running the models in stan via R, and getting the posteriors for the guppy and killifish IPM# R"""
# source("R/MainScript.R")
# """


# #include("Julia/IPM_G.jl") # run from IED
# include("IPM_G.jl") # run from terminal
# Guppy_IPM(DataFrame(post[1:5,:]); nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0, area = 0.0, canopy =0.0)
   

 a = Vector(-1.0:0.5:1.0)
# G_sim =[Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0, area = 0.0 , canopy = a[1]),
# Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0, area = 0.0 , canopy = a[2]),
# Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0, area = 0.0 , canopy = a[3]),
# Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0, area = 0.0 , canopy = a[4]),
# Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0, area = 0.0 , canopy = a[5])]
# @save "Guppy_IPMs.jld2" G_sim  
# #@load "Guppy_IPMs.jld2" G_sim  


#include("Julia/IPM_K.jl")
include("IPM_K.jl")
Killifish_IPM(DataFrame(post[1,:]); nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0, area = 0.0 , canopy = a[3])


K_sim = [Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0, area = 0.0 , canopy = a[1]),
Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0, area = 0.0 , canopy = a[2]),
Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0, area = 0.0 , canopy = a[3]),
Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0, area = 0.0 , canopy = a[4]),
Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0, area = 0.0 , canopy = a[5])]
@save "Killifish_IPMs.jld2" K_sim  
# #@load "Guppy_IPMs.jld2" G_sim  



# size(G_sim)[1]

# G_sim[1][1]

# summ_tab = DataFrame()
# df = G_sim[1][1]
# ci = (mapcols(x -> HDI(x, credible_mass=0.95), df))
# # p_val = (mapcols(x -> boot_p(x), Δ13C_net))
# summ_tab1 = DataFrame(parameter = names(df), mean = round.(mean.(eachcol(df)), digits=3),
# lc = round.(Vector(ci[1, :]), digits =3), up = round.(Vector(ci[2, :]), digits =3) );
# summ_tab1.Canopy = zeros(size(summ_tab1)[1]) .+ a[1]
# append!(summ_tab, summ_tab1)
                                              
# df = G_sim[2][1]
# ci = (mapcols(x -> HDI(x, credible_mass=0.95), df))
# # p_val = (mapcols(x -> boot_p(x), Δ13C_net))
# summ_tab2 = DataFrame(parameter = names(df), mean = round.(mean.(eachcol(df)), digits=3),
# lc = round.(Vector(ci[1, :]), digits =3), up = round.(Vector(ci[2, :]), digits =3) );
# summ_tab2.Canopy = zeros(size(summ_tab2)[1]) .+ a[2]
# append!(summ_tab, summ_tab2)

# df = G_sim[3][1]
# ci = (mapcols(x -> HDI(x, credible_mass=0.95), df))
# # p_val = (mapcols(x -> boot_p(x), Δ13C_net))
# summ_tab3 = DataFrame(parameter = names(df), mean = round.(mean.(eachcol(df)), digits=3),
# lc = round.(Vector(ci[1, :]), digits =3), up = round.(Vector(ci[2, :]), digits =3) );
# summ_tab3.Canopy = zeros(size(summ_tab3)[1]) .+ a[3]
# append!(summ_tab, summ_tab3)


# df = G_sim[4][1]
# ci = (mapcols(x -> HDI(x, credible_mass=0.95), df))
# # p_val = (mapcols(x -> boot_p(x), Δ13C_net))
# summ_tab4 = DataFrame(parameter = names(df), mean = round.(mean.(eachcol(df)), digits=3),
# lc = round.(Vector(ci[1, :]), digits =3), up = round.(Vector(ci[2, :]), digits =3) );
# summ_tab4.Canopy = zeros(size(summ_tab4)[1]) .+ a[4]
# append!(summ_tab, summ_tab4)


# df = G_sim[5][1]
# ci = (mapcols(x -> HDI(x, credible_mass=0.95), df))
# # p_val = (mapcols(x -> boot_p(x), Δ13C_net))
# summ_tab5 = DataFrame(parameter = names(df), mean = round.(mean.(eachcol(df)), digits=3),
# lc = round.(Vector(ci[1, :]), digits =3), up = round.(Vector(ci[2, :]), digits =3) );
# summ_tab5.Canopy = zeros(size(summ_tab5)[1]) .+ a[5]
# append!(summ_tab, summ_tab5)

# # CSV.write("G_IPM_sum.csv", summ_tab )
# df = filter(:parameter => x -> x == "lam_GR" || x == "lam_NR" || x == "delta_lam"  , summ_tab)

# println(df)
