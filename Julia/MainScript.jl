using Distributions
using CSV
using DataFrames
using LinearAlgebra
using StatsPlots
using JLD2
using RCall
using Plots
using Plots.Measures


# Set the working directory to my folder
cd("$(homedir())")
cd("Dropbox/Jaime M/Projects_JM/FSU/Pool_manipulation/KG_git/")
pwd()

# load my functions
include("Functions.jl")




## Run the Bayesian model via rstan.
R"""
    source("R/MainScript.R")
"""
@rget post


# Load the posteriors from the model
post = CSV.read("outputs/Posteriors.csv", DataFrame)
# I am running the models in stan via R, and getting the posteriors for the guppy and killifish IPM# 


 

# Test the IPM functions, 
# This also helps compiling the functions early and then the calculation faster. 
# In Julia the first time a function is used it takes a bit of time, after that it becomes faster and faster
# Run both IPM models

@time Guppy_IPM(post[1:3,:]; nBigMatrix = 100, min_size = 2, max_size = 45, size_cen = 18.0) # 
@time Killifish_IPM(post[1:3,:]; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0)


 nBigMatrix = 100
IPMs = [Guppy_IPM(post; nBigMatrix = nBigMatrix, min_size = 2, max_size = 45, size_cen = 18.0),
Killifish_IPM(post; nBigMatrix = nBigMatrix, min_size = 2, max_size = 110, size_cen = 18.0)]

## Save the IPMs
@save "outputs/IPMs.jld2" IPMs  

## load the IPMs
@load "outputs/IPMs.jld2" IPMs  


##### Figures vital rates

f1 = include("Man_plots.jl")
savefig(f1, "plots/Figure 1.png")

sTab[1][:,:Species] = fill("Guppy", size(sTab[1])[1])
sTab[2][:,:Species] = fill("Kilifish", size(sTab[2])[1])
sumTab = sTab[1]
append!(sumTab, sTab[2])
println(sumTab)
CSV.write("outputs/Summary_statsIPM.csv", sumTab )


include("Figure_2.jl")
ylabel!("LOS (%) \n ΔV > 0 ")
savefig("plots/Figure 2.png")


include("Figure S.jl")
savefig("plots/Figure S2.png")


plot(pSG, pSK, layout = (2,1))
xlabel!("Size (mm)")
ylabel!("Frequency (N)")
savefig("plots/Figure S1.png")
