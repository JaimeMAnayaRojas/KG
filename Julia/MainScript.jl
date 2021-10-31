using Distributions
using CSV
using DataFrames
using LinearAlgebra
# using Turing
# using DynamicHMC
using StatsPlots
using JLD2
using RCall

cd("$(homedir())")
cd("Dropbox/Projects_JM/FSU/Pool_manipulation/KG_git/")
pwd()

# Statistical analyses
# I am running the models in stan via R, and getting the posteriors for the guppy and killifish IPM

@time begin
    R"""
    source("R/MainScript.R")
    """
    include("IPM_G.jl")
    include("IPM_K.jl")

end