using Distributions
using CSV
using DataFrames
using LinearAlgebra
using StatsPlots
using JLD2
using RCall

function HDI(samples; credible_mass=0.95)
	# Computes highest density interval from a sample of representative values,
	# estimated as the shortest credible interval
	# Takes Arguments posterior_samples (samples from posterior) and credible mass (normally .95)
	# Originally from https://stackoverflow.com/questions/22284502/highest-posterior-density-region-and-central-credible-region
	# Adapted to Julialang
	sorted_points = sort(samples)
	ciIdxInc = Int(ceil(credible_mass * length(sorted_points)))
	nCIs = length(sorted_points) - ciIdxInc
	ciWidth = repeat([0.0],nCIs)
	for i in range(1, stop=nCIs)
		ciWidth[i] = sorted_points[i + ciIdxInc] - sorted_points[i]
	end
	HDImin = sorted_points[findfirst(isequal(minimum(ciWidth)),ciWidth)]
	HDImax = sorted_points[findfirst(isequal(minimum(ciWidth)),ciWidth)+ciIdxInc]
	return([HDImin, HDImax])
end


function LOS(v, b = 0)
	return 100*length(findall(v .> b)) ./length(v)
end


cd("$(homedir())")
cd("Dropbox/Projects_JM/FSU/Pool_manipulation/KG_git/")
pwd()


# R"""
#     source("R/MainScript.R")
# """
# @rget post
# I am running the models in stan via R, and getting the posteriors for the guppy and killifish IPM# 



include("IPM_G.jl") # run from IED
include("IPM_K.jl")
 

# Run both IPM models

# IPMs = [Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0),
# Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, size_cen = 18.0)]

## Save the IPMs
# @save "IPMs.jld2" IPMs  
@load "IPMs.jld2" IPMs  


##### Figures vital rates

f1 = include("Man_plots.jl")

savefig(f1, "Figure 1.png")

sTab[1][:,:Species] = fill("Guppy", size(sTab[1])[1])
sTab[2][:,:Species] = fill("Kilifish", size(sTab[2])[1])

sumTab = sTab[1]

append!(sumTab, sTab[2])

println(sumTab)
CSV.write("Summary_statsIPM.csv", sumTab )


include("Figure_2.jl")
savefig("Figure 2.png")


include("Figure S.jl")
savefig("Figure S1.png")