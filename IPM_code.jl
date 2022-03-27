### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 283e2bac-aabf-11ec-1e3e-ed84e95588ca
md"""
# Supplementary Material for:
## Do interspecific interactions still play a structuring role in natural fish communities after coexistence?
"""

# ╔═╡ 89453239-97e9-4d6a-89a0-905118574bf4
md"""
## Estimation of the demographic rates 

We estimated the demographic rates used for the vital rate functions of the IPMs using bayesian methods via stan as follow:
"""


# ╔═╡ b788ec99-b506-4f8d-b725-e1278b457d92


# ╔═╡ 02227c46-5a15-484d-a302-0f2a4a51c749
md"""
## Model construction and analysis

### Integral Projection model for guppies
The model for guppies is a single-sex (female only) integral projection model (IPM) written:

$n(z',t+1) = \int [G(z',z)S(z) + D(z',z)F(z)S(z)]n(z,t)dz$

The model is a discrete-time population projection model with the population structured by body size ($$z$$). The function $$n(z',t+1)$$ describes the density of individuals in the population at the beginning of the experiment, such that $$\int_{4}^{35} n(z,t)dz$$ is the number of individuals between standard length 4 mm and 35mm. No guppies ever reach these extreme sizes and these extreme values are used to ensure no boundary effects in the models.  $$S(z)$$ is a continuous function describing the probability of an individual with body length $$z$$ at the beginning of the experiment surviving to the end experiment.  $$F(z)$$ is a continuous function describing the mean number of offspring produced by an individual with body length $$z$$ at reproduction. The model assumes that all reproducing individuals do so at the end of the interval, immediately preceding the next census.  $$G(z',z)$$   is the conditional probability density function describing transitions from length $$z$$ at the beginning of the experiment to length $$z'$$ at the end of the experiment. D(z',z) is the conditional probability density function describing the distribution of offspring body length $$z'$$ at the end of the experiment produced by parents with body length $$z$$ at the beginning of the experiment. The demographic rates used for these functions were estimated using regression methods (see main text).
"""

# ╔═╡ 20235581-b113-4445-a8b2-c9407f145a1c
# Vital rate functions
	## Surival function
	function s_z(df::AbstractDataFrame, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		α= df.α_surv[row]
		β= df.βz_surv[row]
		linear_p = α .+ β * (z .- size_cen)  # linear predictor
		p = 1 ./(1 .+ exp.(-linear_p))
		p = diagm(p)
		return(p)
	end

# ╔═╡ b60e91fc-fb4c-417f-8713-9e21a33d6825
function g_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		α= df.α_grow[row]
		β= df.βz_grow[row]		
		σ= df.σ_grow[row]
		p_den_grow = zeros(size(z)[1],size(z)[1])
		μ = α .+ β * (z .- size_cen ) 
		for i in 1:nBigMatrix
			p_den_grow[:,i] = pdf.(Normal(μ[i], σ), z1).*h
		end
		return(p_den_grow)
	end

# ╔═╡ 2cb8caf7-d2be-4cd6-aadc-67e7ccc7d58a
## Recruitment function (N.B - from birth in spring to first summer), logistic regression
	function pr_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
		α= df.α_fec[row]
		β= df.βz_fec[row]
		linear_p = α .+ β * (z .- size_cen)     # linear predictor
		p = exp.(linear_p)*(1/2)
		p = diagm(p)
		return(p)
	end

# ╔═╡ 60303de8-71ff-455c-9a3c-e2f450fdab7c
## Recruit size function
	function c_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		#α= df.rcz_int[row]
		#β= df.rcz_z[row]
		σ= 0.4 
		p_den_grow = zeros(nBigMatrix,nBigMatrix)
		μ = 7 .+ 0*(z .- size_cen) 
		for i in 1:nBigMatrix
			p_df = pdf.(Normal(μ[i], σ), z1)*h
			for j in 1:nBigMatrix
			p_den_grow[j,i] = p_df[j]
			end
		end
		return(p_den_grow)
	end

# ╔═╡ fca2ea39-40b7-42a7-9cfe-b1b9aeb8e2f0
md"""
### Integral Projection model for killifish
The model for killifish is also a single-sex (female only) IPM, but is slightly different than the one for guppies. Because killifish lay eggs that hatch at two weeks, two broods are produced per time-step. This means that the overall projection kernel for a month is a composite of two, two-week projection kernels. Each two-week projection equation is: 
"""

# ╔═╡ 648dcd5b-c356-409e-aab8-a19475734ae7


# ╔═╡ Cell order:
# ╟─283e2bac-aabf-11ec-1e3e-ed84e95588ca
# ╟─89453239-97e9-4d6a-89a0-905118574bf4
# ╠═b788ec99-b506-4f8d-b725-e1278b457d92
# ╟─02227c46-5a15-484d-a302-0f2a4a51c749
# ╠═20235581-b113-4445-a8b2-c9407f145a1c
# ╠═b60e91fc-fb4c-417f-8713-9e21a33d6825
# ╠═2cb8caf7-d2be-4cd6-aadc-67e7ccc7d58a
# ╠═60303de8-71ff-455c-9a3c-e2f450fdab7c
# ╟─fca2ea39-40b7-42a7-9cfe-b1b9aeb8e2f0
# ╠═648dcd5b-c356-409e-aab8-a19475734ae7
