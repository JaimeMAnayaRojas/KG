### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ 5736dcb2-0be4-4be8-9fcb-21d679caf199
# Load the necessary modules
begin
	using Distributions
	using CSV
	using DataFrames
	using LinearAlgebra
	using StatsPlots
	using JLD2
	using RCall
end

# ╔═╡ 283e2bac-aabf-11ec-1e3e-ed84e95588ca
md"""
# Details on the data analyses and implementation of the IMP models:
In this notebook, we provide a detailed description of our modeling framework and the running of the code:

#### Part I: Overview of the Integral Projection Models (IPM)
	1. Running and construction of the Bayesian model and extracting the posterior
	2. Constructing the IPMs
	3. Vital functions and validation
	4. Fitness calculations
	5. Life Table Response Experiment (LTRE) Analyses

#### Part II: Analyses for the manuscript
	1. Running the full IPMS and LTREs
	2. Plots and interpretation
"""

# ╔═╡ b6c40268-0a24-4691-aff9-a487eb9846e8
md"""
## Part I: Overview of the Integral Projection Models (IPM)
"""

# ╔═╡ 89453239-97e9-4d6a-89a0-905118574bf4
md"""
### 1. Running and construction of the Bayesian model and extracting the posterior
##### Estimation of the demographic rates 

We estimated the demographic rates used for the vital rate functions of the IPMs using bayesian methods via Stan as follow:

##### Survival
$S_{i} \approx Binomial(1, p_{i})$
$p_{i} = \alpha_{S(z)} + \beta_{S(z),NK}NK_{i} + \beta_{S(z),z}(z_{i} - \bar{z}) +$  
$\beta_{S(z), zNK}NK_{i} (z_{i} - \bar{z}) +  \beta_{S(z), area}Area_{i} + \beta_{S(z), canopy} + Canopy_{i} + v_{S(z),[stream,i]}$

Priors:

$\sigma_{stream} \approx Cauchy(0, 1)$
$v_{S(z)} \approx Normal(0, \sigma_{stream})$
$\beta_{S(z),NK} \approx Normal(0, 1)$
$\beta_{S(z),z} \approx Normal(0, 1)$
$\beta_{S(z),zNK} \approx Normal(0, 1)$
$\beta_{S(z),area} \approx Normal(0, 1)$
$\beta_{S(z),canopy} \approx Normal(0, 1)$
$\alpha_{S(z)} \approx Normal( 2.6 , 3 )$ 

##### Somathic growth
$z'_{i} \approx Normal(\mu_{i}, \sigma)$
$\mu_{i} = \alpha_{G(z',z)} + \beta_{G(z',z),NK}NK_{i} + \beta_{G(z',z),z}(z_{i} - \bar{z}) +$
$\beta_{G(z',z), zNK}NK_{i} (z_{i} - \bar{z}) + \beta_{G(z',z), area}Area_{i} + \beta_{G(z',z), canopy} + Canopy_{i} + v_{G(z',z),[stream,i]}$

Priors:

$\sigma_{G(z',z)} \approx Cauchy(0, 2)$
$v_{G(z',z)} \approx Normal(0, \sigma_{stream})$
$\beta_{G(z',z),NK} \approx Normal(0, 1)$
$\beta_{G(z',z),z} \approx Normal(0, 1)$
$\beta_{G(z',z),zNK} \approx Normal(0, 1)$
$\beta_{G(z',z),area} \approx Normal(0, 1)$
$\beta_{G(z',z),canopy} \approx Normal(0, 1)$
$\alpha_{G(z',z)} \approx Normal(18, 10)$ 



##### Fecundity
$R_{i} \approx Poisson(\lambda_{i})$
$\lambda_{i} = \alpha_{G(z',z)} + \beta_{G(z',z),NK}NK_{i} + \beta_{G(z',z),z}(z_{i} - \bar{z}) +$
$\beta_{R(z), zNK}NK_{i} (z_{i} - \bar{z}) + \beta_{R(z), area}Area_{i} + \beta_{R(z), canopy} + Canopy_{i} + v_{R(z),[stream,i]}$

Priors:

$\sigma_{R(z)} \approx Cauchy(0, 2)$
$v_{R(z)} \approx Normal(0, \sigma_{stream})$
$\beta_{R(z),NK} \approx Normal(0, 1)$
$\beta_{R(z),z} \approx Normal(0, 1)$
$\beta_{R(z),zNK} \approx Normal(0, 1)$
$\beta_{R(z),area} \approx Normal(0, 1)$
$\beta_{R(z),canopy} \approx Normal(0, 1)$
$\alpha_{R(z)} \approx Normal(18, 10)$ 

"""

# ╔═╡ dc2034eb-cc01-4032-824b-9e0ace72e329
md"""
To combine demographic rates and to control for any covariation among demographic rates, we use the Bayesian Hierarchical Centering (BHC) method following Elderd and Miller (2016) by simultaneously running all vital rate models for both species and estimating a global random stream drainage effect parameter for each species at the intercept ($$\sigma_{stream}$$).  We use the same model structure for both species, but only the priors for survival were different between guppies and killifish.
"""

# ╔═╡ 560f7b39-731b-4a8e-a16b-b822d551073f
md"""

### Model executions
In the next cell, we run the Bayesian model via R, using the rstan and rethinking packages as interface.
"""

# ╔═╡ 162e9c0c-e0f3-4e1d-9a4e-a9d3694d361b
md"""
##### If you wish to run the model on your computer run the cell below:
"""

# ╔═╡ b12084bb-8040-40d1-9657-2800e5b6ccb2
# R"""
#source("R/MainScript.R") # Runs the model via R
# """

# ╔═╡ 6bb7c116-8c16-4df6-8f34-501f652faa12
md"""
#### Alternatively, you can load an already run model and start using the posteriors. Next lines...
"""

# ╔═╡ ce97c58b-c9c5-4aa3-962a-7e19e1fb26c1
R"""
library("rethinking")
mod <- readRDS("Model_KG.RDS")
sum = as.data.frame(precis(mod, digits = 5, prob = .95, depth = 2))
names = rownames(sum)
post <- as.data.frame(extract(mod))
""";

# ╔═╡ ad66cf19-db0f-4f54-9c48-deb72543e7ad
# get the summary table, the names of the rows, and the posterios
begin
	@rget sum names post;
	sumTab = DataFrame(sum)
	sumTab.Pars = names
	sumTab[:, [7,1,2,3,4,5,6] ]
end


# ╔═╡ 02227c46-5a15-484d-a302-0f2a4a51c749
md"""
## 2. Constructing the IPMs

### Integral Projection model for guppies
The model for guppies is a single-sex (female only) integral projection model (IPM) written:

$n(z',t+1) = \int_L^U K(z',z)n(z,t)dz,$

The model is a discrete-time population projection model with the population structured by body size ($$z$$). The function $$n(z',t+1)$$ describes the density of individuals in the population at the beginning of the experiment, such that $$\int_{4}^{45} n(z,t)dz$$ is the number of individuals between standard length 4 mm and 45mm. $K(z',z)$ is the complete IPM kernel, which is composed of the survival and fecundity kernels:

$n(z',t+1) = \int_L^U [P(z',z)+F(z'z)] n(z,t)dz$


The kernel, $P(z',z)h$, is the probability that an individual of size $z$ survives between $t$ and $t+1$ and grows to size z'. $F(z',z)h$ is the number of offspring of size $z'$ produced by an individual of size $z$ survives between $t$ and $t+1$. Where $h$ is the interval between size classes and if $h \to 0$ the approximation becomes exact. Hence, the survival kernel, $P(z',z)$ describes the survival and possible growth or shrinkage of individuals of size $z$ within the census interval (here, 25 days). This kernel is composed of two vital functions: $P(z',z)$ =  $$S(z) G(z',z)$$. Where   $$S(z)$$ is a continuous function describing the probability of an individual with body length $$z$$ at the beginning of the experiment surviving to the end experiment. $$G(z',z)$$, on the other hand,  is the conditional probability density function describing transitions from length $$z$$ at the beginning of the experiment to length $$z'$$ at the end of the experiment.

The fecundity or reproduction kernel, $F(z',z) = C(z',z)R(z)S(z)$ is composed of the survival function $S(z)$ and two extra vital functions. $C(z',z)$ is the conditional probability density function describing the distribution of offspring of body length $$z'$$ at the end of the experiment produced by parents with body length $$z$$ at the beginning of the experiment. $$R(z)$$ is a continuous function describing the mean number of offspring produced by an individual with body length $$z$$ at reproduction. The model assumes that all reproducing individuals do so at the end of the interval, immediately preceding the next census. The IMP model can therefore be written as:

$n(z',t+1) = \int_L^U [G(z',z)S(z) + C(z',z)R(z)S(z)]n(z,t)dz$

"""

# ╔═╡ ed107163-cd0d-4f5b-8840-3859f2363cf9
md"""
### 3. Vital functions and validations
"""

# ╔═╡ 95769a16-6e17-4d10-92d4-ec38f82f7aad
md"""
#### Survival function, S(z):
$S(z) = \frac{1}{1 + 10^{-[\alpha + \beta_{z}(z-\bar{z})]}}$
"""

# ╔═╡ 1c358e1c-a072-4273-bf7f-9d24ab84ea46
function s_z(df::AbstractDataFrame, z::AbstractVector, 
	size_cen::AbstractFloat, row::Integer)
	α= df.α_surv[row]
	β= df.βz_surv[row]
	linear_p = α .+ β * (z .- size_cen)  # linear predictor
	p = 1 ./(1 .+ exp.(-linear_p))
	p = diagm(p)
	return(p)
end

# ╔═╡ a4bcbf5d-2cdf-444c-9c1f-2acd8151250a
md"""
#### Growth function, G(z',z):
Probability of growint to $$z'$$ given $$z$$

$G(z',z) = \frac{1}{\sqrt[2]2\pi\sigma^2}e^{\frac{-[z'- (\alpha + \beta_{z}(z-\bar{z}))]^2}{2\sigma^2}}$
"""

# ╔═╡ c60cd474-89af-439c-b457-b302ccb1e2d8
md"""
### Recruitment function, R(z):
Number of offspring for an individual of size $$z$$, divided by two to account for females

$$R(z) = e^{\alpha + \beta_{z} (z-\bar{z})}  \frac{1}{2}$$ 
"""

# ╔═╡ d3fedac0-0a07-4f2d-bda9-ecc5643e23b9
function R_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
	α= df.α_fec[row]
	β= df.βz_fec[row]
	linear_p = α .+ β * (z .- size_cen)     # linear predictor
	p = exp.(linear_p)*(1/2)
	p = diagm(p)
	return(p)
end

# ╔═╡ bf727831-4cef-4e7e-adda-6c47ab7c6ac2
md"""
#### Offspring size, C(z',z):
Probability of having offspring with mean size $$z'$$ given the size of the mother $$z$$

$C(z',z) = \frac{1}{\sqrt[2]2\pi\sigma^2}e^{\frac{-[z'- (\alpha + \beta_{z}(z-\bar{z}))]^2}{2\sigma^2}}$

Because we have no good estimation of the parameters of this function, we assumed that size has no effect on the size of offspring and when born they are approximately $7.0 mm$ $(\pm 0.4)$
"""

# ╔═╡ fed2e2fd-bce9-46ff-8f48-360a86136209
md"""
### Kernel calculations
"""

# ╔═╡ 7b44a540-00b2-4db9-91e4-2e782e5b77d4
md"""
### 4. Fitness calculations
#### IPM construction
The continuous size IPM’s are converted to matrix models prior to analysis using the midpoint rule of calculus (Easterling et al. 2000). The guppy model is a standard pre-breeding census model with the matrix model approximation: 

$n(t+1) = An(t),$

where $A$ and $n(t)$ are structured by body size ($z$). The minimum (L) and maximum (U) size ranges for guppies were 2mm and 45mm with 100 mesh points. The projection matrix, A, is made up of matrices that describe the individual demographic rates:

$A = P + CRS$

"""

# ╔═╡ 13db1664-0f55-4e32-b59c-6f25309e93b6
begin
	U= 45
	L=2
	nBigMatrix = 100
	m= nBigMatrix
	h = (U - L)/m
	z1 = zeros(nBigMatrix)
	z1 =  L .+ (collect(1:m) .- 0.5) * h
	z = z1
end;

# ╔═╡ 87d2c9f8-b26c-4c5e-9f59-33c7f73971fc
function g_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
	size_cen::AbstractFloat, row::Integer)
	α= df.α_grow[row]
	β= df.βz_grow[row]
	σ= df.σ_grow[row]
	p = zeros(size(z)[1],size(z)[1])
	μ = α .+ β * (z .- size_cen ) 
	for i in 1:nBigMatrix
		p[:,i] = pdf.(Normal(μ[i], σ), z1).*h
	end
	return(p)
end

# ╔═╡ f504b0a1-fb69-4069-a8ca-76e815438c5d
## Survival kernel
function P_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
	size_cen::AbstractFloat, row::Integer)
	out = g_z1z(df, z1, z, size_cen, row) * s_z(df, z, size_cen, row)
	return(out)

end

# ╔═╡ 9d241282-e6f2-4d63-8cf2-6d49dde03ac0
function c_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
	#α= df.rcz_int[row]
	#β= df.rcz_z[row]
	σ= 0.4 #pars.rcz_sd[row]
	#   βa= df.β_area_rcz[row]
	#   βc= df.β_canopy_rcz[row]
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

# ╔═╡ e242e6dc-f7f2-4663-a893-a6bfbcdafd8f
# Reproduction kernel
function F_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector,
	size_cen::AbstractFloat, row::Integer)
	out = c_z1z(df, z1, z, size_cen, row) * R_z(df, z, size_cen, row) * s_z(df, z, size_cen, row)
	return(out)
end


# ╔═╡ b58511e7-0703-4e60-8119-7c9a8e01fbd1
# IPM kernel. Put together the kernels
function mk_K(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector,
	size_cen::AbstractFloat, row::Integer)
	F = F_z1z(df, z1, z, size_cen, row)
	P = P_z1z(df, z1, z, size_cen, row)
	K = P + F
	out = Dict("K"=> K, "meshpts" => z, "P" => P, "F"=> F)
	return(out)
end

# ╔═╡ fe287182-f606-43f3-8c30-def2d342af54
md"""
To convert the functions into matrices, we calculate the matrix for each vital rate for both NK and KG guppies using the posterior parameters. Therefore, we got the posterior for each treatment as separate dataframes: 
"""

# ╔═╡ 5467e312-3522-4343-9fcb-fe1c3fe925dd
md"""
##### Getting the parameter values
We have 6000 posterior samples from the Bayesian model to run the IPMs 6000 times. To run the guppy IPM, for now, we just use one of those 6000 posteriors:
"""

# ╔═╡ 33ec9eea-ef8e-4cd3-ad0d-ca340a1613c4
begin
	pars_KG_G = select(post, # Posteriors for the Killifish-Guppy fish
		:Intercept_survG => :"α_surv",
    	:b_z_survG => :"βz_surv",
    	:Intercept_growG => :"α_grow",
    	:b_z_growG => :"βz_grow",
    	:sigma_growG => :σ_grow,
    	:Intercept_recrG => :"α_fec",
    	:b_z_recrG => :"βz_fec"
	)		
	pars_NK_G = DataFrame(  # posteriors for the No-Killifish fish
		α_surv = post.Intercept_survG .+ post.b_NK_survG, 
		βz_surv = post.b_z_survG .+ post.b_zNK_survG,
		α_grow = post.Intercept_growG .+ post.b_NK_growG,
		βz_grow= post.b_z_growG .+ post.b_zNK_growG ,
		σ_grow = post.sigma_growG,
		α_fec = post.Intercept_recrG .+ post.b_NK_recrG,
		βz_fec = post.b_z_recrG + post.b_zNK_recrG
	)

end;

# ╔═╡ 00bbac7d-8da9-4d3f-8937-3b581f7e9d9d
md"""
#### IPM for KG and NK guppies
In the next lines, we will build an IPM model for guppies
"""

# ╔═╡ 9be4047b-e7fa-4872-9a2c-a40986a1db65
# Here, we choose the first posterior values to make and IMP and test our code
# we center our analyses at 18.0 mm
begin
		size_cen = 18.0
		row = 1
		IPM_KG = mk_K(pars_KG_G, z1, z, size_cen, row)
		IPM_NK = mk_K(pars_NK_G, z1, z, size_cen, row)
end;

# ╔═╡ c54efaad-e70d-4f41-ba82-e1dca6987afe
# see the kernel (K)
IPM_KG["K"]

# ╔═╡ 953d638a-f757-4b6a-9e78-eea9f4ae9729
md"""
After constructing the IMP, we can estimate some interesting values. Such as the asymptotic growth rate (λ), our measurement of fitness.
"""

# ╔═╡ b35fc72e-2b81-4fc2-8657-ff375beaf906
begin
	# calculate the population growth rate (λ)
	vv_KG = eigen(IPM_KG["K"])
	vv_NK = eigen(IPM_NK["K"])
	λ_KG = real(vv_KG.values[end])
	λ_NK = real(vv_NK.values[end])
end;

# ╔═╡ 8643c13c-1550-4622-8b4d-d51577c4fb7a
# let's check our results
[λ_KG,λ_NK]

# ╔═╡ 3608e3a1-de75-4c74-a92d-e9f676889b2d
md"""
This suggest that the fitness of KG guppies is slightly higher than that of NK guppies.
"""

# ╔═╡ 02879703-37d6-45f8-8fef-b90f6da4e157
md"""
### Life Table Response Experiment (LTRE) Analyses

We used LTRE analyses to decompose the fitness contrasts into effects caused by the individual demographic rates (Caswell 2001). The first-order approximation of the decomposition for each contrast for each demographic rate is:

$\lambda_V \approx  (V_i - V_j) \frac{d\lambda}{d^T\bar{V}}$ 

where, $e^T$ is a vector of ones, $vec V_i$ and $vec V_j$ are the vector transformed demographic rates (e.g., S(z)) for treatments $i$ and $j$. $\frac{d\lambda}{d vec^T \bar{V}}$ is the sensitivity of absolute fitness (λ) to the vector transformed matrix. The sensitivity of guppy fitness to the demographic rates in the guppy model was calculated as: 

$\frac{d\lambda}{d ^T \bar{V}}= \frac{d\lambda}{d^T A} \frac{d A}{d^T V}$

where $\frac{d\lambda}{d^T A}$ is the vector transformed sensitivity of absolute fitness to the matrix elements of $A$. $\frac{d A}{d ^T V}$ is the sensitivity of the elements of A to each demographic rate (V). 

"""

# ╔═╡ 3a6bf87d-bd2a-412c-900b-34f72b79b83d
## LTRE
begin		
	## calculate average matrix
	K_avg = (IPM_NK["K"] + IPM_KG["K"])./2
	vv_avg = eigen(K_avg)
	
	# Normalize stable size distribution
	ω_avg = real(vv_avg.vectors[:, end]) 
	
	# Reproductive value
	a_avg = eigen(transpose(K_avg))
	v_avg = real(a_avg.vectors[:, end]) ./ Base.sum(real(a_avg.vectors[:, end]))
	v_avg = v_avg ./ dot(transpose(v_avg), ω_avg)
	
	## Sensitivity matrix
	sens_avg = v_avg * ω_avg' # this equivalent to outer() in R
	ΔK = IPM_NK["K"] - IPM_KG["K"]
	λ_eff = ΔK .* sens_avg
	Δλ = Base.sum(λ_eff)
	[Base.sum(λ_eff),Δλ]
end

# ╔═╡ 5224d66c-7210-409d-bb1c-884f1e6b7e3f
begin
	## Make life-response table
	one_mat = ones(nBigMatrix, nBigMatrix)
	
	# Function differences
	Δ_grow = g_z1z(pars_NK_G, z1, z, size_cen, row) - g_z1z(pars_KG_G, z1, z, size_cen, row)
	Δ_fec = one_mat*(R_z(pars_NK_G, z, size_cen, row) - R_z(pars_KG_G, z, size_cen, row))
	Δ_rcz = (c_z1z(pars_NK_G, z1, z, size_cen, row) - c_z1z(pars_KG_G, z1, z, size_cen, row))
	Δ_sur = one_mat*(s_z(pars_NK_G, z, size_cen, row) - s_z(pars_KG_G, z, size_cen, row))
	
	# Function averages
	grow_avg = (g_z1z(pars_NK_G, z1, z, size_cen, row) + g_z1z(pars_KG_G, z1, z, size_cen, row))/2
	fec_avg = (one_mat*(R_z(pars_NK_G, z, size_cen, row) + R_z(pars_KG_G, z, size_cen, row)))/2
	rcz_avg = ((c_z1z(pars_NK_G, z1, z, size_cen, row) + c_z1z(pars_KG_G, z1, z, size_cen, row)))/2
	sur_avg = (one_mat*(s_z(pars_NK_G, z, size_cen, row) + s_z(pars_KG_G, z, size_cen, row)))/2
	
	# derivates
	𝛿_grow= sur_avg
	𝛿_fec=  rcz_avg .* sur_avg
	𝛿_sur = grow_avg .+   fec_avg .* rcz_avg
	𝛿_rcz = fec_avg .* sur_avg
	
	λ_grow = Δ_grow .* sens_avg .* 𝛿_grow
	λ_fec = Δ_fec .* sens_avg .* 𝛿_fec
	λ_rcz = Δ_rcz .* sens_avg .* 𝛿_rcz
	λ_sur = Δ_sur .* sens_avg .* 𝛿_sur
	
	# to put in a DataFrame
	sur_con = Base.sum(λ_sur)
	grow_con = Base.sum(λ_grow)
	fec_con = Base.sum(λ_fec)
	rcz_con = Base.sum(λ_rcz)
	sum_con = sur_con + grow_con + fec_con + rcz_con

end


# ╔═╡ fca2ea39-40b7-42a7-9cfe-b1b9aeb8e2f0
md"""
### Integral Projection model for killifish
The model for killifish is also a single-sex (female only) IPM but is slightly different than the one for guppies. Because killifish lay eggs that hatch at two weeks, two broods are produced per time-step. This means that the overall projection kernel for a month is a composite of two, two-week projection kernels. Each two-week projection equation is: 

$n(t+1/2)_E = \int R(z)S(z)n(z,t)_A dz,$

$n(z, t+1/2)_A = \int C(z',z)vn(t)_E + \int G(z',z)S(z)n(z,t)_A dz,$

the subscripts $E$ and $A$ denote which equation describes the production of eggs or the transition from eggs to adults or adults to adults. $v$ is the probability of egg survival and hatching, here set to 0.75 for all treatments as we do not have information about this transition. Additionally, because killifish have a larger size range than guppies we set the L and U size range between 2 mm and 110 mm.  The composite functions for the entire 28 days are then:

$n(t+1)_E = \int R(z)S(z)[C(z',z) v n(t)_E + \int G(z',z)S(z)n(z,t)_A dz ]dz,$

$n(z, t+1)_A = \int C(z',z) v  \int R(z)S(z)n(z,t)_A dz + \int G(z',z)S(z)[C(z',z)S(z) v n(t)_E + \int G(z',z)S(z)n(z,t)_A dz]dz,$
"""

# ╔═╡ f0afef35-7cab-467b-a6a1-33947e6f2faa
begin
  # nBigMatrix = 100
  Uk= 110.0
  Lk= 2.0
  mk = nBigMatrix
  hk = (Uk - Lk)/mk

  z1k =  Lk .+ (collect(1:mk) .- 0.5) .* hk
  z1k = round.(z1k, digits = 6)
  zk = z1k
  meshptsk = z1k
end

# ╔═╡ a8227963-838b-4a92-ae8b-9995440b7ed4
begin
	pars_KG_K = select(post, # Posteriors for the Killifish-Guppy fish
		:Intercept_survK => :"α_surv",
    	:b_z_survK => :"βz_surv",
    	:Intercept_growK => :"α_grow",
    	:b_z_growK => :"βz_grow",
    	:sigma_growK => :σ_grow,
    	:Intercept_recrK => :"α_fec",
    	:b_z_recrK => :"βz_fec"
	)		
	pars_NG_K = DataFrame(  # posteriors for the No-Killifish fish
		α_surv = post.Intercept_survK .+ post.b_NG_survK, 
		βz_surv = post.b_z_survK .+ post.b_zNG_survK,
		α_grow = post.Intercept_growK .+ post.b_NG_growK,
		βz_grow= post.b_z_growK .+ post.b_zNG_growK ,
		σ_grow = post.sigma_growK,
		α_fec = post.Intercept_recrK .+ post.b_NG_recrK,
		βz_fec = post.b_z_recrK + post.b_zNG_recrK
	)

end;

# ╔═╡ 15538fe8-3248-4832-8fbd-671e781a6585
md"""
To acomodate the egg stage in the killfish matrices, we added an extra column and extra row. See the structure of each matrix.
"""

# ╔═╡ cae2a45d-debf-46c9-a710-e9194f1b4bd4
function g_z1zK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
    α= df.α_grow[row]
    β= df.βz_grow[row]
    σ= df.σ_grow[row]
    p_den_grow = zeros(size(z)[1],size(z)[1])
    μ = ((α .+ β .* (z .- size_cen ) .- (z))./2 .+ z) # average growth in two weeks
    for i in 1:length(z)
		p_den_grow[:,i] = pdf.(Normal(μ[i], σ), z1).*h
    end
    matex = zeros(length(z)+1, length(z)+1)
    matex[2:end,2:end] = p_den_grow
    return(matex)
  end

# ╔═╡ 9ef5ad8b-d034-4e13-af2b-025b2fc946f8
g_z1zK(pars_KG_K, z1k, zk, 18.0, 1)

# ╔═╡ 6a835c44-dfd4-4138-9e9c-a22b37ea4b2c
## Surival function
function s_zK(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, row::Integer)
    α= df.α_surv[row]
    β= df.βz_surv[row]
    linear_p = α .+ β * (z .- size_cen)      # linear predictor
    p = 1 ./(1 .+ exp.(-linear_p))
    p = diagm(sqrt.(p))
    matex = zeros(length(z)+1, length(z)+1)
    matex[2:end,2:end] = p
    return(matex)
end

# ╔═╡ a03b064c-da03-4ef6-ba0b-5bca3eed79df
s_zK(pars_KG_K, zk, 18.0, 1)

# ╔═╡ 28eab40c-b147-4905-b743-0fd31df65553
## Recruitment function (N.B - from birth in spring to first summer), logistic regression
function R_zK(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
    α= df.α_fec[row]
    β= df.βz_fec[row]
    linear_p = α .+ β * (z .- size_cen)   # linear predictor  # linear predictor
    p = exp.(linear_p).*(1/2)
    matex = zeros(length(z)+1, length(z)+1)
    matex[1,2:end] = p
    return(matex)
end

# ╔═╡ 522006b1-7c16-4439-a335-f312ea22a2ab
R_zK(pars_KG_K, zk, 18.0, 1)

# ╔═╡ be507426-cf42-4665-bd6d-3bb214b8f921
## Recruit size function
function c_z1zK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
    size_cen::AbstractFloat, row::Integer)
	
	σ= 0.6 #pars.rcz_sd[row]
    p_den_grow = zeros(length(z)+1,length(z)+1)
    μ = 4.4 .+ 0#*(z .- size_cen)
    p_den_grow[2:end, 1] =  pdf.(Normal(μ, σ), z).*h
    return(p_den_grow)
end

# ╔═╡ 41285ce5-682d-48f4-bd34-b8caadc6b79d
c_z1zK(pars_KG_K,z, z, 18.0, 1)

# ╔═╡ 3d685bd5-39d1-4700-91c9-4ce7b7ff60d8
## Functions to build IPM kernels P, F, and K
function mk_KK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, size_cen::AbstractFloat, row::Integer, V::AbstractFloat)
	F = R_zK(df, z, size_cen, row) * s_zK(df, z, size_cen, row)
	P = g_z1zK(df, z, z, size_cen, row) * s_zK(df, z, size_cen, row)
	A = (P + F) + (V .* c_z1zK(df, z1, z, size_cen, row))
	K = A * A
	out = Dict("K"=> K, "meshpts" => z, "P" => P, "F"=> F, "A" => A)
	return(out)
end

# ╔═╡ f7e048f3-bf0e-488a-bdf8-ae8aa0c13d6c
md"""
#### IPM for KG and NK killifish
In the next lines, we will build an IPM model for killifish
"""

# ╔═╡ af4ef0da-7d00-452c-aa6d-ad663a14ea69
begin
# Make projection kernels
	IPM_GR = mk_KK(pars_KG_K, z1, z, size_cen, row, 0.7)
	IPM_NG = mk_KK(pars_NG_K, z1, z, size_cen, row, 0.7)
	# calculate the population growth rate (λ)
	vv_GR = eigen(IPM_GR["K"])
	vv_NG = eigen(IPM_NG["K"])
	λ_GR = real(vv_GR.values[end])
	λ_NG = real(vv_NG.values[end])
	[λ_NG, λ_GR, λ_NG-λ_GR] # fitness differences
end

# ╔═╡ 0d11ad15-5c5e-43f5-b2eb-2560e5785419
md"""
### Life Table Response Experiment (LTRE) Analyses for killifish

The first-order approximation of the decomposition for each contrast for each demographic rate is:

$\lambda_V \approx  e^T (vec V_i - vec V_j) \frac{d\lambda}{d vec ^T\bar{V}}$ 

where, $e^T$ is a vector of ones and $vec$ represents  vector transformed matrices  (e.g., S(z)). The sensitivity for killifish requires an extra derivate to relate how changes in $A$ alters $\bar{A}$. This derivate was calculated:

$\frac{d vec \bar{A}}{d vec^T A} = (I \otimes A) + (\bar{A}^T \otimes I)$

where $I$ is an identity matrix and $otimes$ means to take the Kronecker product. This derivate is inserted into the sensitivity equation:

$\frac{d\lambda}{d vec^T \bar{V}}= \frac{d\lambda}{d vec^T \bar{A}} \frac{d vec \bar{A}}{d vec^T A} \frac{d vec A}{d vec^T V}$


"""

# ╔═╡ 07a56b90-33e9-428b-9f13-e1658e0f0ba7
# ## Killifish LTRE
# begin
# 	## calculate average matrix
# 	K_avg = (IPM_NG["K"] + IPM_GR["K"])./2
# 	vv_avg = eigen(K_avg)
	
# 	# Normalize stable size distribution
# 	ω_GR = real(vv_GR.vectors[:, end]) ./ Base.sum(real(vv_GR.vectors[:, end]))
# 	ω_NG = real(vv_NG.vectors[:, end]) ./  Base.sum(real(vv_NG.vectors[:, end]))
# 	ω_avg = real(vv_avg.vectors[:, end]) ./ Base.sum(real(vv_avg.vectors[:, end]))
	
# 	# Reproductive value
# 	a_avg = eigen(K_avg')
# 	v_avg = real(a_avg.vectors[:, end]) 
# 	v_avg = v_avg ./ dot(transpose(v_avg), ω_avg)
	
# 	## Sensitivity matrix
# 	sens_avg = v_avg * ω_avg'  # this equivalent to outer() in R
# 	ΔK = IPM_NG["K"] .- IPM_GR["K"]
# 	λ_eff = ΔK .* sens_avg
# 	Δλ = sum(λ_eff)
	    
# 	## Make life-response table
	
# 	# Function differences
# 	Δ_grow = g_z1zK(pars_NG_K, z1, z, size_cen, row) .- g_z1zK(pars_GR_K, z1, z, size_cen, row)
# 	Δ_fec = (pr_zK(pars_NG_K, z, size_cen, row) .- pr_zK(pars_GR_K, z, size_cen, row))
# 	    Δ_rcz = (c_z1zK(pars_NG_K, z1, z, size_cen, row) .- c_z1zK(pars_GR_K, z1, z, size_cen, row))
# 	Δ_sur = (s_zK(pars_NG_K, z, size_cen, row) .- s_zK(pars_GR_K, z, size_cen, row))
	
# 	# Function averages
# 	grow_avg = (g_z1zK(pars_NG_K, z1, z, size_cen, row) + g_z1zK(pars_GR_K, z1, z, size_cen, row))/2
	
# 	fec_avg = ((pr_zK(pars_NG_K, z, size_cen, row) + pr_zK(pars_GR_K, z, size_cen, row)))/2
	
# 	rcz_avg = ((c_z1zK(pars_NG_K, z1, z, size_cen, row) + c_z1zK(pars_GR_K, z1, z, size_cen, row)))/2
	
# 	sur_avg = ((s_zK(pars_NG_K, z, size_cen, row) + s_zK(pars_GR_K, z, size_cen, row)))/2
	
# 	# derivates
# 	I_mat = diagm(ones(nBigMatrix+1))
# 	𝛿_grow = kron(transpose(sur_avg), I_mat)
# 	𝛿_sur  = kron(transpose(I_mat),grow_avg) +  kron(transpose(I_mat),(fec_avg * rcz_avg))
# 	𝛿_fec  = kron(transpose(sur_avg), rcz_avg)
# 	𝛿_rcz  = kron(transpose(fec_avg * sur_avg), I_mat)
	
# 	𝛿A = kron(I_mat, K_avg) + kron(transpose(K_avg), I_mat)
	
# 	λ_grow = Δ_grow .* reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_grow), (nBigMatrix+1,nBigMatrix+1))
	
# 	λ_fec = Δ_fec .*  reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_fec), (nBigMatrix+1,nBigMatrix+1))
	
# 	λ_rcz = Δ_rcz .*  reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_rcz), (nBigMatrix+1,nBigMatrix+1)) 
	
# 	λ_sur = Δ_sur .* reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_sur), (nBigMatrix+1,nBigMatrix+1)) 
	
	
	
# 	# to put in a DataFrame
# 	sur_con = Base.sum(λ_sur)
# 	grow_con = Base.sum(λ_grow)
# 	fec_con = Base.sum(λ_fec)
# 	rcz_con = Base.sum(λ_rcz)
# 	sum_con = sur_con + grow_con + fec_con + rcz_con

# end

# ╔═╡ ce9b44ee-3a44-493b-ab2d-51335a7d0fd3
md"""
# Putting things together

To use the IPMs and make them meaningful, we have to run them with each iteration from the posterior samples. To facilitate that process you only have to run the *".../Julia/MainScript.jl"* file. This script will call other scripts where the IMP for guppies and killifish have been compiled into functions. Thes analyses will basically run the analyses 6000 times. The *"MainScript.jl"* file will also produce the figures and estimate all the necessary stats. Writing these IPMs in the Julia language has been easy and straightforward, but the major advantage is its speed. If you prefer to run the analyses in "R" you can also do so, but running the script "R/MainScript.R" the results are basically the same, there are some small differences but it is more because of the differences in the decimal digits use by the two programing lenguages. 

Running the killifish IPM in R took around 3 days on my computer, during that time I translated the code to Julia. Then, in Julia, the same code run in 4 hours !!! 




"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
CSV = "~0.10.3"
DataFrames = "~1.3.2"
Distributions = "~0.25.52"
JLD2 = "~0.4.22"
RCall = "~0.13.13"
StatsPlots = "~0.14.33"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "91ca22c4b8437da89b030f08d71db55a379ce958"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.3"

[[Arpack_jll]]
deps = ["Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "e214a9b9bd1b4e1b4f15b22c0994862b66af7ff7"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.0+3"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "9310d9495c1eb2e4fa1955dd478660e2ecab1fbb"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.3"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "109664d3a6f2202b1225478335ea8fea3cd8706b"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.5"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "c43e992f186abaf9965cc45e372f4693b7754b22"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.52"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "4f00cc36fede3c04b8acf9b2e2763decfdcecfa6"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.13"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "7008a3412d823e29d370ddc77411d593bd8a3d03"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.9.1"

[[NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "16baacfdc8758bc374882566c9187e785e85c2f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.9"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e8185b83b9fc56eb6456200e873ce598ebc7f262"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.7"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "90021b03a38f1ae9dbd7bf4dc5e3dcb7676d302c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.2"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "28ef6c7ce353f0b35d0df0d5930e0d072c1f5b9b"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[RCall]]
deps = ["CategoricalArrays", "Conda", "DataFrames", "DataStructures", "Dates", "Libdl", "Missings", "REPL", "Random", "Requires", "StatsModels", "WinReg"]
git-tree-sha1 = "72fddd643785ec1f36581cbc3d288529b96e99a7"
uuid = "6f49c342-dc21-5d91-9882-a32aef131414"
version = "0.13.13"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "6976fab022fea2ffea3d945159317556e5dad87c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.2"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "03c99c7ef267c8526953cafe3c4239656693b8ab"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.29"

[[StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "4d9c69d65f1b270ad092de0abe13e859b8c55cad"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.14.33"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "505c31f585405fc375d99d02588f6ceaba791241"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.5"

[[WinReg]]
deps = ["Test"]
git-tree-sha1 = "808380e0a0483e134081cc54150be4177959b5f4"
uuid = "1b915085-20d7-51cf-bf83-8f477d6f5128"
version = "0.3.1"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─283e2bac-aabf-11ec-1e3e-ed84e95588ca
# ╟─b6c40268-0a24-4691-aff9-a487eb9846e8
# ╟─89453239-97e9-4d6a-89a0-905118574bf4
# ╟─dc2034eb-cc01-4032-824b-9e0ace72e329
# ╟─560f7b39-731b-4a8e-a16b-b822d551073f
# ╠═5736dcb2-0be4-4be8-9fcb-21d679caf199
# ╟─162e9c0c-e0f3-4e1d-9a4e-a9d3694d361b
# ╠═b12084bb-8040-40d1-9657-2800e5b6ccb2
# ╟─6bb7c116-8c16-4df6-8f34-501f652faa12
# ╠═ce97c58b-c9c5-4aa3-962a-7e19e1fb26c1
# ╠═ad66cf19-db0f-4f54-9c48-deb72543e7ad
# ╟─02227c46-5a15-484d-a302-0f2a4a51c749
# ╟─ed107163-cd0d-4f5b-8840-3859f2363cf9
# ╟─95769a16-6e17-4d10-92d4-ec38f82f7aad
# ╠═1c358e1c-a072-4273-bf7f-9d24ab84ea46
# ╟─a4bcbf5d-2cdf-444c-9c1f-2acd8151250a
# ╠═87d2c9f8-b26c-4c5e-9f59-33c7f73971fc
# ╟─c60cd474-89af-439c-b457-b302ccb1e2d8
# ╠═d3fedac0-0a07-4f2d-bda9-ecc5643e23b9
# ╟─bf727831-4cef-4e7e-adda-6c47ab7c6ac2
# ╠═9d241282-e6f2-4d63-8cf2-6d49dde03ac0
# ╟─fed2e2fd-bce9-46ff-8f48-360a86136209
# ╠═f504b0a1-fb69-4069-a8ca-76e815438c5d
# ╠═e242e6dc-f7f2-4663-a893-a6bfbcdafd8f
# ╟─7b44a540-00b2-4db9-91e4-2e782e5b77d4
# ╠═13db1664-0f55-4e32-b59c-6f25309e93b6
# ╠═b58511e7-0703-4e60-8119-7c9a8e01fbd1
# ╟─fe287182-f606-43f3-8c30-def2d342af54
# ╟─5467e312-3522-4343-9fcb-fe1c3fe925dd
# ╠═33ec9eea-ef8e-4cd3-ad0d-ca340a1613c4
# ╟─00bbac7d-8da9-4d3f-8937-3b581f7e9d9d
# ╠═9be4047b-e7fa-4872-9a2c-a40986a1db65
# ╠═c54efaad-e70d-4f41-ba82-e1dca6987afe
# ╟─953d638a-f757-4b6a-9e78-eea9f4ae9729
# ╠═b35fc72e-2b81-4fc2-8657-ff375beaf906
# ╠═8643c13c-1550-4622-8b4d-d51577c4fb7a
# ╟─3608e3a1-de75-4c74-a92d-e9f676889b2d
# ╟─02879703-37d6-45f8-8fef-b90f6da4e157
# ╠═3a6bf87d-bd2a-412c-900b-34f72b79b83d
# ╠═5224d66c-7210-409d-bb1c-884f1e6b7e3f
# ╠═fca2ea39-40b7-42a7-9cfe-b1b9aeb8e2f0
# ╠═f0afef35-7cab-467b-a6a1-33947e6f2faa
# ╠═a8227963-838b-4a92-ae8b-9995440b7ed4
# ╟─15538fe8-3248-4832-8fbd-671e781a6585
# ╠═cae2a45d-debf-46c9-a710-e9194f1b4bd4
# ╠═9ef5ad8b-d034-4e13-af2b-025b2fc946f8
# ╠═6a835c44-dfd4-4138-9e9c-a22b37ea4b2c
# ╠═a03b064c-da03-4ef6-ba0b-5bca3eed79df
# ╠═28eab40c-b147-4905-b743-0fd31df65553
# ╠═522006b1-7c16-4439-a335-f312ea22a2ab
# ╠═be507426-cf42-4665-bd6d-3bb214b8f921
# ╠═41285ce5-682d-48f4-bd34-b8caadc6b79d
# ╠═3d685bd5-39d1-4700-91c9-4ce7b7ff60d8
# ╟─f7e048f3-bf0e-488a-bdf8-ae8aa0c13d6c
# ╠═af4ef0da-7d00-452c-aa6d-ad663a14ea69
# ╟─0d11ad15-5c5e-43f5-b2eb-2560e5785419
# ╠═07a56b90-33e9-428b-9f13-e1658e0f0ba7
# ╟─ce9b44ee-3a44-493b-ab2d-51335a7d0fd3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
