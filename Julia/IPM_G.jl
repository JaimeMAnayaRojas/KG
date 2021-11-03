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
  




function Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0)
	
	# nBigMatrix = 100
	# min_size = 4
	# max_size = 35
	# size_cen = 18.0
	# area = 0.0
	# canopy =0.0

	U= Float64(max_size)
	L=Float64(min_size)
	m = nBigMatrix
	h = (U - L)/m
	z1 = zeros(nBigMatrix)
	z1 =  L .+ (collect(1:m) .- 0.5) * h


	z = z1
	meshpts = z1

	pars_GR_G = select(post,
    	:Intercept_survG => :"α_surv",
    	:b_z_survG => :"βz_surv",
    	:b_area_survG => :"β_area_surv",
    	:b_canopy_survG => :"β_canopy_surv",

    	:Intercept_growG => :"α_grow",
    	:b_z_growG => :"βz_grow",
    	:b_area_growG => :"β_area_grow",
    	:b_canopy_growG => :"β_canopy_grow",
		:sigma_growG => :σ_grow,

    	:Intercept_recrG => :"α_fec",
    	:b_z_recrG => :"βz_fec",
    	:b_area_recrG => :"β_area_fec",
    	:b_canopy_recrG => :"β_canopy_fec"
	)		


	pars_NR_G = select(post,
		:Intercept_survG => :"α_surv",
		:b_z_survG => :"βz_surv",
		:b_area_survG => :"β_area_surv",
		:b_canopy_survG => :"β_canopy_surv",

		:Intercept_growG => :"α_grow",
		:b_z_growG => :"βz_grow",
		:b_area_growG => :"β_area_grow",
		:b_canopy_growG => :"β_canopy_grow",
		:sigma_growG => :σ_grow,

		:Intercept_recrG => :"α_fec",
		:b_z_recrG => :"βz_fec",
		:b_area_recrG => :"β_area_fec",
		:b_canopy_recrG => :"β_canopy_fec"


	)


	pars_NR_G.α_surv = pars_NR_G.α_surv .+ post.b_NK_survG
	pars_NR_G.βz_surv = pars_NR_G.βz_surv .+ post.b_zNK_survG 


	pars_NR_G.α_grow = pars_NR_G.α_grow .+ post.b_NK_growG
	pars_NR_G.βz_grow = pars_NR_G.βz_grow .+ post.b_zNK_growG 


	pars_NR_G.α_fec = pars_NR_G.α_fec .+ post.b_NK_recrG
	pars_NR_G.βz_fec = pars_NR_G.βz_fec .+ post.b_zNK_recrG 


	df = nothing
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


	#@time g = g_z1z(pars_GR_G, z1, z, size_cen, row, 0.0, 0.0)

	# columns should sum to 1
	#sum.(eachcol(g))


	## Surival function
	#row = 1
	function s_z(df::AbstractDataFrame, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		α= df.α_surv[row]
		β= df.βz_surv[row]
		
		linear_p = α .+ β * (z .- size_cen)  # linear predictor
		p = 1 ./(1 .+ exp.(-linear_p))
		p = diagm(p)
		return(p)
	end

	#s_z(pars_GR_G, z, size_cen, 1, 0.0, 0.)

	## Reproduction function, logistic regression
	# function pb_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
	#   α= df.α_rep[row]
	#   β= df.βz_rep[row]
	#   linear_p = α .+ β * (z .- size_cen)       # linear predictor
	#   p = 1 ./(1 .+ exp.(-linear_p))
	#   p = diagm(p)
	#   return(p)
	# end



	## Recruitment function (N.B - from birth in spring to first summer), logistic regression
	function pr_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
		α= df.α_fec[row]
		β= df.βz_fec[row]
		

		linear_p = α .+ β * (z .- size_cen)     # linear predictor
		p = exp.(linear_p)*(1/2)
		p = diagm(p)
		return(p)
	end


	## Recruit size function
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


	##----------------------------------------------------
	## Functions to build IPM kernels P, F, and K

	function P_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)

		out = g_z1z(df, z1, z, size_cen, row) * s_z(df, z, size_cen, row)
		return(out)

	end

	function F_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		out1 = c_z1z(df, z1, z, size_cen, row)
		out2 = pr_z(df, z, size_cen, row) * s_z(df, z, size_cen, row)
		out = out1 * out2
		return(out)
	
	end


	function mk_K(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		F = F_z1z(df, z1, z, size_cen, row)
		P = P_z1z(df, z1, z, size_cen, row)
		K = P + F
		out = Dict("K"=> K, "meshpts" => z, "P" => P, "F"=> F)
		return(out)
	end

	#mat = IPM_GR["K"]



	Gsurv_mat = zeros(size(post)[1], nBigMatrix)
	Ggrow_mat = zeros(size(post)[1], nBigMatrix)
	#Grep_mat = zeros(size(pars_GR_G)[1], nBigMatrix)
	Gfec_mat = zeros(size(post)[1], nBigMatrix)
	Grcz_mat = zeros(size(post)[1], nBigMatrix)

	## make DataFrame to store the results


	Gres_IPM = DataFrame(zeros(size(post)[1], 15), :auto)
	Gres_IPM = select(Gres_IPM, :x1 => "lam_GR", :x2 => "lam_NR", :x3 => "delta_lam",
						:x4 => "sum_lam_eff", :x5 => "grow_con", :x6 => "fec_con", 
						:x7 => "rcz_con", :x8 => "sur_con",
						:x9 => "sum_con")

    
	for row in 1:size(post)[1]
		# Make projection kernels
		IPM_GR = mk_K(pars_GR_G, z1, z, size_cen, row)
		IPM_NR = mk_K(pars_NR_G, z1, z, size_cen, row)
		# calculate the population growth rate (λ)

		vv_GR = eigen(IPM_GR["K"])
		vv_NR = eigen(IPM_NR["K"])

		λ_GR = real(vv_GR.values[end])
		λ_NR = real(vv_NR.values[end])
		Gres_IPM.lam_GR[row] = λ_GR
		Gres_IPM.lam_NR[row] = λ_NR

		## calculate average matrix
		K_avg = (IPM_NR["K"] + IPM_GR["K"])./2
		vv_avg = eigen(K_avg)

		# Normalize stable size distribution
		ω_GR = real(vv_GR.vectors[:, end]) 
		ω_NR = real(vv_NR.vectors[:, end]) 
		ω_avg = real(vv_avg.vectors[:, end]) 

		# Reproductive value
		a_GR = eigen(transpose(IPM_GR["K"]))
		a_NR = eigen(transpose(IPM_NR["K"]))
		a_avg = eigen(transpose(K_avg))

		v_GR = real(a_GR.vectors[:, end]) 
		v_NR = real(a_NR.vectors[:, end])
		v_avg = real(a_avg.vectors[:, end]) / sum(real(a_avg.vectors[:, end]))
		v_avg = v_avg / dot(transpose(v_avg), ω_avg)

		## Sensitivity matrix

		sens_avg = v_avg * ω_avg' # this equivalent to outer() in R
		ΔK = IPM_NR["K"] - IPM_GR["K"]
		λ_eff = ΔK .* sens_avg
		Δλ = sum(λ_eff)
		

		Gres_IPM.sum_lam_eff[row] = Δλ
		Gres_IPM.delta_lam[row] = λ_NR - λ_GR

		## Make life-response table

		one_mat = ones(nBigMatrix, nBigMatrix)

		# Function differences
		Δ_grow = g_z1z(pars_NR_G, z1, z, size_cen, row) - g_z1z(pars_GR_G, z1, z, size_cen, row)
		#Δ_rep = one_mat*(pb_z(pars_NR_G, z, size_cen, row) - pb_z(pars_GR_G, z, size_cen, row))
		Δ_fec = one_mat*(pr_z(pars_NR_G, z, size_cen, row) - pr_z(pars_GR_G, z, size_cen, row))
		Δ_rcz = (c_z1z(pars_NR_G, z1, z, size_cen, row) - c_z1z(pars_GR_G, z1, z, size_cen, row))
		Δ_sur = one_mat*(s_z(pars_NR_G, z, size_cen, row) - s_z(pars_GR_G, z, size_cen, row))

		# Function averages
		grow_avg = (g_z1z(pars_NR_G, z1, z, size_cen, row) + g_z1z(pars_GR_G, z1, z, size_cen, row))/2
		#rep_avg = (one_mat*(pb_z(pars_NR_G, z, size_cen, row) + pb_z(pars_GR_G, z, size_cen, row)))/2
		fec_avg = (one_mat*(pr_z(pars_NR_G, z, size_cen, row) + pr_z(pars_GR_G, z, size_cen, row)))/2
		rcz_avg = ((c_z1z(pars_NR_G, z1, z, size_cen, row) + c_z1z(pars_GR_G, z1, z, size_cen, row)))/2
		sur_avg = (one_mat*(s_z(pars_NR_G, z, size_cen, row) + s_z(pars_GR_G, z, size_cen, row)))/2

		# derivates

		# 𝛿_grow= sur_avg
		# 𝛿_rep= fec_avg*rcz_avg*sur_avg
		# 𝛿_fec= rep_avg * fec_avg * sur_avg
		# 𝛿_sur = grow_avg +  rep_avg * fec_avg * rcz_avg
		# 𝛿_rcz = rep_avg * fec_avg * sur_avg

		# λ_grow = Δ_grow .* sens_avg .* 𝛿_grow
		# λ_rep = Δ_rep .* sens_avg .* 𝛿_rep
		# λ_fec = Δ_fec .* sens_avg .* 𝛿_fec
		# λ_rcz = Δ_rcz .* sens_avg .* 𝛿_rcz
		# λ_sur = Δ_sur .* sens_avg .* 𝛿_sur

		𝛿_grow= sur_avg
		𝛿_fec=  rcz_avg .* sur_avg
		𝛿_sur = grow_avg .+   fec_avg .* rcz_avg
		𝛿_rcz = fec_avg .* sur_avg

		λ_grow = Δ_grow .* sens_avg .* 𝛿_grow
		λ_fec = Δ_fec .* sens_avg .* 𝛿_fec
		λ_rcz = Δ_rcz .* sens_avg .* 𝛿_rcz
		λ_sur = Δ_sur .* sens_avg .* 𝛿_sur



		# to put in a DataFrame
		sur_con = sum(λ_sur)
		grow_con = sum(λ_grow)
		#rep_con = sum(λ_rep)
		fec_con = sum(λ_fec)
		rcz_con = sum(λ_rcz)
		sum_con = sur_con + grow_con + fec_con + rcz_con

		Gres_IPM.sum_con[row] = sum_con
		Gres_IPM.sur_con[row] = sur_con
		Gres_IPM.grow_con[row] = grow_con
		#Gres_IPM.rep_con[row] = rep_con
		Gres_IPM.fec_con[row] = fec_con
		Gres_IPM.rcz_con[row] = rcz_con

		Gsurv_mat[row, : ] =  sum.(eachcol(λ_sur)) 
		Ggrow_mat[row, : ] =  sum.(eachcol(λ_grow)) 
		#rep_mat[row, : ] =  sum.(eachcol(λ_rep)) 
		Gfec_mat[row, : ] =  sum.(eachcol(λ_fec)) 
		Grcz_mat[row, : ] =  sum.(eachcol(λ_rcz)) 

	end

	return [Gres_IPM, Gsurv_mat, Ggrow_mat, Gfec_mat, Grcz_mat]

end


