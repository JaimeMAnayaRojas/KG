## Funcitons for analyses


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


function my_summary(df, digits = 3)
    vm= DataFrame(df, :auto)
    ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
    ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
    ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))
    df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=digits),
    l99 = round.(Vector(ci99[1, :]), digits =digits),
    l95 = round.(Vector(ci95[1, :]), digits =digits), 
    l68 = round.(Vector(ci68[1, :]), digits =digits),
    u68 = round.(Vector(ci68[2, :]), digits =digits),
    u95 = round.(Vector(ci95[2, :]), digits =digits),
    u99 = round.(Vector(ci99[2, :]), digits =digits))

    return df
end


function My_Logit(z)
    l = log(z / (1 - (z)))
    return l 
end 




function My_Logistic(x)
    l = 1/(1+exp(-x))
    return l 
end 
  
function Guppy_IPM(post; nBigMatrix = 100, min_size = 4, max_size = 35, size_cen = 18.0)
	
	# nBigMatrix = 100
	# min_size = 4
	# max_size = 35
	# size_cen = 18.0
	# BiomassK = 0.0
	# canopy =0.0

	U= Float64(max_size)
	L=Float64(min_size)
	m = nBigMatrix
	h = (U - L)/m
	z1 = zeros(nBigMatrix)
	z1 =  L .+ (collect(1:m) .- 0.5) * h
	z = z1
	

	pars_KG_G = select(post,
    	:Intercept_survG => :"??_surv",
    	:b_z_survG => :"??z_surv",
    	:Intercept_growG => :"??_grow",
    	:b_z_growG => :"??z_grow",
    	:sigma_growG => :??_grow,
    	:Intercept_recrG => :"??_fec",
    	:b_z_recrG => :"??z_fec",
	)		


	pars_NK_G = select(post,
		:Intercept_survG => :"??_surv",
		:b_z_survG => :"??z_surv",
		:Intercept_growG => :"??_grow",
		:b_z_growG => :"??z_grow",
		:sigma_growG => :??_grow,
		:Intercept_recrG => :"??_fec",
		:b_z_recrG => :"??z_fec"
	)

	pars_NK_G.??_surv = pars_NK_G.??_surv .+ post.b_NK_survG
	pars_NK_G.??z_surv = pars_NK_G.??z_surv .+ post.b_zNK_survG 
	pars_NK_G.??_grow = pars_NK_G.??_grow .+ post.b_NK_growG
	pars_NK_G.??z_grow = pars_NK_G.??z_grow .+ post.b_zNK_growG 
	pars_NK_G.??_fec = pars_NK_G.??_fec .+ post.b_NK_recrG
	pars_NK_G.??z_fec = pars_NK_G.??z_fec .+ post.b_zNK_recrG 


	function g_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		??= df.??_grow[row]
		??= df.??z_grow[row]
		
		??= df.??_grow[row]


		p_den_grow = zeros(size(z)[1],size(z)[1])
		?? = ?? .+ ?? * (z .- size_cen ) 
		for i in 1:nBigMatrix
			p_den_grow[:,i] = pdf.(Normal(??[i], ??), z1).*h
		end
		return(p_den_grow)
	end


	#@time g = g_z1z(pars_KG_G, z1, z, size_cen, row, 0.0, 0.0)
	# columns should sum to 1
	#sum.(eachcol(g))


	## Surival function
	#row = 1
	function s_z(df::AbstractDataFrame, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		??= df.??_surv[row]
		??= df.??z_surv[row]
		
		linear_p = ?? .+ ?? * (z .- size_cen)  # linear predictor
		p = 1 ./(1 .+ exp.(-linear_p))
		p = diagm(p)
		return(p)
	end

	s_z(pars_KG_G, z, size_cen, 1)

	## Reproduction function, logistic regression
	# function pb_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
	#   ??= df.??_rep[row]
	#   ??= df.??z_rep[row]
	#   linear_p = ?? .+ ?? * (z .- size_cen)       # linear predictor
	#   p = 1 ./(1 .+ exp.(-linear_p))
	#   p = diagm(p)
	#   return(p)
	# end



	## Recruitment function (N.B - from birth in spring to first summer), logistic regression
	function pr_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
		??= df.??_fec[row]
		??= df.??z_fec[row]
		

		linear_p = ?? .+ ?? * (z .- size_cen)     # linear predictor
		p = exp.(linear_p)*(1/2)
		p = diagm(p)
		return(p)
	end


	## Recruit size function
	function c_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
		size_cen::AbstractFloat, row::Integer)
		#??= df.rcz_int[row]
		#??= df.rcz_z[row]
		??= 0.4 #pars.rcz_sd[row]
		#   ??a= df.??_area_rcz[row]
		#   ??c= df.??_canopy_rcz[row]

		p_den_grow = zeros(nBigMatrix,nBigMatrix)
		?? = 7 .+ 0*(z .- size_cen) 
		for i in 1:nBigMatrix
			p_df = pdf.(Normal(??[i], ??), z1)*h
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

	#mat = IPM_KG["K"]



	Gsurv_mat = zeros(size(post)[1], nBigMatrix)
	Ggrow_mat = zeros(size(post)[1], nBigMatrix)
	#Grep_mat = zeros(size(pars_KG_G)[1], nBigMatrix)
	Gfec_mat = zeros(size(post)[1], nBigMatrix)
	Grcz_mat = zeros(size(post)[1], nBigMatrix)

	## make DataFrame to store the results


	Gres_IPM = DataFrame(zeros(size(post)[1], 15), :auto)
	Gres_IPM = select(Gres_IPM, :x1 => "lam_KG", :x2 => "lam_NK", :x3 => "delta_lam",
						:x4 => "sum_lam_eff", :x5 => "grow_con", :x6 => "fec_con", 
						:x7 => "rcz_con", :x8 => "sur_con",
						:x9 => "sum_con")

    
	for row in 1:size(post)[1]
		# Make projection kernels
		IPM_KG = mk_K(pars_KG_G, z1, z, size_cen, row)
		IPM_NK = mk_K(pars_NK_G, z1, z, size_cen, row)
		# calculate the population growth rate (??)

		vv_KG = eigen(IPM_KG["K"])
		vv_NK = eigen(IPM_NK["K"])

		??_KG = real(vv_KG.values[end])
		??_NK = real(vv_NK.values[end])
		Gres_IPM.lam_KG[row] = ??_KG
		Gres_IPM.lam_NK[row] = ??_NK

		## calculate average matrix
		K_avg = (IPM_NK["K"] + IPM_KG["K"])./2
		vv_avg = eigen(K_avg)

		# Normalize stable size distribution
		??_KG = real(vv_KG.vectors[:, end]) 
		??_NK = real(vv_NK.vectors[:, end]) 
		??_avg = real(vv_avg.vectors[:, end]) 

		# Reproductive value
		a_KG = eigen(transpose(IPM_KG["K"]))
		a_NK = eigen(transpose(IPM_NK["K"]))
		a_avg = eigen(transpose(K_avg))

		v_KG = real(a_KG.vectors[:, end]) 
		v_NK = real(a_NK.vectors[:, end])
		v_avg = real(a_avg.vectors[:, end]) / sum(real(a_avg.vectors[:, end]))
		v_avg = v_avg / dot(transpose(v_avg), ??_avg)

		## Sensitivity matrix

		sens_avg = v_avg * ??_avg' # this equivalent to outer() in R
		??K = IPM_NK["K"] - IPM_KG["K"]
		??_eff = ??K .* sens_avg
		???? = sum(??_eff)
		

		Gres_IPM.sum_lam_eff[row] = ????
		Gres_IPM.delta_lam[row] = ??_NK - ??_KG

		## Make life-response table

		one_mat = ones(nBigMatrix, nBigMatrix)

		# Function differences
		??_grow = g_z1z(pars_NK_G, z1, z, size_cen, row) - g_z1z(pars_KG_G, z1, z, size_cen, row)
		#??_rep = one_mat*(pb_z(pars_NK_G, z, size_cen, row) - pb_z(pars_KG_G, z, size_cen, row))
		??_fec = one_mat*(pr_z(pars_NK_G, z, size_cen, row) - pr_z(pars_KG_G, z, size_cen, row))
		??_rcz = (c_z1z(pars_NK_G, z1, z, size_cen, row) - c_z1z(pars_KG_G, z1, z, size_cen, row))
		??_sur = one_mat*(s_z(pars_NK_G, z, size_cen, row) - s_z(pars_KG_G, z, size_cen, row))

		# Function averages
		grow_avg = (g_z1z(pars_NK_G, z1, z, size_cen, row) + g_z1z(pars_KG_G, z1, z, size_cen, row))/2
		#rep_avg = (one_mat*(pb_z(pars_NK_G, z, size_cen, row) + pb_z(pars_KG_G, z, size_cen, row)))/2
		fec_avg = (one_mat*(pr_z(pars_NK_G, z, size_cen, row) + pr_z(pars_KG_G, z, size_cen, row)))/2
		rcz_avg = ((c_z1z(pars_NK_G, z1, z, size_cen, row) + c_z1z(pars_KG_G, z1, z, size_cen, row)))/2
		sur_avg = (one_mat*(s_z(pars_NK_G, z, size_cen, row) + s_z(pars_KG_G, z, size_cen, row)))/2

		# derivates

		# ????_grow= sur_avg
		# ????_rep= fec_avg*rcz_avg*sur_avg
		# ????_fec= rep_avg * fec_avg * sur_avg
		# ????_sur = grow_avg +  rep_avg * fec_avg * rcz_avg
		# ????_rcz = rep_avg * fec_avg * sur_avg

		# ??_grow = ??_grow .* sens_avg .* ????_grow
		# ??_rep = ??_rep .* sens_avg .* ????_rep
		# ??_fec = ??_fec .* sens_avg .* ????_fec
		# ??_rcz = ??_rcz .* sens_avg .* ????_rcz
		# ??_sur = ??_sur .* sens_avg .* ????_sur

		????_grow= sur_avg
		????_fec=  rcz_avg .* sur_avg
		????_sur = grow_avg .+   fec_avg .* rcz_avg
		????_rcz = fec_avg .* sur_avg

		??_grow = ??_grow .* sens_avg .* ????_grow
		??_fec = ??_fec .* sens_avg .* ????_fec
		??_rcz = ??_rcz .* sens_avg .* ????_rcz
		??_sur = ??_sur .* sens_avg .* ????_sur



		# to put in a DataFrame
		sur_con = sum(??_sur)
		grow_con = sum(??_grow)
		fec_con = sum(??_fec)
		rcz_con = sum(??_rcz)
		sum_con = sur_con + grow_con + fec_con + rcz_con

		Gres_IPM.sum_con[row] = sum_con
		Gres_IPM.sur_con[row] = sur_con
		Gres_IPM.grow_con[row] = grow_con
		Gres_IPM.fec_con[row] = fec_con
		Gres_IPM.rcz_con[row] = rcz_con

		Gsurv_mat[row, : ] =  sum.(eachcol(??_sur)) 
		Ggrow_mat[row, : ] =  sum.(eachcol(??_grow)) 
		Gfec_mat[row, : ] =  sum.(eachcol(??_fec)) 
		Grcz_mat[row, : ] =  sum.(eachcol(??_rcz)) 

	end

	return [Gres_IPM, Gsurv_mat, Ggrow_mat, Gfec_mat, Grcz_mat]

end


# these dataframes are from previous analyses in stan
function Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, 
    size_cen = 18.0)
      
	# size_cen = 18.0
    # nBigMatrix = 100
    # min_size = 2.0
    # max_size = 120

	U= max_size
    L= min_size
    m = nBigMatrix
    h = (U - L)/m
  
    z1 =  L .+ (collect(1:m) .- 0.5) .* h
    z1 = round.(z1, digits = 6)
    z = z1
    meshpts = z1
  
  
    # Get the parameters
  

    pars_GR_K = select(post,
		:Intercept_survK => :"??_surv",
        :b_z_survK => :"??z_surv",
        :Intercept_growK => :"??_grow",
        :b_z_growK => :"??z_grow",
		:b_z2_growK => :"??z2_grow",
        :sigma_growK => ByRow(x-> sqrt(x)) =>:??_grow,
        :Intercept_recrK => :"??_fec",
        :b_z_recrK => :"??z_fec"
  	)
  
  
  
    pars_NG_K = select(post,
        :Intercept_survK => :"??_surv",
        :b_z_survK => :"??z_surv",
        :Intercept_growK => :"??_grow",
        :b_z_growK => :"??z_grow",
		:b_z2_growK => :"??z2_grow",
       	:sigma_growK => ByRow(x-> sqrt(x)) =>:??_grow,
        :Intercept_recrK => :"??_fec",
        :b_z_recrK => :"??z_fec"
	)
  
  
  
    pars_NG_K.??_surv = pars_NG_K.??_surv .+ post.b_NG_survK
    pars_NG_K.??z_surv = pars_NG_K.??z_surv .+ post.b_zNG_survK 
    pars_NG_K.??_grow = pars_NG_K.??_grow .+ post.b_NG_growK
    pars_NG_K.??z_grow = pars_NG_K.??z_grow .+ post.b_zNG_growK 
    pars_NG_K.??_fec = pars_NG_K.??_fec .+ post.b_NG_recrK
    pars_NG_K.??z_fec = pars_NG_K.??z_fec .+ post.b_zNG_recrK 
  
  
  
    function g_z1zK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
      z2 = z .* z .- size_cen*size_cen  
	  zc = z .- size_cen 
		
	  ??= df.??_grow[row]
      ??= df.??z_grow[row]
	  ??2= df.??z2_grow[row]
     
      ??= df.??_grow[row]
      p_den_grow = zeros(size(z)[1],size(z)[1])
      ?? = (((?? .+ ?? .* zc .+ ??2 .* z2) .- z)./2) .+ z # average growth in two weeks
      #?? = round.(??, digits = 10)
      
      for i in 1:nBigMatrix
          
        p_den_grow[:,i] =   (pdf.(Normal(??[i], ??), z1).*h)
        
      end
      matex = zeros(nBigMatrix+1, nBigMatrix+1)
      matex[2:end,2:end] = p_den_grow
      return(matex)
    end
  
  
    # @time g = g_z1zK(pars_GR_K, z1, z, size_cen, 1)
  
  
	# sum.(eachcol(g)) 

    ## Surival function
  
  
    function s_zK(df::AbstractDataFrame, z::AbstractVector, 
      size_cen::AbstractFloat, row::Integer)
      ??= df.??_surv[row]
      ??= df.??z_surv[row]
     
      linear_p = ?? .+ ?? * (z .- size_cen)      # linear predictor
      p = 1 ./(1 .+ exp.(-linear_p))
      p = diagm(sqrt.(p))
      matex = zeros(nBigMatrix+1, nBigMatrix+1)
      matex[2:end,2:end] = p
      return(matex)
    end
  
  
  
    ## Reproduction function, logistic regression
    # function pb_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
    #   ??= df.??_rep[row]
    #   ??= df.??z_rep[row]
    #   linear_p = ?? .+ ?? * (z .- size_cen)       # linear predictor
    #   p = 1 ./(1 .+ exp.(-linear_p))
    #   p = diagm(p)
    #   return(p)
    # end
  
  
  
    ## Recruitment function (N.B - from birth in spring to first summer), logistic regression
    function pr_zK(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
      ??= df.??_fec[row]
      ??= df.??z_fec[row]
      linear_p = ?? .+ ?? * (z .- size_cen)   # linear predictor  # linear predictor
      p = exp.(linear_p).*(1/2)
      #p = diagm(p)
      matex = zeros(nBigMatrix+1, nBigMatrix+1)
      matex[1,2:end] = p
      return(matex)
  
    end
  
    ## Recruit size function
    function c_z1zK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
      size_cen::AbstractFloat, row::Integer)
      #??= df.rcz_int[row]
      #??= df.rcz_z[row]
      #??a= df.??_area_rcz[row]
      #??c= df.??_canopy_rcz[row]
      ??= 0.6 #pars.rcz_sd[row]
      p_den_grow = zeros(nBigMatrix+1,nBigMatrix+1)
  
      ?? = 4.4 .+ 0#*(z .- size_cen)
        #   for i in 1:nBigMatrix
        #     p_df = pdf.(Normal(??[i], ??), z1)*h
        #     for j in 1:nBigMatrix
        #       p_den_grow[j,i] = p_df[j]
        #     end
        #   end
      p_den_grow[2:end, 1] =  pdf.(Normal(??, ??), z).*h
      return(p_den_grow)
    end
   
  
    ##----------------------------------------------------
    ## Functions to build IPM kernels P, F, and K
  
  
  
    function mk_KK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, size_cen::AbstractFloat, 
      row::Integer, V::AbstractFloat)
        F = pr_zK(df, z, size_cen, row) * s_zK(df, z, size_cen, row)
        P = g_z1zK(df, z, z, size_cen, row) * s_zK(df, z, size_cen, row)
        A = (P + F) + (V .* c_z1zK(df, z1, z, size_cen, row))
        K = A * A
        out = Dict("K"=> K, "meshpts" => z, "P" => P, "F"=> F, "A" => A)
        return(out)
    end
  
    # Matrices for vital rates
    Ksurv_mat = zeros(size(post)[1], nBigMatrix+1)
    Kgrow_mat = zeros(size(post)[1], nBigMatrix+1)
    Krep_mat = zeros(size(post)[1], nBigMatrix+1)
    Kfec_mat = zeros(size(post)[1], nBigMatrix+1)
    Krcz_mat = zeros(size(post)[1], nBigMatrix+1)
  
    ## make DataFrame to store the results
  
  
    Kres_IPM = DataFrame(zeros(size(post)[1], 15), :auto)
  
  
    Kres_IPM = select(Kres_IPM, :x1 => "lam_GR", :x2 => "lam_NG", :x3 => "delta_lam",
              :x4 => "sum_lam_eff", :x5 => "grow_con", :x6 => "fec_con", 
              :x7 => "rcz_con", :x8 => "sur_con",
              :x9 => "sum_con")
  
  
  
    for row in 1:size(pars_GR_K)[1]
        # Make projection kernels
      IPM_GR = mk_KK(pars_GR_K, z1, z, size_cen, row, 0.7)
      IPM_NG = mk_KK(pars_NG_K, z1, z, size_cen, row, 0.7)
        # calculate the population growth rate (??)
  
      vv_GR = eigen(IPM_GR["K"])
      vv_NG = eigen(IPM_NG["K"])
  
      ??_GR = real(vv_GR.values[end])
      ??_NG = real(vv_NG.values[end])
      Kres_IPM.lam_GR[row] = ??_GR
      Kres_IPM.lam_NG[row] = ??_NG
        ??_NG -??_GR
  
        ## calculate average matrix
      K_avg = (IPM_NG["K"] + IPM_GR["K"])./2
      
      vv_avg = eigen(K_avg)
  
        # Normalize stable size distribution
      ??_GR = real(vv_GR.vectors[:, end]) ./ sum(real(vv_GR.vectors[:, end]))
      ??_NG = real(vv_NG.vectors[:, end]) ./  sum(real(vv_NG.vectors[:, end]))
      ??_avg = real(vv_avg.vectors[:, end]) ./ sum(real(vv_avg.vectors[:, end]))
  
        # Reproductive value
  
      a_avg = eigen(K_avg')
      v_avg = real(a_avg.vectors[:, end]) 
      v_avg = v_avg ./ dot(transpose(v_avg), ??_avg)
  
  
        ## Sensitivity matrix
  
      sens_avg = v_avg * ??_avg'  # this equivalent to outer() in R
      ??K = IPM_NG["K"] .- IPM_GR["K"]
      ??_eff = ??K .* sens_avg
      ???? = sum(??_eff)
      
  
      Kres_IPM.sum_lam_eff[row] = ????
      Kres_IPM.delta_lam[row] = ??_NG - ??_GR
  
      ## Make life-response table
  
      #one_mat = ones(nBigMatrix+1, nBigMatrix+1)
  
        # Function differences
      ??_grow = g_z1zK(pars_NG_K, z1, z, size_cen, row) .- g_z1zK(pars_GR_K, z1, z, size_cen, row)
      #??_rep = one_mat*(pb_z(pars_NG_K, z, size_cen, row) - pb_z(pars_GR_K, z, size_cen, row))
      ??_fec = (pr_zK(pars_NG_K, z, size_cen, row) .- pr_zK(pars_GR_K, z, size_cen, row))
      ??_rcz = (c_z1zK(pars_NG_K, z1, z, size_cen, row) .- c_z1zK(pars_GR_K, z1, z, size_cen, row))
      ??_sur = (s_zK(pars_NG_K, z, size_cen, row) .- s_zK(pars_GR_K, z, size_cen, row))
  
        # Function averages
      grow_avg = (g_z1zK(pars_NG_K, z1, z, size_cen, row) + g_z1zK(pars_GR_K, z1, z, size_cen, row))/2
      #rep_avg = (one_mat*(pb_z(pars_NG_K, z, size_cen, row) + pb_z(pars_GR_K, z, size_cen, row)))/2
      fec_avg = ((pr_zK(pars_NG_K, z, size_cen, row) + pr_zK(pars_GR_K, z, size_cen, row)))/2
      rcz_avg = ((c_z1zK(pars_NG_K, z1, z, size_cen, row) + c_z1zK(pars_GR_K, z1, z, size_cen, row)))/2
      sur_avg = ((s_zK(pars_NG_K, z, size_cen, row) + s_zK(pars_GR_K, z, size_cen, row)))/2
  
        # derivates
      I_mat = diagm(ones(nBigMatrix+1))
  
  
      ????_grow = kron(transpose(sur_avg), I_mat)
      ????_sur  = kron(transpose(I_mat),grow_avg) +  kron(transpose(I_mat),(fec_avg * rcz_avg))
      ????_fec  = kron(transpose(sur_avg), rcz_avg)
      ????_rcz  = kron(transpose(fec_avg * sur_avg), I_mat)
  
      ????A = kron(I_mat, K_avg) + kron(transpose(K_avg), I_mat)
  
      ??_grow = ??_grow .* reshape((([transpose(sens_avg[:]) ; ] * ????A )* ????_grow), (nBigMatrix+1,nBigMatrix+1))
      ??_fec = ??_fec .*  reshape((([transpose(sens_avg[:]) ; ] * ????A )* ????_fec), (nBigMatrix+1,nBigMatrix+1))
      ??_rcz = ??_rcz .*  reshape((([transpose(sens_avg[:]) ; ] * ????A )* ????_rcz), (nBigMatrix+1,nBigMatrix+1)) 
      ??_sur = ??_sur .* reshape((([transpose(sens_avg[:]) ; ] * ????A )* ????_sur), (nBigMatrix+1,nBigMatrix+1)) 
  
  
  
        # to put in a DataFrame
      sur_con = sum(??_sur)
      grow_con = sum(??_grow)
      #rep_con = sum(??_rep)
      fec_con = sum(??_fec)
      rcz_con = sum(??_rcz)
      sum_con = sur_con + grow_con + fec_con + rcz_con
      Kres_IPM.sum_con[row] = sum_con
  
  
      Kres_IPM.sur_con[row] = sur_con
      Kres_IPM.grow_con[row] = grow_con
      #Kres_IPM.rep_con[row] = rep_con
      Kres_IPM.fec_con[row] = fec_con
      Kres_IPM.rcz_con[row] = rcz_con
  
  
      Ksurv_mat[row, : ] =  sum.(eachcol(??_sur)) 
      Kgrow_mat[row, : ] =  sum.(eachcol(??_grow)) 
      #Krep_mat[row, : ] =  sum.(eachcol(??_rep)) 
      Kfec_mat[row, : ] =  sum.(eachcol(??_fec)) 
      Krcz_mat[row, : ] =  sum.(eachcol(??_rcz)) 
  
    end
    return [Kres_IPM, Ksurv_mat, Kgrow_mat, Kfec_mat, Krcz_mat]
  
end