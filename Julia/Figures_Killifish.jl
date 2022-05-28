
# survival
KG = zeros(size(post_SurvK)[1], length(z));
NG = zeros(size(post_SurvK)[1], length(z));

size_cen = 18.0

z = collect(5:0.25:100);


for j in 1:size(post_SurvK)[1]
    KG[j,:] =  (p_linkK(post_SurvK, z, size_cen, 0, j))
    NG[j,:] =  (p_linkK(post_SurvK, z, size_cen, 1, j))
end

my_summary(KG)

# Plot NG
p = My_Logistic.(my_summary(KG)).*100
Fig_1A =plot(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), 
    linewidth = 5, label = "KG", title = "a)", titleloc = :left, legend = :bottomright    
)

p = My_Logistic.(my_summary(NG)).*100
plot!(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), linewidth = 5, label = "NG")
#xlabel!("Initial size (mm)")
ylabel!("Survival (%)")
ylims!((0,100))
xlims!(5,95)
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, label= :false)

    # estimation graphs

lab= round(LOS(post_SurvK.b_NG), digits = 2)


Fig_1A_α = plot(kde(post_SurvK.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "            Statistical test \n                KG > NG  \n $(lab)%", titleloc = :left, titlefontsize = 9)
plot!(kde(post_SurvK.b_Intercept .+ post_SurvK.b_NG), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

lab= round(LOS(post_SurvK."b_z.NG"), digits = 2)
Fig_1A_β = plot(kde(post_SurvK.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n  $(lab)%", titlefontsize = 9, titleloc= :left, ticks =:false)
plot!(kde(post_SurvK.b_z .+ post_SurvK."b_z.NG"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope")

l = @layout [
    a [b{0.4h}
       c{0.4h}]
]

FA = plot(Fig_1A, Fig_1A_α, Fig_1A_β,
    layout = l
)

# growth
KG = zeros(size(post_GrowthK)[1], length(z));
NG = zeros(size(post_GrowthK)[1], length(z));


function p_linkK_growth(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NG::Integer, row::Integer)
    zc = z .- size_cen
    α= df.b_Intercept[row]
    βz= df.b_z[row]
    βNG= df."b_NG"[row]
    βzNG= df."b_z.NG"[row]
    μ = α .+ βNG .* NG .+  (βz .+ βzNG .* NG) .* zc  # linear predictor
    return(μ)
end



for i in 1:size(post_GrowthK)[1]
    KG[i,:] =  (p_linkK(post_GrowthK, z, size_cen, 0, i))
    NG[i,:] =  (p_linkK(post_GrowthK, z, size_cen, 1, i))
end


# Plot NG
p = my_summary(KG)
Fig_2A =plot(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), 
    linewidth = 5, label = "KG", title = "b)", titleloc = :left, legend = :bottomright    
)

p = my_summary(NG)
plot!(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), linewidth = 5, label = "NG")
#xlabel!("Initial size (mm)")
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, legend = :false)

ylabel!("Final size (mm)")
ylims!((5,95))
xlims!(5,95)



scatter!(DataK.SL1_mm, DataK.SL2_mm, groups = DataK.NG, c= [:lightskyblue, :red], alpha = 0.8, label =false)

xlabel!("Killifish size (mm)")
    # estimation graphs



lab= round(LOS(post_GrowthK.b_NG), digits = 2)

Fig_2A_α = plot(kde(post_GrowthK.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "\n $(lab)%", titlefontsize = 9, titleloc = :left)
plot!(kde(post_GrowthK.b_Intercept .+ post_GrowthK.b_NG), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

lab= round(LOS(post_GrowthK."b_z.NG"), digits=2)
Fig_2A_β = plot(kde(post_GrowthK.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n $(lab)%", titlefontsize = 9, titleloc = :left, ticks =:false)
plot!(kde(post_GrowthK.b_z .+ post_GrowthK."b_z.NG"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope")


l = @layout [
    a [b{0.4h}
       c{0.4h}]
]


FB = plot(Fig_2A, Fig_2A_α, Fig_2A_β,
    layout = l
)



# growth
KG = zeros(size(post_ReprK)[1], length(z));
NG = zeros(size(post_ReprK)[1], length(z));



for i in 1:size(post_ReprK)[1]
    KG[i,:] =  (p_linkK(post_ReprK, z, size_cen, 0, i))
    NG[i,:] =  (p_linkK(post_ReprK, z, size_cen, 1, i))
end


# Plot NG
p = exp.(my_summary(KG))
Fig_3A =plot(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), 
    linewidth = 5, label = "KG", title = "c)", titleloc = :left, legend = :bottomright    
)

p = exp.(my_summary(NG))
plot!(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), linewidth = 5, label = "NG")
#xlabel!("Initial size (mm)")
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, legend = :false)

ylabel!("Offspring (N)")

scatter!(DataK.SL1_mm, DataK.Recr, groups = DataK.NG, c= [:lightskyblue, :red], alpha = 0.8, label =false)


ylims!((-1,20))

xlims!(5,90)


    # estimation graphs



lab= round(LOS(post_ReprK.b_NG), digits = 2)

Fig_3A_α = plot(kde(post_ReprK.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "\n $(lab)%", titlefontsize = 9, titleloc = :left)
plot!(kde(post_ReprK.b_Intercept .+ post_ReprK.b_NG), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

lab= round(LOS(post_GrowthK."b_z.NG"), digits=2)
Fig_3A_β = plot(kde(post_ReprK.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n $(lab)%", titlefontsize = 9, titleloc = :left, ticks =:false)
plot!(kde(post_ReprK.b_z .+ post_ReprK."b_z.NG"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope")


l = @layout [
    a [b{0.4h}
       c{0.4h}]
]


FC = plot(Fig_3A, Fig_3A_α, Fig_3A_β,
    layout = l
)


plot(FA, FB, layout = (2,1), size = (500, 700))



