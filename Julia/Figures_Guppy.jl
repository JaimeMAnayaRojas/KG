# size clases
size_cen = 18.0
z =  collect(5:0.1:30)
z1 = z
# survival
KG = zeros(size(post_SurvG)[1], length(z));
NK = zeros(size(post_SurvG)[1], length(z));


for i in 1:size(post_SurvG)[1]
    KG[i,:] =  (p_link(post_SurvG, z, size_cen, 0, i))
    NK[i,:] =  (p_link(post_SurvG, z, size_cen, 1, i))
end

my_summary(KG)

# Plot NG
p = My_Logistic.(my_summary(KG)).*100
Fig_1A =plot(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), 
    linewidth = 5, label = "KG", title = "a)", titleloc = :left, legend = :bottomright    
)

p = My_Logistic.(my_summary(NK)).*100
plot!(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), linewidth = 5, label = "NK")
#xlabel!("Initial size (mm)")
ylabel!("Survival (%)")
ylims!((0,100))
xlims!(8,30)
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, label= :false)

    # estimation graphs

lab= round(LOS(post_SurvG.b_NK), digits = 2)


Fig_1A_α = plot(kde(post_SurvG.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "            Statistical test \n                KG > NK  \n $(lab)%", titleloc = :left, titlefontsize = 9)
plot!(kde(post_SurvG.b_Intercept .+ post_SurvG.b_NK), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

lab= round(LOS(post_SurvG."b_z.NK"), digits = 2)
Fig_1A_β = plot(kde(post_SurvG.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n  $(lab)%", titlefontsize = 9, titleloc= :left, ticks =:false)
plot!(kde(post_SurvG.b_z .+ post_SurvG."b_z.NK"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope")


l = @layout [
    a [b{0.4h}
       c{0.4h}]
]

FA = plot(Fig_1A, Fig_1A_α, Fig_1A_β,
    layout = l
)

# growth
KG = zeros(size(post_GrowthG)[1], length(z));
NK = zeros(size(post_GrowthG)[1], length(z));

for i in 1:size(post_GrowthG)[1]
    KG[i,:] =  (p_link(post_GrowthG, z, size_cen, 0, i))
    NK[i,:] =  (p_link(post_GrowthG, z, size_cen, 1, i))
end


# Plot NG
p = my_summary(KG)
Fig_2A =plot(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), 
    linewidth = 5, label = "KG", title = "b)", titleloc = :left, legend = :bottomright    
)

p = my_summary(NK)
plot!(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), linewidth = 5, label = "NK")
#xlabel!("Initial size (mm)")


plot!([18,18], [minimum(p.l95),maximum(p.u95)], linestyle = :dash, colour = :gray, legend = :false)

ylabel!("Growth ln(z₁/z)")
#ylims!((5,30))
xlims!(8,30)


DataG.growth = log.(DataG.SL2_mm ./ DataG.SL1_mm)
scatter!(DataG.SL1_mm, DataG.growth, groups = DataG.NK, c= [:lightskyblue, :red], alpha = 0.8, label =false)

    # estimation graphs



lab= round(LOS(post_GrowthG.b_NK), digits = 2)

Fig_2A_α = plot(kde(post_GrowthG.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "\n $(lab)%", titlefontsize = 9, titleloc = :left)
plot!(kde(post_GrowthG.b_Intercept .+ post_GrowthG.b_NK), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

lab= round(LOS(post_GrowthG."b_z.NK"), digits=2)
Fig_2A_β = plot(kde(post_GrowthG.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n $(lab)%", titlefontsize = 9, titleloc = :left, ticks =:false)
plot!(kde(post_GrowthG.b_z .+ post_GrowthG."b_z.NK"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope")



l = @layout [
    a [b{0.4h}
       c{0.4h}]
]
FB = plot(Fig_2A, Fig_2A_α, Fig_2A_β,
    layout = l
)



# growth
KG = zeros(size(post_ReprG)[1], length(z));
NK = zeros(size(post_ReprG)[1], length(z));



for i in 1:size(post_ReprG)[1]
    KG[i,:] =  (p_link(post_ReprG, z, size_cen, 0, i))
    NK[i,:] =  (p_link(post_ReprG, z, size_cen, 1, i))
end


# Plot NG
p = exp.(my_summary(KG))
Fig_3A =plot(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), 
    linewidth = 5, label = "KG", title = "c)", titleloc = :left, legend = :bottomright    
)

p = exp.(my_summary(NK))
plot!(z, p[:,:median], ribbon = (p.median .- p.l95, p.u95 .- p.median), linewidth = 5, label = "NK")
#xlabel!("Initial size (mm)")
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, legend = :false)

ylabel!("Offspring (N)")

scatter!(DataG.SL1_mm, DataG.Recr, groups = DataG.NK, c= [:lightskyblue, :red], alpha = 0.8, label =false)


ylims!((-1,20))

xlims!(5,30)

xlabel!("Guppy size (mm)")

    # estimation graphs



lab= round(LOS(post_ReprG.b_NK), digits = 2)

Fig_3A_α = plot(kde(post_ReprG.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "\n $(lab)%", titlefontsize = 9, titleloc = :left)
plot!(kde(post_ReprG.b_Intercept .+ post_ReprG.b_NK), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

lab= round(LOS(post_GrowthG."b_z.NK"), digits=2)
Fig_3A_β = plot(kde(post_ReprG.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n $(lab)%", titlefontsize = 9, titleloc = :left, ticks =:false)
plot!(kde(post_ReprG.b_z .+ post_ReprG."b_z.NK"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope")


l = @layout [
    a [b{0.4h}
       c{0.4h}]
]

FC = plot(Fig_3A, Fig_3A_α, Fig_3A_β,
    layout = l
)


plot(FA, FB, FC, layout = (3,1), size = (500, 700))


