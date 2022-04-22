using Distributions
using CSV
using DataFrames
using LinearAlgebra
using StatsPlots
using JLD2
using RCall
using Plots
using Plots.Measures
using CategoricalArrays

Gz = CSV.read("data/01_GuppyField_data.csv", DataFrame)
Kz = CSV.read( "data/01_KillifishField_data.csv", DataFrame)


println(describe(Gz))
replace!(Gz.SL1_mm, missing=>NaN)
replace!(Gz.SL2_mm, missing=>NaN)



replace!(Kz.SL1_mm, missing=>NaN)
replace!(Kz.SL2_mm, missing=>NaN)




filter!(:Location => x -> x !="JOE", Gz)
filter!(:Location => x -> x !="JOE", Kz)

Gz.Location = categorical(Gz.Location)



# Caigual
p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "CAI", Gz)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1a = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "t₀ (105)")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "t₁ (109)", 
        titlefont = font(10),  
        title = "a) Caigual-KG", titleloc = :left
)
ylabel!("Frequency (N)")


p = filter([:NK, :Location] => (x,y) -> x == 1 && y == "CAI", Gz)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))


p1b = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "103")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "170", 
        titlefont = font(10),  
        title = "b) Caigual-NK", titleloc = :left
)

# Naranjo (NAR)


p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "NAR", Gz)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1c = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "90")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "93", 
        titlefont = font(10),  
        title = "c) Naranjo-KG", titleloc = :left
)
ylabel!("Frequency (N)")


p = filter([:NK, :Location] => (x,y) -> x == 1 && y == "NAR", Gz)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1d = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "34")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "84", 
        titlefont = font(10),  
        title = "d) Naranjo-NK", titleloc = :left
)



# Quare 1 (QUA)

p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "QUA", Gz)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))


p1e = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "60")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "95", 
        titlefont = font(10),  
        title = "e) Quare 1-KG", titleloc = :left
)
ylabel!("Frequency (N)")


p = filter([:NK, :Location] => (x,y) -> x == 1 && y == "QUA", Gz)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1f = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "83")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "132", 
        titlefont = font(10),  
        title = "f) Quare 1-NK", titleloc = :left
)



# Quare 2 (QUA2)


p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "QUA2", Gz)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1g = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "16")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "10", 
        titlefont = font(10),  
        title = "g) Quare 2-KG", titleloc = :left
)
ylabel!("Frequency (N)")

xlabel!("Size (mm")

p = filter([:NK, :Location] => (x,y) -> x == 1 && y == "QUA2", Gz)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))


p1h = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "39")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "28", 
        titlefont = font(10),  
        title = "f) Quare 1-NK", titleloc = :left
)
xlabel!("Size (mm")


plot(p1a, p1b, p1c, p1d,p1e,p1f,p1g,p1h, bins= 8, layout = (4,2), size = (600, 700))
ylims!((0,25))
xlims!((5,32))
savefig("plots/Figure_S3.png")


plot(p1a, p1b, p1c, p1d,p1e,p1f,p1g,p1h, bins= 8, layout = (4,2), size = (600, 700))
#ylims!((0,25))
xlims!((5,32))
savefig("plots/Figure_S3b.png")




# Caigual
p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "CAI", Kz)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1a = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "t₀ (12)")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "t₁ (10)", 
        titlefont = font(10),  
        title = "a) Caigual-KG", titleloc = :left
)
ylabel!("Frequency (N)")


p = filter([:NG, :Location] => (x,y) -> x == 1 && y == "CAI", Kz)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))


p1b = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "17")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "18", 
        titlefont = font(10),  
        title = "b) Caigual-NG", titleloc = :left
)

# Naranjo (NAR)


unique(Kz.Location)

p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "NAR", Kz)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1c = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "22")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "26", 
        titlefont = font(10),  
        title = "c) Naranjo-KG", titleloc = :left
)
ylabel!("Frequency (N)")


p = filter([:NG, :Location] => (x,y) -> x == 1 && y == "NAR", Kz)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1d = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "28")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "36", 
        titlefont = font(10),  
        title = "d) Naranjo-NG", titleloc = :left
)



# Quare 1 (QUA)


p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "QUA", Kz)


length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))


p1e = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "108")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "174", 
        titlefont = font(10),  
        title = "e) Quare 1-KG", titleloc = :left
)
ylabel!("Frequency (N)")

# WHY ARE THERE NOT KILLIFISH IN THIS POOL?
p = filter([:NG, :Location] => (x,y) -> x == 1 && y == "QUA", Kz)


length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1f = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "66")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "108", 
        titlefont = font(10),  
        title = "f) Quare 1-NG", titleloc = :left
)



# Quare 2 (QUA2)


p = filter([:KG, :Location] => (x,y) -> x == 1 && y == "QUA2", Kz)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))

p1g = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "49")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "30", 
        titlefont = font(10),  
        title = "g) Quare 2-KG", titleloc = :left
)
ylabel!("Frequency (N)")

xlabel!("Size (mm")

p = filter([:NG, :Location] => (x,y) -> x == 1 && y == "QUA2", Kz)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))


p1h = histogram(p.SL1_mm, bins = 30, alpha=0.4, color = :black, label = "25")

histogram!(p.SL2_mm, bins = 30, alpha=0.5, color = :orange, label = "17", 
        titlefont = font(10),  
        title = "f) Quare 1-NG", titleloc = :left
)
xlabel!("Size (mm")


plot(p1a, p1b, p1c, p1d,p1e,p1f,p1g,p1h, bins= 12, layout = (4,2), size = (600, 700))
ylims!((0,40))
xlims!((5,100))
savefig("plots/Figure_S4.png")

plot(p1a, p1b, p1c, p1d,p1e,p1f,p1g,p1h, bins= 12, layout = (4,2), size = (600, 700))
#ylims!((0,40))
xlims!((5,100))
savefig("plots/Figure_S4b.png")