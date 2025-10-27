# this script pulls data and runs Bayesian fits of TPC models to little r data
# it uses statsplots for the Bayes fitting figures like plotting chains 
# it calls convenience functions to plot the fits along with priors/posteriors

using Turing, MCMCChains, StatsPlots, Distributions, CSV, DataFrames

# first try some single curves to test

# test data 1
dir = "C:/Users/jdelong2/OneDrive - University of Nebraska/Projects in progress/Paramecium caudatum TPC evolution/"
cd(dir)

# if you can find your file, try this and see if it's there
#allfiles = readdir(dir)

tpcdata = CSV.read("P_caudatum_rTPC_spring_2020.csv",DataFrame)
    # pull out data for a single tpc
    indices = findall(tpcdata.Population .== 4 .&& tpcdata.Scenario .== "Acute" .&& tpcdata.Home_Temp .== 20)
    temps = tpcdata.Exp_Temp[indices]
    rint = tpcdata.r[indices]

# test data 2
dir = "C:/Users/jdelong2/OneDrive - University of Nebraska/Projects in progress/Julia based fitting/"
cd(dir)
#tpcdata = CSV.read("TPC_for_Paramecium_aurelia_bypredation.csv",DataFrame)

tpcdata = CSV.read("TPC_paramecium_caudatum_evolved_lines_2021.csv",DataFrame)


    # pull out data for a single tpc
    indices = findall(tpcdata.Treatment .== "Predator")
    temps = 0.0 .+ tpcdata.Temp[indices]
    rint = tpcdata.r[indices]

##########################################
###### non-hierarchical, single-curve fitting
##########################################

# Lactin 2 thermal performance curve model
@model function tpc_lactin2(temperature,rint) 
    ρ ~ truncated(Normal(0.05,0.1),0,1)
    ΔT ~ truncated(Normal(2,10),0,10)
    λ ~ Normal(-2,3)
    tmax ~ truncated(Normal(40,5),minimum(temperature),1.5*maximum(temperature))

	for i in 1:length(rint)
        rint[i] ~ Normal(exp(ρ*temperature[i]) - exp(ρ*tmax - (tmax - temperature[i]) / ΔT) + λ)
    end
end

model = tpc_lactin2(temps,rint)

# call the fitting
chain_tpc = sample(
    model,
    NUTS(5000,0.65),
    MCMCSerial(),
    2000,
    #init_params = [(starting_a,starting_h,10)],
    1)

plot(chain_tpc)
summarystats(chain_tpc)

# call the custom fitting function
include("plot_tpc_fit_single_curve.jl")
fig1 = plot_tpc_fit(temps,rint,chain_tpc,model)
savefig(fig1,"single_fit_tpc.png")


##########################################################
###### hierarchical fitting for more than one curve
##########################################################

# test data 1
indices = findall(tpcdata.Population .== 4 .&& tpcdata.Home_Temp .== 20)
temps = tpcdata.Exp_Temp[indices]
rint = tpcdata.r[indices]
treatments = tpcdata.Scenario[indices]
levels = zeros(Int64,size(treatments))

# turn the string scenarios into integer levels for hierarchies
indices = findall(treatments .== "Acute")
levels[indices] .= 1
indices = findall(treatments .== "CG")
levels[indices] .= 2

# test data 2 - pull out two levels
indices2 = findall(tpcdata.Treatment .== "Start") # identify the set we are not using
tpcdata = tpcdata[Not(indices2),:] # drop the starts from the dataset
indices = findall(tpcdata.Treatment .== "Predator" .|| tpcdata.Treatment .== "No predator")
treatments = tpcdata.Treatment[indices]
temps = 0.0 .+ tpcdata.Temp[indices]
rint = tpcdata.r[indices]
levels = zeros(Int64,size(treatments))

# test data 2 - use all three levels
treatments = tpcdata.Treatment
temps = 0.0 .+ tpcdata.Temp
rint = tpcdata.r

levels = zeros(Int64,size(tpcdata,1))
# turn the string scenarios into integer levels for hierarchies
indices = findall(tpcdata.Treatment .== "Predator")
levels[indices] .= 1
indices = findall(tpcdata.Treatment .== "No predator")
levels[indices] .= 2
indices = findall(tpcdata.Treatment .== "Start")
levels[indices] .= 3

# test data 3 - use all three levels
treatments = tpcdata.PopID
temps = 0.0 .+ tpcdata.Temp
rint = tpcdata.r2
levels = tpcdata.Replicate

# Lactin 2 thermal performance curve model - hierarchical
@model function tpc_lactin2_h(temperature,rint,levels) 
    n_rep = length(unique(levels)) # how many sublevel parameters, one for each scenario

    # model parameters, upper level in the hieararchy
    ρ ~ truncated(Normal(0.05,0.05),0,1)
    ΔT ~ truncated(Normal(2,10),0,7)
    λ ~ Normal(-2,3)
    tmax ~ truncated(Normal(40,5),minimum(temperature),1.5*maximum(temperature))

    # model parameters, lower level in the hierarchy
    σ_ρ ~ LogNormal(1/0.05)
    ρ_rep ~ filldist(truncated(Normal(ρ,σ_ρ),0,1), n_rep) # group-level parameters

    σ_ΔT ~ LogNormal(1/2)
    ΔT_rep ~ filldist(truncated(Normal(ΔT,σ_ΔT),0,7), n_rep) # group-level parameters

    σ_λ ~ LogNormal(1/1)
    λ_rep ~ filldist(Normal(λ,σ_λ), n_rep) # group-level parameters

    σ_tmax ~ LogNormal(1/40)
    tmax_rep ~ filldist(truncated(Normal(tmax,σ_tmax),minimum(temperature),1.5*maximum(temperature)), n_rep) # group-level parameters

    for i in 1:length(rint)
        rint[i] ~ Normal(exp(ρ_rep[levels[i]]*temperature[i]) - exp(ρ_rep[levels[i]]*tmax_rep[levels[i]] - (tmax_rep[levels[i]] - temperature[i]) / ΔT_rep[levels[i]]) + λ_rep[levels[i]])
    end
end

model = tpc_lactin2_h(temps,rint,levels) 

# call the fitting
chain_tpc = sample(
    model,
    NUTS(20000,0.75),
    MCMCSerial(),
    10000,
    #init_params = [(starting_a,starting_h,10)],
    1)

plot(chain_tpc)

# call the custom fitting function
include("plot_tpc_fit_multiple_curves.jl")
fig2 = plot_tpc_fit_multiple_curves(temps,rint,chain_tpc,model,treatments)

savefig(fig2,"tpc_hot_warm_adapted.png")
