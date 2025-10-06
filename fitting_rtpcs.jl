dir = "C:/Users/jdelong2/OneDrive - University of Nebraska/Projects in progress/Paramecium caudatum TPC evolution/"
cd(dir)

#allfiles = readdir(dir)
using Turing, MCMCChains, StatsPlots, Distributions, CSV, DataFrames

tpcdata = CSV.read("P_caudatum_rTPC_spring_2020.csv",DataFrame)
temperature = tpcdata.Exp_Temp
rint = tpcdata.r
pop = tpcdata.Population
scenario = tpcdata.Scenario
    




scatter(temps,rint)

##########################################
###### non-hierarchical fitting
##########################################

# pull out data for a single tpc
indices = findall(tpcdata.Population .== 4 .&& tpcdata.Scenario .== "Acute" .&& tpcdata.Home_Temp .== 20)
temps = tpcdata.Exp_Temp[indices]
rint = tpcdata.r[indices]

# Lactin 2 thermal performance curve model
@model function tpc_lactin2(temperature,rint) 
    ρ ~ truncated(Normal(0.05,0.05),0,1)
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

# Plot the fitted curve
p1 = scatter(temps,rint)
xrange = LinRange(10, maximum(temps), 50)
@. lactin2(x, p) = exp.(p[1].*x) - exp.(p[1].*p[4] .- (p[4] - x) ./ p[2]) + p[3]

# plot sample lines from posterior
for i in 1:100 # Plot 100 random samples for clarity
    # Get parameter values for the current sample
    ρ_samp = chain_tpc[:ρ][i]
    ΔT_samp = chain_tpc[:ΔT][i]
    λ_samp = chain_tpc[:λ][i]
    tmax_samp = chain_tpc[:tmax][i]

    y_fit = lactin2(xrange, [ρ_samp,ΔT_samp,λ_samp,tmax_samp])

    # Plot the fitted line with a transparent line style
    plot!(p1, xrange, y_fit, label=nothing, color=:gray, alpha=0.2)
end

# plot the curve for the median values
# pull the media values
fitted_ρ = median(chain_tpc[:ρ])
fitted_ΔT = median(chain_tpc[:ΔT])
fitted_λ = median(chain_tpc[:λ])
fitted_tmax = median(chain_tpc[:tmax])

y_fit = lactin2(xrange, [fitted_ρ,fitted_ΔT,fitted_λ,fitted_tmax])
plot!(p1,xrange, y_fit, color=:red, label="Fit", linewidth=2)






##########################################
###### hierarchical fitting
##########################################


indices = findall(tpcdata.Population .== 4 .&& tpcdata.Home_Temp .== 20)
temps = tpcdata.Exp_Temp[indices]
rint = tpcdata.r[indices]
scenario = tpcdata.Scenario[indices]
levels = zeros(Int64,size(scenario))

# turn the string scenarios into integer levels for hierarchies
indices = findall(scenario .== "Acute")
levels[indices] .= 1
indices = findall(scenario .== "CG")
levels[indices] .= 2

# Lactin 2 thermal performance curve model - hierarchical
@model function tpc_lactin2_h(temperature,rint,levels) 
    n_rep = length(unique(levels)) # how many sublevel parameters, one for each scenario

    # model parameters, upper level in the hiararchy
    ρ ~ truncated(Normal(0.05,0.05),0,1)
    ΔT ~ truncated(Normal(2,10),0,10)
    λ ~ Normal(-2,3)
    tmax ~ truncated(Normal(40,5),minimum(temperature),1.5*maximum(temperature))

    # model parameters, lower level in the hierarchy
    σ_ρ ~ LogNormal(1/0.05)
    ρ_rep ~ filldist(truncated(Normal(ρ,σ_ρ),0,1), n_rep) # group-level parameters

    σ_ΔT ~ LogNormal(1/2)
    ΔT_rep ~ filldist(truncated(Normal(ΔT,σ_ΔT),0,10), n_rep) # group-level parameters

    σ_λ ~ LogNormal(1/1)
    λ_rep ~ filldist(Normal(λ,σ_λ), n_rep) # group-level parameters

    σ_tmax ~ LogNormal(1/40)
    tmax_rep ~ filldist(truncated(Normal(tmax,σ_tmax),minimum(temperature),1.5*maximum(temperature)), n_rep) # group-level parameters

    # replicate specific parameter set
    #p = [ρ_rep[j], ΔT_rep[j], λ_rep[j], tmax_rep[j]]

    for i in 1:length(rint)
        rint[i] ~ Normal(exp(ρ_rep[levels[i]]*temperature[i]) - exp(ρ_rep[levels[i]]*tmax_rep[levels[i]] - (tmax_rep[levels[i]] - temperature[i]) / ΔT_rep[levels[i]]) + λ_rep[levels[i]])
    end

    #=for j in 1:n_rep
        # replicate specific parameter set
        p = [ρ_rep[j], ΔT_rep[j], λ_rep[j], tmax_rep[j]]

        for i in 1:length(rint)
            rint[i,j] ~ Normal(exp(ρ*temperature[i]) - exp(ρ*tmax - (tmax - temperature[i]) / ΔT) + λ)
        end
    end=#
end

model = tpc_lactin2_h(temps,rint,levels) 

# call the fitting
chain_tpc = sample(
    model,
    NUTS(5000,0.65),
    MCMCSerial(),
    2000,
    #init_params = [(starting_a,starting_h,10)],
    1)

plot(chain_tpc)

# Plot the fitted curve
p1 = scatter(temps,rint,group=scenario)
xrange = LinRange(10, maximum(temps), 50)
@. lactin2(x, p) = exp.(p[1].*x) - exp.(p[1].*p[4] .- (p[4] - x) ./ p[2]) + p[3]

# grab parameters for first group
fitted_ρ = median(chain_tpc[:"ρ_rep[1]"])
fitted_ΔT = median(chain_tpc[:"ΔT_rep[1]"])
fitted_λ = median(chain_tpc[:"λ_rep[1]"])
fitted_tmax = median(chain_tpc[:"tmax_rep[1]"])
y_fit = lactin2(xrange, [fitted_ρ,fitted_ΔT,fitted_λ,fitted_tmax])
plot!(p1,xrange, y_fit, color=:blue, label="Fit", linewidth=2)


# grab parameters for second group
fitted_ρ = median(chain_tpc[:"ρ_rep[2]"])
fitted_ΔT = median(chain_tpc[:"ΔT_rep[2]"])
fitted_λ = median(chain_tpc[:"λ_rep[2]"])
fitted_tmax = median(chain_tpc[:"tmax_rep[2]"])
y_fit = lactin2(xrange, [fitted_ρ,fitted_ΔT,fitted_λ,fitted_tmax])
plot!(p1,xrange, y_fit, color=:red, label="Fit", linewidth=2)

savefig(p1,"hierarchical_tpc.png")

