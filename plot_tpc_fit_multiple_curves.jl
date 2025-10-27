function plot_tpc_fit_multiple_curves(temps,rint,chain_tpc,model,treatments)

# get some information to automate how many sub parameters there are in the hieararchy
chain_headers = string.(names(chain_tpc))
indices_ρ = findall(x -> occursin("ρ", x), chain_headers)
indices_ΔT = findall(x -> occursin("ΔT", x), chain_headers)
indices_λ = findall(x -> occursin("λ", x), chain_headers)
indices_tmax = findall(x -> occursin("tmax", x), chain_headers)
numb_sub_params = length(indices_ρ)-2

# pick some colors
colors_to_use = cgrad(:berlin,numb_sub_params,categorical = true) # create a color gradient

# Plot the data
p1 = scatter(temps,rint,group=treatments,color_palette=[colors_to_use[1],colors_to_use[2],colors_to_use[3],colors_to_use[4],
    colors_to_use[5],colors_to_use[6],colors_to_use[7],],ylims=(-1, 4))


#p1 = scatter(temps,rint,group=treatments,color_palette=colors_to_use,ylims=(-1, 4))

# generate a new list of x values
xrange = LinRange(10, maximum(temps), 50)

# here's the model
@. lactin2(x, p) = exp.(p[1].*x) - exp.(p[1].*p[4] .- (p[4] - x) ./ p[2]) + p[3]

# sample the posteriors
chain_sample = sample(chain_tpc, 100; replace=false)

# convert the sample chain into a data frame
chn_tpc_sample = DataFrame(chain_sample)

# convert the whole chain into a data frame
chn_tpc_whole = DataFrame(chain_tpc)

# cycle through the number of sub parameters j
for j = 1:numb_sub_params
    # plot sample lines from posterior
    for i in 1:100 # Plot random samples
        # Get parameter values for the current sample
        ρ_samp = chn_tpc_sample[i,indices_ρ[j+2]+2] # for ρ
        ΔT_samp = chn_tpc_sample[i,indices_ΔT[j+2]+2] # for ΔT
        λ_samp = chn_tpc_sample[i,indices_λ[j+2]+2] # for λ
        tmax_samp = chn_tpc_sample[i,indices_tmax[j+2]+2] # for tmax

        y_fit = lactin2(xrange, [ρ_samp,ΔT_samp,λ_samp,tmax_samp])

        # Plot the fitted line with a transparent line style
        plot!(p1, xrange, y_fit, label=nothing, color=colors_to_use[j], linewidth=0.5, alpha=0.1)
    end

    # plot the curve for the median values
    fitted_ρ = median(chn_tpc_whole[:,indices_ρ[j+2]+2])
    fitted_ΔT = median(chn_tpc_whole[:,indices_ΔT[j+2]+2])
    fitted_λ = median(chn_tpc_whole[:,indices_λ[j+2]+2])
    fitted_tmax = median(chn_tpc_whole[:,indices_tmax[j+2]+2])

    y_fit = lactin2(xrange, [fitted_ρ,fitted_ΔT,fitted_λ,fitted_tmax])
    plot!(p1,xrange, y_fit, color=colors_to_use[j], label="", linewidth=3)

end

# sample the priors
chain = sample(model, Prior(), 100)

p2 = density(chain[:ρ], label = "Prior", color=:black, linewidth = 2, xlabel="ρ", xlims=(0.005, 0.18), legend=:best)
density!(p2, chain_tpc[:ρ], label = "Posterior", color=:gray, linewidth = 2)

p3 = density(chain[:ΔT], label = "Prior", color=:black, linewidth = 2, xlabel="ΔT", xlims=(4, 8), legend=false)
density!(p3, chain_tpc[:ΔT], label = "Posterior", color=:gray, linewidth = 2)

p4 = density(chain[:λ], label = "Prior", color=:black, linewidth = 2, xlabel="λ", xlims=(-2.5, 0), legend=false)
density!(p4, chain_tpc[:λ], label = "Posterior", color=:gray, linewidth = 2)

p5 = density(chain[:tmax], label = "Prior", color=:black, linewidth = 2, xlabel="tmax", xlims=(38, 41), legend=false)
density!(p5, chain_tpc[:tmax], label = "Posterior", color=:gray, linewidth = 2)

for j = 1:numb_sub_params
    density!(p2, chn_tpc_whole[:,indices_ρ[j+2]+2], label = "", color=colors_to_use[j], linewidth = 2)
    density!(p3, chn_tpc_whole[:,indices_ΔT[j+2]+2], label = "", color=colors_to_use[j], linewidth = 2)
    density!(p4, chn_tpc_whole[:,indices_λ[j+2]+2], label = "", color=colors_to_use[j], linewidth = 2)
    density!(p5, chn_tpc_whole[:,indices_tmax[j+2]+2], label = "", color=colors_to_use[j], linewidth = 2)
end

# put the plots together

l = @layout[a grid(2,2)]
p6 = plot(p1,p2,p3,p4,p5, layout = l, size=(1000,600), dpi=600)

    return p6
end