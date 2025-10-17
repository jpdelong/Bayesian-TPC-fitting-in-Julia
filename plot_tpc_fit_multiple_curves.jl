function plot_tpc_fit_multiple_curves(temps,rint,chain_tpc,model,scenario,lvls)
#=
chain_headers = names(chain_tpc)
findall(':ρ', chain_headers)
chain_parameters = chain_tpc.name_map.parameters
chain_tpc[:"ρ_rep[]"]

par1 = chain_parameters[5 + j]
chain_parameters[5 + lvls + 1 + j]
chain_parameters[5 + 2*lvls + 2 + j]
chain_parameters[5 + 3*lvls + 3 + j]

chain_tpc.string(chain_parameters[1])

chain_tpc[:, r'ρ']

julia> df[:, Between(begin, "x1")]

chain_tpc[:, Between(begin, "ρ")]

symbol("ρ_rep[",string(j),"]")
=#

# sample the priors
chain = sample(model, Prior(), 100)

# Plot the fitted curve
p1 = scatter(temps,rint,group=scenario,color_palette=[:blue,:red])

xrange = LinRange(10, maximum(temps), 50)
@. lactin2(x, p) = exp.(p[1].*x) - exp.(p[1].*p[4] .- (p[4] - x) ./ p[2]) + p[3]

################################
#################### first level
################################

# plot sample lines from posterior
for i in 1:100 # Plot 100 random samples for clarity
    # Get parameter values for the current sample
    ρ_samp = chain_tpc[:"ρ_rep[1]"][i]
    ΔT_samp = chain_tpc[:"ΔT_rep[1]"][i]
    λ_samp = chain_tpc[:"λ_rep[1]"][i]
    tmax_samp = chain_tpc[:"tmax_rep[1]"][i]

    y_fit = lactin2(xrange, [ρ_samp,ΔT_samp,λ_samp,tmax_samp])

    # Plot the fitted line with a transparent line style
    plot!(p1, xrange, y_fit, label=nothing, color=:blue, linewidth=0.5, alpha=0.2)
end

# plot the curve for the median values
fitted_ρ = median(chain_tpc[:"ρ_rep[1]"])
fitted_ΔT = median(chain_tpc[:"ΔT_rep[1]"])
fitted_λ = median(chain_tpc[:"λ_rep[1]"])
fitted_tmax = median(chain_tpc[:"tmax_rep[1]"])

y_fit = lactin2(xrange, [fitted_ρ,fitted_ΔT,fitted_λ,fitted_tmax])
plot!(p1,xrange, y_fit, color=:blue, label="Fit", linewidth=3)

################################
#################### second level
################################

# plot sample lines from posterior
for i in 1:100 # Plot 100 random samples for clarity
    # Get parameter values for the current sample
    ρ_samp = chain_tpc[:"ρ_rep[2]"][i]
    ΔT_samp = chain_tpc[:"ΔT_rep[2]"][i]
    λ_samp = chain_tpc[:"λ_rep[2]"][i]
    tmax_samp = chain_tpc[:"tmax_rep[2]"][i]

    y_fit = lactin2(xrange, [ρ_samp,ΔT_samp,λ_samp,tmax_samp])

    # Plot the fitted line with a transparent line style
    plot!(p1, xrange, y_fit, label=nothing, color=:red, linewidth=0.5, alpha=0.2)
end

# pull the median values
fitted_ρ = median(chain_tpc[:"ρ_rep[2]"])
fitted_ΔT = median(chain_tpc[:"ΔT_rep[2]"])
fitted_λ = median(chain_tpc[:"λ_rep[2]"])
fitted_tmax = median(chain_tpc[:"tmax_rep[2]"])
y_fit = lactin2(xrange, [fitted_ρ,fitted_ΔT,fitted_λ,fitted_tmax])
plot!(p1,xrange, y_fit, color=:red, label="Fit", linewidth=3)

p2 = density(chain[:ρ], label = "Prior", color=:black, linewidth = 2, xlabel="ρ", xlims=(0.005, 0.2), legend=:best)
density!(p2, chain_tpc[:ρ], label = "Posterior", color=:gray, linewidth = 2)
density!(p2, chain_tpc[:"ρ_rep[1]"], label = "Post_acute", color=:blue, linewidth = 2)
density!(p2, chain_tpc[:"ρ_rep[2]"], label = "Post_cg", color=:red, linewidth = 2)

p3 = density(chain[:ΔT], label = "Prior", color=:black, linewidth = 2, xlabel="ΔT", legend=false)
density!(p3, chain_tpc[:ΔT], label = "Posterior", color=:gray, linewidth = 2)
density!(p3, chain_tpc[:"ΔT_rep[1]"], label = "Post1", color=:blue, linewidth = 2)
density!(p3, chain_tpc[:"ΔT_rep[2]"], label = "Post2", color=:red, linewidth = 2)

p4 = density(chain[:λ], label = "Prior", color=:black, linewidth = 2, xlabel="λ", xlims=(-4, 4), legend=false)
density!(p4, chain_tpc[:λ], label = "Posterior", color=:gray, linewidth = 2)
density!(p4, chain_tpc[:"λ_rep[1]"], label = "Post1", color=:blue, linewidth = 2)
density!(p4, chain_tpc[:"λ_rep[2]"], label = "Post2", color=:red, linewidth = 2)

p5 = density(chain[:tmax], label = "Prior", color=:black, linewidth = 2, xlabel="tmax", xlims=(35, 45), legend=false)
density!(p5, chain_tpc[:tmax], label = "Posterior", color=:gray, linewidth = 2)
density!(p5, chain_tpc[:"tmax_rep[1]"], label = "Post1", color=:blue, linewidth = 2)
density!(p5, chain_tpc[:"tmax_rep[2]"], label = "Post2", color=:red, linewidth = 2)

#= Create a dummy plot for the legend, with the labels
p_legend = plot(
    [], [], label="Prior", linecolor=:black,
    [], [], label="Posterior", linecolor=:gray,
    [], [], label="Post1", linecolor=:blue,
    [], [], label="Post2", linecolor=:red,
    grid=false, showaxis=false, legend=:topright,
    framestyle=:none, background_color_subplot=nothing
)=#

l = @layout[a grid(2,2)]
p6 = plot(p1,p2,p3,p4,p5, layout = l, size=(1000,600), dpi=600)

    return p6
end