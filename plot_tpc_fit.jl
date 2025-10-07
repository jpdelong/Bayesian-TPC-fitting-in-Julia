function plot_tpc_fit(temps,rint,chain_tpc,model)

# sample the priors
chain = sample(model, Prior(), 100)

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

p2 = density(chain[:ρ], label = "Prior", linewidth = 2)
density!(p2, chain_tpc[:ρ], label = "Posterior", color=:red, linewidth = 2, xlabel="ρ")

p3 = density(chain[:ΔT], label = "Prior", linewidth = 2)
density!(p3, chain_tpc[:ΔT], label = "Posterior", color=:red, linewidth = 2, xlabel="ΔT")

p4 = density(chain[:λ], label = "Prior", linewidth = 2)
density!(p4, chain_tpc[:λ], label = "Posterior", color=:red, linewidth = 2, xlabel="λ")

p5 = density(chain[:tmax], label = "Prior", linewidth = 2)
density!(p5, chain_tpc[:tmax], label = "Posterior", color=:red, linewidth = 2, xlabel="tmax")

l = @layout[a grid(2,2)]
p6 = plot(p1,p2,p3,p4,p5, layout = l, size=(600,400), dpi=600)

    return p6
end