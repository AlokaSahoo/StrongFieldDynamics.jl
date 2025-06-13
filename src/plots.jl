using CairoMakie

export plot_energy_distribution, plot_angular_distribution

"""
    plot_energy_distribution(ed::EnergyDistribution; title::String="", xlabel::String="Energy (a.u.)", 
                            ylabel::String="Differential Probability d²P/dΩdE", save_path::String="")

Plot the energy distribution of photoelectrons with proper labels.

# Arguments
- `ed::EnergyDistribution`: Energy distribution data to plot
- `title::String=""`: Plot title (auto-generated if empty)
- `xlabel::String="Energy (a.u.)"`: X-axis label
- `ylabel::String="Differential Probability d²P/dΩdE"`: Y-axis label  
- `save_path::String=""`: Path to save the plot (optional)

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_energy_distribution(ed::EnergyDistribution; title::String="", xlabel::String="Energy (a.u.)", 
                                 ylabel::String="Differential Probability d²P/dΩdE", save_path::String="")
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    
    # Auto-generate title if not provided
    plot_title = if isempty(title)
        "Energy Distribution at θ = $(round(ed.θ * 180/π, digits=1))°, φ = $(round(ed.ϕ * 180/π, digits=1))°"
    else
        title
    end
    ax.title = plot_title
    
    # Plot the energy spectrum
    lines!(ax, ed.energies, ed.spectrum, linewidth=2, color=:blue, label="Photoelectron Spectrum")
    
    # Add grid and formatting
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    
    # Set axis limits with some padding
    if !isempty(ed.spectrum) && maximum(ed.spectrum) > 0
        ylims!(ax, 0, maximum(ed.spectrum) * 1.1)
    end
    xlims!(ax, minimum(ed.energies), maximum(ed.energies))
    
    # Add legend if needed
    # axislegend(ax, position=:rt)
    
    # Save if path provided
    if !isempty(save_path)
        save(save_path, fig)
        println("Plot saved to: $save_path")
    end
    
    return fig
end

"""
    plot_angular_distribution(ad::AngularDistribution; title::String="", colormap=:viridis, save_path::String="")

Plot the angular distribution of photoelectrons as a polar heatmap.

# Arguments
- `ad::AngularDistribution`: Angular distribution data to plot
- `title::String=""`: Plot title (auto-generated if empty)
- `colormap=:viridis`: Colormap for the heatmap
- `save_path::String=""`: Path to save the plot (optional)

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_angular_distribution(ad::AngularDistribution; title::String="", colormap=:viridis, save_path::String="")
    fig = Figure(resolution=(800, 800))
    ax = PolarAxis(fig[1, 1], title = isempty(title) ? 
        "Angular Distribution at E = $(round(ad.energy, digits=3)) a.u." : title)
    
    # Create meshgrid for polar coordinates
    θ_grid = repeat(ad.θ', length(ad.ϕ), 1)
    ϕ_grid = repeat(ad.ϕ, 1, length(ad.θ))
    
    # Plot as polar heatmap
    surface!(ax, θ_grid, ϕ_grid, ad.distribution', colormap=colormap)
    
    # Add colorbar
    Colorbar(fig[1, 2], limits=extrema(ad.distribution), colormap=colormap, 
             label="Probability Density")
    
    # Save if path provided
    if !isempty(save_path)
        save(save_path, fig)
        println("Plot saved to: $save_path")
    end
    
    return fig
end


# f = Figure();

# ax = Axis(f[1,1],
#     title = "Title of the Plot",
#     xlabel = "Time (seconds)",
#     ylabel = "Value",
# );

# lines!(ax, a, b, linestyle = :dash, label= "label 1", color = :tomato)
# lines!(ax, a, b, linestyle=:dashdot, label= "label 2", color = :blue)

# axislegend(position = :rb)

# display(f)

# savefig("figurename.png", f)
