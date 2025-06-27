using CairoMakie

export plot_energy_distribution, plot_angular_distribution, plot_momentum_distribution

"""
    plot_energy_distribution(ed::EnergyDistribution; title::String="", xlabel::String="Energy (a.u.)", 
                            ylabel::String="Probability", save_path::String="")

Plot the energy distribution of photoelectrons with proper labels.

# Arguments
- `ed::EnergyDistribution`: Energy distribution data to plot
- `title::String=""`: Plot title (auto-generated if empty)
- `xlabel::String="Energy (a.u.)"`: X-axis label
- `ylabel::String="Probability"`: Y-axis label  
- `save_path::String=""`: Path to save the plot (optional)

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_energy_distribution(ed::EnergyDistribution; title::String="", xlabel::String="Energy (a.u.)", 
                                 ylabel::String="Probability", save_path::String="")
    fig = Figure(size=(800, 600))
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
    axislegend(ax, position=:rt)
    
    # Save if path provided
    if !isempty(save_path)
        save(save_path, fig)
        println("Plot saved to: $save_path")
    end
    
    return fig
end

"""
    plot_angular_distribution(ad::AngularDistribution; title::String="", save_path::String="", 
                             plot_type::Symbol=:polar)

Plot the angular distribution of photoelectrons at fixed energy and theta.

# Arguments
- `ad::AngularDistribution`: Angular distribution data to plot
- `title::String=""`: Plot title (auto-generated if empty)
- `save_path::String=""`: Path to save the plot (optional)
- `plot_type::Symbol=:polar`: Plot type (:polar for polar plot, :cartesian for line plot)

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_angular_distribution(ad::AngularDistribution; title::String="", save_path::String="", 
                                  plot_type::Symbol=:polar)
    
    # Auto-generate title if not provided
    plot_title = if isempty(title)
        "Angular Distribution at E = $(round(ad.energy, digits=3)) a.u., θ = $(round(ad.θ * 180/π, digits=1))°"
    else
        title
    end
    
    if plot_type == :polar
        # Create polar plot
        fig = Figure(size=(800, 800))
        ax = PolarAxis(fig[1, 1], title=plot_title)
        
        # Convert to radial coordinates for polar plot
        # Since we have fixed θ and varying φ, we plot as a radial function
        r_values = ad.distribution ./ maximum(ad.distribution)  # Normalize for better visualization
        
        # Plot as polar line
        lines!(ax, ad.ϕ, r_values, linewidth=3, color=:blue)
        
    else  # cartesian plot
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], 
                 xlabel="Azimuthal Angle φ (radians)", 
                 ylabel="Probability Density",
                 title=plot_title)
        
        # Plot as line plot
        lines!(ax, ad.ϕ, ad.distribution, linewidth=2, color=:blue, label="P(φ)")
        
        # Add grid and formatting
        ax.xgridvisible = true
        ax.ygridvisible = true
        ax.xminorgridvisible = true
        ax.yminorgridvisible = true
        
        # Set axis limits with some padding
        if !isempty(ad.distribution) && maximum(ad.distribution) > 0
            ylims!(ax, 0, maximum(ad.distribution) * 1.1)
        end
        xlims!(ax, minimum(ad.ϕ), maximum(ad.ϕ))
        
        # Add x-axis tick labels in terms of π
        ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
        
        # Add legend
        axislegend(ax, position=:rt)
    end
    
    # Save if path provided
    if !isempty(save_path)
        save(save_path, fig)
        println("Plot saved to: $save_path")
    end
    
    return fig
end

"""
    plot_momentum_distribution(md::MomentumDistribution; title::String="", save_path::String="", 
                               slice_type::Symbol=:p_slice, slice_value::Float64=1.0, 
                               colormap=:viridis)

Plot the momentum distribution of photoelectrons with different visualization options.

# Arguments
- `md::MomentumDistribution`: Momentum distribution data to plot
- `title::String=""`: Plot title (auto-generated if empty)
- `save_path::String=""`: Path to save the plot (optional)
- `slice_type::Symbol=:p_slice`: Type of slice (:p_slice for fixed p, :theta_slice for fixed θ, :phi_slice for fixed φ)
- `slice_value::Float64=1.0`: Value at which to take the slice
- `colormap=:viridis`: Colormap for heatmaps

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_momentum_distribution(md::MomentumDistribution; title::String="", save_path::String="", 
                                   slice_type::Symbol=:theta_slice, slice_value::Float64=1.0,
                                   colormap=:viridis)
    
    fig = Figure(size=(800, 600))
    
    if slice_type == :p_slice
        # Plot at fixed momentum magnitude
        p_idx = argmin(abs.(md.p .- slice_value))
        actual_p = md.p[p_idx]
        
        # Create 2D slice: θ vs φ
        slice_data = md.distribution[p_idx, :, :]
        
        ax = Axis(fig[1, 1], 
                 xlabel="Azimuthal Angle φ (radians)", 
                 ylabel="Polar Angle θ (radians)",
                 title=isempty(title) ? "Momentum Distribution at p = $(round(actual_p, digits=3)) a.u." : title)
        
        hm = heatmap!(ax, md.φ, md.θ, slice_data', colormap=colormap)
        
        # Add colorbar
        Colorbar(fig[1, 2], hm, label="Probability Density")
        
        # Set axis ticks in terms of π
        ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
        ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        
    elseif slice_type == :theta_slice
        # Plot at fixed polar angle
        # θ_idx = argmin(abs.(md.θ .- slice_value))
        # actual_θ = md.θ[θ_idx]
        
        # Create 2D slice: p vs φ
        slice_data = [ md.distribution[p, 1, phi] for phi in eachindex(md.φ), p in eachindex(md.p)]
        
        ax = PolarAxis(fig[1, 1], 
                 xlabel="Azimuthal Angle φ (radians)", 
                 ylabel="Momentum p (a.u.)",
                 title=isempty(title) ? "Momentum Distribution at θ = $(round(actual_θ * 180/π, digits=1))°" : title)
        
        hm = surface!(ax, md.φ, md.p, zeros(size(dist)), color = slice_data, shading = NoShading, colormap = :coolwarm)
        
        # Add colorbar
        # Colorbar(fig[1, 2], hm, label="Probability Density")
        Colorbar(f[2, 1], p, vertical = false, flipaxis = false, label="Probability Density")
        
        # Set x-axis ticks in terms of π
        ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
        
    elseif slice_type == :phi_slice
        # Plot at fixed azimuthal angle
        φ_idx = argmin(abs.(md.φ .- slice_value))
        actual_φ = md.φ[φ_idx]
        
        # Create 2D slice: p vs θ
        slice_data = md.distribution[:, :, φ_idx]
        
        ax = Axis(fig[1, 1], 
                 xlabel="Polar Angle θ (radians)", 
                 ylabel="Momentum p (a.u.)",
                 title=isempty(title) ? "Momentum Distribution at φ = $(round(actual_φ * 180/π, digits=1))°" : title)
        
        hm = heatmap!(ax, md.θ, md.p, slice_data, colormap=colormap)
        
        # Add colorbar
        Colorbar(fig[1, 2], hm, label="Probability Density")
        
        # Set x-axis ticks in terms of π
        ax.xticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        
    else
        throw(ArgumentError("slice_type must be :p_slice, :theta_slice, or :phi_slice"))
    end
    
    # Add grid
    ax.xgridvisible = true
    ax.ygridvisible = true
    
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
