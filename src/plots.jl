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
function plot_energy_distribution(ed::StrongFieldDynamics.EnergyDistribution; title::String="", xlabel::String="Energy (a.u.)", 
                                 ylabel::String="Normalized Probability", save_path::String="")
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    
    # Auto-generate title if not provided
    plot_title = if isempty(title)
        # Convert intensity back to W/cm²
        I_Wcm2 = ed.pulse.I₀ * 3.51e16
        pulse_info = "I₀ = $(round(I_Wcm2/1e14, digits=1))×10¹⁴ W/cm², $(ed.pulse.np) cycles, ϵ = $(ed.pulse.ϵ), helicity = $(ed.pulse.helicity)"
        "Energy Distribution at θ = $(round(ed.θ * 180/π, digits=1))°, φ = $(round(ed.ϕ * 180/π, digits=1))°, $(pulse_info)"
    else
        title
    end
    ax.title = plot_title
    
    # Normalize the distribution for plotting
    normalized_distribution = if !isempty(ed.distribution) && maximum(ed.distribution) > 0
        ed.distribution ./ maximum(ed.distribution)
    else
        ed.distribution
    end
    
    # Plot the normalized energy distribution
    lines!(ax, ed.energies, normalized_distribution, linewidth=2, color=:blue, label="Photoelectron Distribution")
    
    # Add grid and formatting
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    
    # Set axis limits with some padding
    if !isempty(normalized_distribution) && maximum(normalized_distribution) > 0
        ylims!(ax, 0, maximum(normalized_distribution) * 1.1)
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
        # Convert intensity back to W/cm²
        I_Wcm2 = ad.pulse.I₀ * 3.51e16
        pulse_info = "I₀ = $(round(I_Wcm2/1e14, digits=1))×10¹⁴ W/cm², $(ad.pulse.np) cycles, ϵ = $(ad.pulse.ϵ), helicity = $(ad.pulse.helicity)"
        "Angular Distribution at E = $(round(ad.energy, digits=3)) a.u., θ = $(round(ad.θ * 180/π, digits=1))°, $(pulse_info)"
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
                               colormap=nothing)

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
                                   colormap=nothing)
    
    # Define custom colormap
    custom_colormap = cgrad([:white, :blue, :cyan, :green, :yellow, :orange, :red])
    colormap = isnothing(colormap) ? custom_colormap : colormap
    
    # Generate pulse info for title
    I_Wcm2 = md.pulse.I₀ * 3.51e16
    pulse_info = "I₀ = $(round(I_Wcm2/1e14, digits=1))×10¹⁴ W/cm², $(md.pulse.np) cycles, ϵ = $(md.pulse.ϵ), helicity = $(md.pulse.helicity)"
    
    if slice_type == :p_slice
        # Plot at fixed momentum magnitude
        p_idx = argmin(abs.(md.p .- slice_value))
        actual_p = md.p[p_idx]
        
        # Create 2D slice: θ vs φ (if multiple theta values exist)
        if length(md.θ) > 1
            slice_data = md.distribution[p_idx, :, :]
            # Normalize the data
            slice_data_norm = slice_data ./ maximum(slice_data)
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Momentum Distribution at p = $(round(actual_p, digits=3)) a.u., $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Azimuthal Angle φ (radians)", 
                     ylabel="Polar Angle θ (radians)",
                     title=plot_title)
            
            hm = heatmap!(ax, md.φ, md.θ, slice_data_norm', colormap=colormap)
            Colorbar(fig[1, 2], hm, label="Normalized Probability Density")
            
            # Set axis ticks in terms of π
            ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
            ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        else
            # Single theta value - create 1D plot vs φ
            slice_data = md.distribution[p_idx, 1, :]
            # Normalize the data
            slice_data_norm = slice_data ./ maximum(slice_data)
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Momentum Distribution at p = $(round(actual_p, digits=3)) a.u., θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Azimuthal Angle φ (radians)", 
                     ylabel="Normalized Probability Density",
                     title=plot_title)
            
            lines!(ax, md.φ, slice_data_norm, linewidth=2, color=:blue)
            ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
        end
        
    elseif slice_type == :theta_slice
        # Plot at fixed polar angle - since we typically have only one theta value
        θ_idx = 1
        
        # Create 2D slice: momentum magnitude vs azimuthal angle P(p, φ) at fixed θ
        slice_data = md.distribution[:, θ_idx, :]  # Shape: (n_p, n_phi)
        # Normalize the data
        slice_data_norm = slice_data ./ maximum(slice_data)
        
        # Create polar plot
        fig = Figure(size=(800, 800))
        plot_title = isempty(title) ? "Momentum Distribution at θ = $(round(md.θ[θ_idx] * 180/π, digits=1))°, $(pulse_info)" : title
        ax = PolarAxis(fig[1, 1], title=plot_title)
        
        # For polar plots, we need to create a surface plot
        # Create meshgrid for polar coordinates
        φ_mesh = repeat(md.φ', length(md.p), 1)
        p_mesh = repeat(md.p, 1, length(md.φ))
        
        # Use surface plot for polar coordinates
        surface!(ax, φ_mesh, p_mesh, slice_data_norm, colormap=colormap, shading=NoShading)
        
        # Add colorbar
        Colorbar(fig[1, 2], colormap=colormap, colorrange=(0, 1), 
                label="Normalized Probability Density")
        
        # # Set radial limits
        # ax.rmax = maximum(md.p)
        # ax.rmin = 0
        
    elseif slice_type == :phi_slice
        # Plot at fixed azimuthal angle
        φ_idx = argmin(abs.(md.φ .- slice_value))
        actual_φ = md.φ[φ_idx]
        
        if length(md.θ) > 1
            # Create 2D slice: p vs θ
            slice_data = md.distribution[:, :, φ_idx]
            # Normalize the data
            slice_data_norm = slice_data ./ maximum(slice_data)
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Momentum Distribution at φ = $(round(actual_φ * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Polar Angle θ (radians)", 
                     ylabel="Momentum p (a.u.)",
                     title=plot_title)
            
            hm = heatmap!(ax, md.θ, md.p, slice_data_norm, colormap=colormap)
            Colorbar(fig[1, 2], hm, label="Normalized Probability Density")
            ax.xticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        else
            # Single theta value - create 1D plot vs p
            slice_data = md.distribution[:, 1, φ_idx]
            # Normalize the data
            slice_data_norm = slice_data ./ maximum(slice_data)
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Momentum Distribution at φ = $(round(actual_φ * 180/π, digits=1))°, θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Momentum p (a.u.)", 
                     ylabel="Normalized Probability Density",
                     title=plot_title)
            
            lines!(ax, md.p, slice_data_norm, linewidth=2, color=:blue)
        end
        
    else
        throw(ArgumentError("slice_type must be :p_slice, :theta_slice, or :phi_slice"))
    end
    
    # Add grid for non-polar plots
    if slice_type != :theta_slice && hasfield(typeof(ax), :xgridvisible)
        ax.xgridvisible = true
        ax.ygridvisible = true
    end
    
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
