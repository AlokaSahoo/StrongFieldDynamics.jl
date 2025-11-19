using CairoMakie
# using Zygote

export plot_vector_potential_cartesian, plot_vector_potential_trajectory, plot_electric_field_trajectory
export plot_energy_distribution, plot_angular_distribution, plot_momentum_distribution

"""
    plot_vector_potential_cartesian(pulse::Pulse, t_range::AbstractVector; 
                                   title::String="", xlabel::String="Time (a.u.)", 
                                   ylabel::String="Vector Potential (a.u.)", 
                                   save_path::String="")

Plot the vector potential Ax and Ay of the pulse in a 2D Cartesian plot.

# Arguments
- `pulse::Pulse`: Pulse object containing Ax and Ay functions
- `t_range::AbstractVector`: Time points to evaluate the vector potential
- `title::String=""`: Plot title (auto-generated if empty)
- `xlabel::String="Time (a.u.)"`: X-axis label
- `ylabel::String="Vector Potential (a.u.)"`: Y-axis label
- `save_path::String=""`: Path to save the plot (optional)

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_vector_potential_cartesian(pulse, t_range::AbstractVector; 
                                         title::String="", 
                                         xlabel::String="Time (a.u.)", 
                                         ylabel::String="Vector Potential (a.u.)", 
                                         save_path::String="")
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    
    plot_title = isempty(title) ? 
        "Vector Potential (Ax, Ay) vs Time\nI₀ = $(round(pulse.I₀ * 3.51e16 / 1e14, digits=1))×10¹⁴ W/cm², $(pulse.np) cycles, ϵ = $(pulse.ϵ), helicity = $(pulse.helicity)" : 
        title
    ax.title = plot_title

    Ax_vals = [pulse.Ax(t) for t in t_range]
    Ay_vals = [pulse.Ay(t) for t in t_range]

    lines!(ax, t_range, Ax_vals, linewidth=2, color=:blue, label="Ax(t)")
    lines!(ax, t_range, Ay_vals, linewidth=2, color=:red, label="Ay(t)")

    axislegend(ax, position=:rt)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true

    if !isempty(save_path)
        save(save_path, fig)
        println("Plot saved to: $save_path")
    end

    return fig
end

"""
    plot_vector_potential_trajectory(pulse::Pulse, t_range::AbstractVector; 
                                    title::String="", xlabel::String="Ax (a.u.)", 
                                    ylabel::String="Ay (a.u.)", save_path::String="")

Plot the trajectory of the vector potential in the x-y plane (Ay vs Ax).

# Arguments
- `pulse::Pulse`: Pulse object containing Ax and Ay functions
- `t_range::AbstractVector`: Time points to evaluate the vector potential
- `title::String=""`: Plot title (auto-generated if empty)
- `xlabel::String="Ax (a.u.)"`: X-axis label
- `ylabel::String="Ay (a.u.)"`: Y-axis label
- `save_path::String=""`: Path to save the plot (optional)

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_vector_potential_trajectory(pulse, t_range::AbstractVector; 
                                          title::String="", 
                                          xlabel::String="Ax (a.u.)", 
                                          ylabel::String="Ay (a.u.)", 
                                          save_path::String="")
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    
    plot_title = isempty(title) ? 
        "Vector Potential Trajectory (Ay vs Ax)\nI₀ = $(round(pulse.I₀ * 3.51e16 / 1e14, digits=1))×10¹⁴ W/cm², $(pulse.np) cycles, ϵ = $(pulse.ϵ), helicity = $(pulse.helicity)" : 
        title
    ax.title = plot_title

    Ax_vals = [pulse.Ax(t) for t in t_range]
    Ay_vals = [pulse.Ay(t) for t in t_range]

    lines!(ax, Ax_vals, Ay_vals, linewidth=2, color=:purple, label="A(t) trajectory")

    axislegend(ax, position=:rt)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true

    if !isempty(save_path)
        save(save_path, fig)
        println("Plot saved to: $save_path")
    end

    return fig
end

"""
    plot_electric_field_trajecto(pulse::Pulse, t_range::AbstractVector; 
                                    title::String="", xlabel::String="Ax (a.u.)", 
                                    ylabel::String="Ay (a.u.)", save_path::String="")

Plot the trajectory of the vector potential in the x-y plane (Ay vs Ax).

# Arguments
- `pulse::Pulse`: Pulse object containing Ax and Ay functions
- `t_range::AbstractVector`: Time points to evaluate the vector potential
- `title::String=""`: Plot title (auto-generated if empty)
- `xlabel::String="Ax (a.u.)"`: X-axis label
- `ylabel::String="Ay (a.u.)"`: Y-axis label
- `save_path::String=""`: Path to save the plot (optional)

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_electric_field_trajectory(pulse, t_range::AbstractVector; 
                                          title::String="", 
                                          xlabel::String="Ex (a.u.)", 
                                          ylabel::String="Ey (a.u.)", 
                                          save_path::String="")
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    
    plot_title = isempty(title) ? 
        "Electric Field Trajectory (Ey vs Ex)\nI₀ = $(round(pulse.I₀ * 3.51e16 / 1e14, digits=1))×10¹⁴ W/cm², $(pulse.np) cycles, ϵ = $(pulse.ϵ), helicity = $(pulse.helicity)" : 
        title
    ax.title = plot_title


    # Ex(t) = - pulse.A₀ * (pulse.ω / 2 / pulse.np) * sin(pulse.ω * t / pulse.np) * cos(pulse.ω * t) + pulse.A₀ * pulse.ω * pulse.f(t) * sin((pulse.ω * t))
    # Ey(t) = - pulse.A₀ * pulse.helicity * pulse.ϵ * ((pulse.ω / 2 / pulse.np) * sin(pulse.ω * t / pulse.np) * sin(pulse.ω * t) + pulse.ω * pulse.f(t) * cos((pulse.ω * t)) )
    # Ex_vals = [Ex(t) for t in t_range]
    # Ey_vals = [Ey(t) for t in t_range]

    Ex(t) = Zygote.gradient(pulse.Ax, t)
    Ey(t) = Zygote.gradient(pulse.Ay, t)
    Ex_vals = [-Ex(t)[1] for t in t_range]
    Ey_vals = [-Ey(t)[1] for t in t_range]

    lines!(ax, Ex_vals, Ey_vals, linewidth=2, color=:purple, label="A(t) trajectory")

    axislegend(ax, position=:rt)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true

    if !isempty(save_path)
        save(save_path, fig)
        println("Plot saved to: $save_path")
    end

    return fig
end


"""
    plot_energy_distribution(ed::EnergyDistribution; title::String="", xlabel::String="Energy (a.u.)", 
                            ylabel::String="", save_path::String="", normalize::Bool=true)

Plot the energy distribution of photoelectrons with proper labels.

# Arguments
- `ed::EnergyDistribution`: Energy distribution data to plot
- `title::String=""`: Plot title (auto-generated if empty)
- `xlabel::String="Energy (a.u.)"`: X-axis label
- `ylabel::String=""`: Y-axis label (auto-set based on normalize when empty)
- `save_path::String=""`: Path to save the plot (optional)
- `normalize::Bool=true`: If true, plot normalized distribution; otherwise, plot absolute scale

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_energy_distribution(ed::StrongFieldDynamics.EnergyDistribution; title::String="", xlabel::String="Energy (a.u.)", 
                                 ylabel::String="", save_path::String="", normalize::Bool=true)
    fig = Figure(size=(800, 600))
    # Determine ylabel if not provided
    _ylabel = isempty(ylabel) ? (normalize ? "Normalized Probability" : "Probability") : ylabel
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=_ylabel)
    
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
    
    # Select distribution for plotting (normalized or absolute)
    plot_distribution = if normalize
        if !isempty(ed.distribution) && maximum(ed.distribution) > 0
            ed.distribution ./ maximum(ed.distribution)
        else
            ed.distribution
        end
    else
        ed.distribution
    end
    
    # Plot the energy distribution
    lines!(ax, ed.energies, plot_distribution, linewidth=2, color=:blue, label="Photoelectron Distribution")
    
    # Add grid and formatting
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    
    # Set axis limits with some padding
    if !isempty(plot_distribution) && maximum(plot_distribution) > 0
        ylims!(ax, 0, maximum(plot_distribution) * 1.1)
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
                             plot_type::Symbol=:auto, normalize::Bool=true)

Plot the angular distribution of photoelectrons at fixed energy.

# Arguments
- `ad::AngularDistribution`: Angular distribution data to plot
- `title::String=""`: Plot title (auto-generated if empty)
- `save_path::String=""`: Path to save the plot (optional)
- `plot_type::Symbol=:auto`: Plot type (:auto for automatic selection, :phi_slice for polar plot, :phi_slice for φ slice at fixed θ, :theta_slice for θ slice at fixed φ, :heatmap for 2D plot)

# Returns
- `Figure`: CairoMakie figure object

# Notes
- When plot_type=:auto (default), the function automatically chooses:
  - Polar plot if either θ or φ has only one value
  - Warning and polar plot at first θ value if both have multiple values
"""
function plot_angular_distribution(ad::AngularDistribution; title::String="", save_path::String="", 
                                  plot_type::Symbol=:auto, normalize::Bool=true)
    
    # Auto-generate title if not provided
    plot_title = if isempty(title)
        # Convert intensity back to W/cm²
        I_Wcm2 = ad.pulse.I₀ * 3.51e16
        pulse_info = "I₀ = $(round(I_Wcm2/1e14, digits=1))×10¹⁴ W/cm², $(ad.pulse.np) cycles, ϵ = $(ad.pulse.ϵ), helicity = $(ad.pulse.helicity)"
        "Angular Distribution at E = $(round(ad.energy, digits=3)) a.u., $(pulse_info)"
    else
        title
    end
    
    # Determine plot type automatically if requested
    if plot_type == :auto
        n_theta = length(ad.θ)
        n_phi = length(ad.ϕ)
        
        if n_theta == 1 && n_phi > 1
            # Single θ, multiple φ - use polar plot
            plot_type = :phi_slice
        elseif n_theta > 1 && n_phi == 1
            # Multiple θ, single φ - use theta slice plot
            plot_type = :theta_slice
        elseif n_theta > 1 && n_phi > 1
            # Both multiple - warn and use polar plot with first θ
            @warn "Both θ and φ have multiple values (θ: $n_theta, φ: $n_phi). Using polar plot with first θ value (θ = $(round(ad.θ[1] * 180/π, digits=1))°)."
            plot_type = :phi_slice
        else
            # Both single - use polar plot
            plot_type = :phi_slice
        end
    end
    
    if plot_type == :phi_slice
        # Create polar plot
        fig = Figure(size=(800, 800))
        ax = PolarAxis(fig[1, 1], title=plot_title)
        
        # Determine which data to plot
        if length(ad.θ) == 1
            # Single θ value, plot vs φ
            r_values = if normalize
                m = maximum(ad.distribution[1, :])
                m > 0 ? ad.distribution[1, :] ./ m : ad.distribution[1, :]
            else
                ad.distribution[1, :]
            end
            angles = ad.ϕ
            actual_θ = ad.θ[1]
            subtitle = "at θ = $(round(actual_θ * 180/π, digits=1))°"
        else
            # Multiple θ values, use first one and plot vs φ
            r_values = if normalize
                m = maximum(ad.distribution[1, :])
                m > 0 ? ad.distribution[1, :] ./ m : ad.distribution[1, :]
            else
                ad.distribution[1, :]
            end
            angles = ad.ϕ
            actual_θ = ad.θ[1]
            subtitle = "at θ = $(round(actual_θ * 180/π, digits=1))° (first value)"
        end
        
        # Update title with specific θ information
        ax.title = plot_title * "\n" * subtitle
        
        # Plot as polar line
        lines!(ax, angles, r_values, linewidth=3, color=:blue)
        
    elseif plot_type == :theta_slice
        # Plot θ distribution at specified or middle φ value using polar plot
        fig = Figure(size=(800, 800))
        φ_idx = length(ad.ϕ) == 1 ? 1 : (length(ad.ϕ) ÷ 2 + 1)
        actual_φ = ad.ϕ[φ_idx]
        
        ax = PolarAxis(fig[1, 1], 
                      title="$plot_title\nSlice at φ = $(round(actual_φ * 180/π, digits=1))°")
        
        # Select distribution for plotting
        r_values = if normalize
            m = maximum(ad.distribution[:, φ_idx])
            m > 0 ? ad.distribution[:, φ_idx] ./ m : ad.distribution[:, φ_idx]
        else
            ad.distribution[:, φ_idx]
        end
        
        # Plot as polar line with theta as angles
        lines!(ax, ad.θ, r_values, linewidth=3, color=:red, label="P(θ)")
        
    elseif plot_type == :heatmap
        # Create 2D heatmap: θ vs φ
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], 
                 xlabel="Azimuthal Angle φ (radians)", 
                 ylabel="Polar Angle θ (radians)",
                 title=plot_title)
        
        # Select distribution for visualization
        normalized_distribution = if normalize
            m = maximum(ad.distribution)
            m > 0 ? ad.distribution ./ m : ad.distribution
        else
            ad.distribution
        end
        
    hm = heatmap!(ax, ad.ϕ, ad.θ, normalized_distribution, colormap=:viridis)
    Colorbar(fig[1, 2], hm, label=(normalize ? "Normalized Ionization Probability" : "Ionization Probability"))
        
        # Set axis ticks in terms of π
        ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
        ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        
        # Add grid
        ax.xgridvisible = true
        ax.ygridvisible = true
        
    else
        throw(ArgumentError("plot_type must be :auto, :phi_slice, :phi_slice, :theta_slice, or :heatmap"))
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
                               colormap=nothing, normalize::Bool=true)

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
                                   colormap=nothing, normalize::Bool=true)
    
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
            # Choose data to plot (normalized or absolute)
            slice_to_plot = if normalize
                m = maximum(slice_data)
                m > 0 ? slice_data ./ m : slice_data
            else
                slice_data
            end
            
            fig = Figure(size=(800, 600), fontsize = 20)
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., \n $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Azimuthal Angle φ (radians)", 
                     ylabel="Polar Angle θ (radians)",
                     title=plot_title)
            
            hm = heatmap!(ax, md.φ, md.θ, slice_to_plot', colormap=colormap)
            Colorbar(fig[1, 2], hm, label=(normalize ? "Normalized Ionization Probability" : "Ionization Probability"))
            
            # Set axis ticks in terms of π
            ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
            ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        else
            # Single theta value - create 1D plot vs φ
            slice_data = md.distribution[p_idx, 1, :]
            # Choose data to plot (normalized or absolute)
            slice_to_plot = if normalize
                m = maximum(slice_data)
                m > 0 ? slice_data ./ m : slice_data
            else
                slice_data
            end
            
            fig = Figure(size=(800, 600), fontsize = 20)
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., \n θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Azimuthal Angle φ (radians)", 
                     ylabel=(normalize ? "Normalized Ionization Probability" : "Ionization Probability"),
                     title=plot_title)
            
            lines!(ax, md.φ, slice_to_plot, linewidth=2, color=:blue)
            ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
        end
        
    elseif slice_type == :theta_slice
        # Plot at fixed polar angle - since we typically have only one theta value
        θ_idx = 1
        
        # Create 2D slice: momentum magnitude vs azimuthal angle P(p, φ) at fixed θ
        slice_data = md.distribution[:, θ_idx, :]  # Shape: (n_p, n_phi)
        # Choose data to plot (normalized or absolute)
        slice_to_plot = if normalize
            m = maximum(slice_data)
            m > 0 ? slice_data ./ m : slice_data
        else
            slice_data
        end
        
        # Create polar plot
        fig = Figure(size=(800, 800), fontsize = 20)
        plot_title = isempty(title) ? "Photoelectron Momentum Distribution \n θ = $(round(md.θ[θ_idx] * 180/π, digits=1))°, $(pulse_info)" : title
        ax = PolarAxis(fig[1, 1], title=plot_title)
        
        # For polar plots, we need to create a surface plot
        # Create meshgrid for polar coordinates
        φ_mesh = repeat(md.φ', length(md.p), 1)
        p_mesh = repeat(md.p, 1, length(md.φ))
        
        # Use surface plot for polar coordinates
    surface!(ax, φ_mesh, p_mesh, slice_to_plot, colormap=colormap, shading=NoShading)
        
        # Add colorbar
    if normalize
        Colorbar(fig[1, 2], colormap=colormap, colorrange=(0, 1), 
            label="Normalized Ionization Probability")
    else
        Colorbar(fig[1, 2], colormap=colormap, colorrange=(minimum(slice_data), maximum(slice_data)), 
            label="Ionization Probability")
    end
        
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
            # Choose data to plot (normalized or absolute)
            slice_to_plot = if normalize
                m = maximum(slice_data)
                m > 0 ? slice_data ./ m : slice_data
            else
                slice_data
            end
            
            fig = Figure(size=(800, 600), fontsize = 20)
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution \n φ = $(round(actual_φ * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Polar Angle θ (radians)", 
                     ylabel="Momentum p (a.u.)",
                     title=plot_title)
            
            hm = heatmap!(ax, md.θ, md.p, slice_to_plot, colormap=colormap)
            Colorbar(fig[1, 2], hm, label=(normalize ? "Normalized Ionization Probability" : "Ionization Probability"))
            ax.xticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        else
            # Single theta value - create 1D plot vs p
            slice_data = md.distribution[:, 1, φ_idx]
            # Choose data to plot (normalized or absolute)
            slice_to_plot = if normalize
                m = maximum(slice_data)
                m > 0 ? slice_data ./ m : slice_data
            else
                slice_data
            end
            
            fig = Figure(size=(800, 600), fontsize = 20)
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution \n φ = $(round(actual_φ * 180/π, digits=1))°, θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Momentum p (a.u.)", 
                     ylabel=(normalize ? "Normalized Ionization Probability" : "Ionization Probability"),
                     title=plot_title)
            
            lines!(ax, md.p, slice_to_plot, linewidth=2, color=:blue)
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

#= Plots without normalizing/absolute scale
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
- `colormap=:thermal`: Colormap for heatmaps

# Returns
- `Figure`: CairoMakie figure object
"""
function plot_momentum_distribution(md::MomentumDistribution; title::String="", save_path::String="", 
                                   slice_type::Symbol=:theta_slice, slice_value::Float64=1.0,
                                   colormap=nothing)
    
    # Define custom colormap - use built-in Makie colormap
    custom_colormap = :thermal  # Options: :viridis, :plasma, :inferno, :turbo, :thermal, etc.
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
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Azimuthal Angle φ (radians)", 
                     ylabel="Polar Angle θ (radians)",
                     title=plot_title)
            
            hm = heatmap!(ax, md.φ, md.θ, slice_data', colormap=colormap)
            Colorbar(fig[1, 2], hm, label="Ionization Probability")
            
            # Set axis ticks in terms of π
            ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
            ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        else
            # Single theta value - create 1D plot vs φ
            slice_data = md.distribution[p_idx, 1, :]
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Azimuthal Angle φ (radians)", 
                     ylabel="Ionization Probability",
                     title=plot_title)
            
            lines!(ax, md.φ, slice_data, linewidth=2, color=:blue)
            ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
        end
        
    elseif slice_type == :theta_slice
        # Plot at fixed polar angle - since we typically have only one theta value
        θ_idx = 1
        
        # Create 2D slice: momentum magnitude vs azimuthal angle P(p, φ) at fixed θ
        slice_data = md.distribution[:, θ_idx, :]  # Shape: (n_p, n_phi)
        
        # Create polar plot
        fig = Figure(size=(800, 800))
        plot_title = isempty(title) ? "Photoelectron Momentum Distribution θ = $(round(md.θ[θ_idx] * 180/π, digits=1))°, $(pulse_info)" : title
        ax = PolarAxis(fig[1, 1], title=plot_title)
        
        # For polar plots, we need to create a surface plot
        # Create meshgrid for polar coordinates
        φ_mesh = repeat(md.φ', length(md.p), 1)
        p_mesh = repeat(md.p, 1, length(md.φ))
        
        # Use surface plot for polar coordinates
        surface!(ax, φ_mesh, p_mesh, slice_data, colormap=colormap, shading=NoShading)
        
        # Add colorbar with actual data range
        Colorbar(fig[1, 2], colormap=colormap, colorrange=(minimum(slice_data), maximum(slice_data)), 
                label="Ionization Probability")
        
    elseif slice_type == :phi_slice
        # Plot at fixed azimuthal angle
        φ_idx = argmin(abs.(md.φ .- slice_value))
        actual_φ = md.φ[φ_idx]
        
        if length(md.θ) > 1
            # Create 2D slice: p vs θ
            slice_data = md.distribution[:, :, φ_idx]
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Polar Angle θ (radians)", 
                     ylabel="Momentum p (a.u.)",
                     title=plot_title)
            
            hm = heatmap!(ax, md.θ, md.p, slice_data, colormap=colormap)
            Colorbar(fig[1, 2], hm, label="Ionization Probability")
            ax.xticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
        else
            # Single theta value - create 1D plot vs p
            slice_data = md.distribution[:, 1, φ_idx]
            
            fig = Figure(size=(800, 600))
            plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
            ax = Axis(fig[1, 1], 
                     xlabel="Momentum p (a.u.)", 
                     ylabel="Ionization Probability",
                     title=plot_title)
            
            lines!(ax, md.p, slice_data, linewidth=2, color=:blue)
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
=#
#=
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
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Azimuthal Angle φ (radians)", 
                            ylabel="Polar Angle θ (radians)",
                            title=plot_title)
                   
                   hm = heatmap!(ax, md.φ, md.θ, slice_data', colormap=colormap)
                   Colorbar(fig[1, 2], hm, label="Ionization Probability")
                   
                   # Set axis ticks in terms of π
                   ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
                   ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
               else
                   # Single theta value - create 1D plot vs φ
                   slice_data = md.distribution[p_idx, 1, :]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Azimuthal Angle φ (radians)", 
                            ylabel="Ionization Probability",
                            title=plot_title)
                   
                   lines!(ax, md.φ, slice_data, linewidth=2, color=:blue)
                   ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
               end
               
           elseif slice_type == :theta_slice
               # Plot at fixed polar angle - since we typically have only one theta value
               θ_idx = 1
               
               # Create 2D slice: momentum magnitude vs azimuthal angle P(p, φ) at fixed θ
               slice_data = md.distribution[:, θ_idx, :]  # Shape: (n_p, n_phi)
               
               # Create polar plot
               fig = Figure(size=(800, 800))
               plot_title = isempty(title) ? "Photoelectron Momentum Distribution θ = $(round(md.θ[θ_idx] * 180/π, digits=1))°, $(pulse_info)" : title
               ax = PolarAxis(fig[1, 1], title=plot_title)
               
               # For polar plots, we need to create a surface plot
               # Create meshgrid for polar coordinates
               φ_mesh = repeat(md.φ', length(md.p), 1)
               p_mesh = repeat(md.p, 1, length(md.φ))
               
               # Use surface plot for polar coordinates
               surface!(ax, φ_mesh, p_mesh, slice_data, colormap=colormap, shading=NoShading)
               
               # Add colorbar with automatic color range
               Colorbar(fig[1, 2], colormap=colormap, colorrange=(minimum(slice_data), maximum(slice_data)), 
                       label="Ionization Probability")
               
           elseif slice_type == :phi_slice
               # Plot at fixed azimuthal angle
               φ_idx = argmin(abs.(md.φ .- slice_value))
               actual_φ = md.φ[φ_idx]
               
               if length(md.θ) > 1
                   # Create 2D slice: p vs θ
                   slice_data = md.distribution[:, :, φ_idx]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Polar Angle θ (radians)", 
                            ylabel="Momentum p (a.u.)",
                            title=plot_title)
                   
                   hm = heatmap!(ax, md.θ, md.p, slice_data, colormap=colormap)
                   Colorbar(fig[1, 2], hm, label="Ionization Probability")
                   ax.xticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
               else
                   # Single theta value - create 1D plot vs p
                   slice_data = md.distribution[:, 1, φ_idx]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Momentum p (a.u.)", 
                            ylabel="Ionization Probability",
                            title=plot_title)
                   
                   lines!(ax, md.p, slice_data, linewidth=2, color=:blue)
               end
               
           else
               throw(ArgumentError("slice_type must be :p_slice, :theta_slice, or :phi_slice"))
           end
julia> function plot_momentum_distribution(md::MomentumDistribution; title::String="", save_path::String="", 
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
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Azimuthal Angle φ (radians)", 
                            ylabel="Polar Angle θ (radians)",
                            title=plot_title)
                   
                   hm = heatmap!(ax, md.φ, md.θ, slice_data', colormap=colormap)
                   Colorbar(fig[1, 2], hm, label="Ionization Probability")
                   
                   # Set axis ticks in terms of π
                   ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
                   ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
               else
                   # Single theta value - create 1D plot vs φ
                   slice_data = md.distribution[p_idx, 1, :]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Azimuthal Angle φ (radians)", 
                            ylabel="Ionization Probability",
                            title=plot_title)
                   
                   lines!(ax, md.φ, slice_data, linewidth=2, color=:blue)
                   ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
               end
               
           elseif slice_type == :theta_slice
               # Plot at fixed polar angle - since we typically have only one theta value
               θ_idx = 1
               
               # Create 2D slice: momentum magnitude vs azimuthal angle P(p, φ) at fixed θ
               slice_data = md.distribution[:, θ_idx, :]  # Shape: (n_p, n_phi)
               
               # Create polar plot
               fig = Figure(size=(800, 800))
               plot_title = isempty(title) ? "Photoelectron Momentum Distribution θ = $(round(md.θ[θ_idx] * 180/π, digits=1))°, $(pulse_info)" : title
               ax = PolarAxis(fig[1, 1], title=plot_title)
               
               # For polar plots, we need to create a surface plot
               # Create meshgrid for polar coordinates
               φ_mesh = repeat(md.φ', length(md.p), 1)
               p_mesh = repeat(md.p, 1, length(md.φ))
               
               # Use surface plot for polar coordinates
               surface!(ax, φ_mesh, p_mesh, slice_data, colormap=colormap, shading=NoShading)
               
               # Add colorbar with automatic color range
               Colorbar(fig[1, 2], colormap=colormap, colorrange=(minimum(slice_data), maximum(slice_data)), 
                       label="Ionization Probability")
               
           elseif slice_type == :phi_slice
               # Plot at fixed azimuthal angle
               φ_idx = argmin(abs.(md.φ .- slice_value))
               actual_φ = md.φ[φ_idx]
               
               if length(md.θ) > 1
                   # Create 2D slice: p vs θ
                   slice_data = md.distribution[:, :, φ_idx]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Polar Angle θ (radians)", 
                            ylabel="Momentum p (a.u.)",
                            title=plot_title)
                   
                   hm = heatmap!(ax, md.θ, md.p, slice_data, colormap=colormap)
                   Colorbar(fig[1, 2], hm, label="Ionization Probability")
                   ax.xticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
               else
                   # Single theta value - create 1D plot vs p
                   slice_data = md.distribution[:, 1, φ_idx]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Momentum p (a.u.)", 
                            ylabel="Ionization Probability",
                            title=plot_title)
                   
                   lines!(ax, md.p, slice_data, linewidth=2, color=:blue)
               end
               
           else
               throw(ArgumentError("slice_type must be :p_slice, :theta_slice, or :phi_slice"))
           end
julia> function plot_momentum_distribution(md::MomentumDistribution; title::String="", save_path::String="", 
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
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Azimuthal Angle φ (radians)", 
                            ylabel="Polar Angle θ (radians)",
                            title=plot_title)
                   
                   hm = heatmap!(ax, md.φ, md.θ, slice_data', colormap=colormap)
                   Colorbar(fig[1, 2], hm, label="Ionization Probability")
                   
                   # Set axis ticks in terms of π
                   ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
                   ax.yticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
               else
                   # Single theta value - create 1D plot vs φ
                   slice_data = md.distribution[p_idx, 1, :]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution p = $(round(actual_p, digits=3)) a.u., θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Azimuthal Angle φ (radians)", 
                            ylabel="Ionization Probability",
                            title=plot_title)
                   
                   lines!(ax, md.φ, slice_data, linewidth=2, color=:blue)
                   ax.xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"])
               end
               
           elseif slice_type == :theta_slice
               # Plot at fixed polar angle - since we typically have only one theta value
               θ_idx = 1
               
               # Create 2D slice: momentum magnitude vs azimuthal angle P(p, φ) at fixed θ
               slice_data = md.distribution[:, θ_idx, :]  # Shape: (n_p, n_phi)
               
               # Create polar plot
               fig = Figure(size=(800, 800))
               plot_title = isempty(title) ? "Photoelectron Momentum Distribution θ = $(round(md.θ[θ_idx] * 180/π, digits=1))°, $(pulse_info)" : title
               ax = PolarAxis(fig[1, 1], title=plot_title)
               
               # For polar plots, we need to create a surface plot
               # Create meshgrid for polar coordinates
               φ_mesh = repeat(md.φ', length(md.p), 1)
               p_mesh = repeat(md.p, 1, length(md.φ))
               
               # Use surface plot for polar coordinates
               surface!(ax, φ_mesh, p_mesh, slice_data, colormap=colormap, shading=NoShading)
               
               # Add colorbar with automatic color range
               Colorbar(fig[1, 2], colormap=colormap, colorrange=(minimum(slice_data), maximum(slice_data)), 
                       label="Ionization Probability")
               
           elseif slice_type == :phi_slice
               # Plot at fixed azimuthal angle
               φ_idx = argmin(abs.(md.φ .- slice_value))
               actual_φ = md.φ[φ_idx]
               
               if length(md.θ) > 1
                   # Create 2D slice: p vs θ
                   slice_data = md.distribution[:, :, φ_idx]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Polar Angle θ (radians)", 
                            ylabel="Momentum p (a.u.)",
                            title=plot_title)
                   
                   hm = heatmap!(ax, md.θ, md.p, slice_data, colormap=colormap)
                   Colorbar(fig[1, 2], hm, label="Ionization Probability")
                   ax.xticks = ([0, π/4, π/2, 3π/4, π], ["0", "π/4", "π/2", "3π/4", "π"])
               else
                   # Single theta value - create 1D plot vs p
                   slice_data = md.distribution[:, 1, φ_idx]
                   
                   fig = Figure(size=(800, 600))
                   plot_title = isempty(title) ? "Photoelectron Momentum Distribution φ = $(round(actual_φ * 180/π, digits=1))°, θ = $(round(md.θ[1] * 180/π, digits=1))°, $(pulse_info)" : title
                   ax = Axis(fig[1, 1], 
                            xlabel="Momentum p (a.u.)", 
                            ylabel="Ionization Probability",
                            title=plot_title)
                   
                   lines!(ax, md.p, slice_data, linewidth=2, color=:blue)
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
=#