#CSD plots 
@doc raw"""
```
plotCSDs(dataset::Dataset, system::Union{Nothing, SystemSpecification} = nothing; 
    ids = 1:length(dataset.T), mode = :combined, figure = (;), axis = (;), plot = (;))
```

Produces a plot of volume-weighted crystal size distributions. Data are plotted as points 
while model outputs are shown as lines.

`dataset` contains the data to be plotted. If providing `system`, the plot will also include 
the model output as comparison to the data. `ids` can be used to specify which datasets are 
included in the plot.

For `mode = :combined` all selected datasets will be shown in a single subfigure whereas 
`mode = :separate` plots each CSD in its own subfigure. 

Figure, Axis and Plot properties can be freely adjusted by providing the keywords `figure`, 
`axis`, and `plot` with a named tuple, e.g., `figure = (fontsize = 20, figure_padding = 30)`. 
Informatoin on valid attributes can be found in the [Makie.jl Documentation](https://docs.makie.org/dev/). 
"""
function plotCSDs(dataset, system = nothing; ids = 1:length(dataset.T), mode = :combined, 
        figure = (fontsize = 20, figure_padding = 30), axis = (;), plot = (;))
    fig = Figure(; figure...)

    if mode == :combined
        ax = Axis(fig[1,1]; xlabel = "particle size L [μm]", ylabel = "volume weighted CSD, L³n(L) × 10⁻⁴ [#/μm]", xticks = 0:400:2000, 
                    yticks = 0:0.4:2.4, aspect = 1, axis...)
    end

    colors = to_colormap(:Paired_12)
    colors[11] = RGBAf(0.0,0.0,0.0,1.0) #replace pale yellow color with black

    (; L, n, T, τ, Cf) = dataset
    for (i, val) in enumerate(ids)
        if mode != :combined
            if i <= 5
                c = i
                r = 1
                xticklabelsvisible = false
                xlabelvisible = false
            else
                c = i - 5
                r = 2
                xticklabelsvisible = true
                xlabelvisible = true
            end
            if i == 1 || i == 6
                    yticklabelsvisible = true
                    ylabelvisible = true
            else
                yticklabelsvisible = false
                ylabelvisible = false
            end
        end
        
        m43 = sum(L.^4 .*n[:,val]) / sum(L.^3 .*n[:,val])
        b = @sprintf "%.0f" m43*1e6

        if mode != :combined
            ax = Axis(fig[r,c], xlabel = "particle size L [μm]", ylabel = "L³n(L) × 10⁻⁴ [#/μm]", 
                    title = "Dataset $val, d₄₃ = $b μm", xticks = 0:400:2000, 
                    yticks = 0:0.4:2.4, yticklabelsvisible = yticklabelsvisible, 
                    ylabelvisible = ylabelvisible, xticklabelsvisible = xticklabelsvisible,
                    xlabelvisible = xlabelvisible, axis...)
        end
        #note the units: n has [m^-3], L^3 has [m^3], therefore expressing L^3n(L) in [um^-1] 
        #means dividing by the numbers by 1e6, then we bring them on a nicer intervall by the 1e4.
        if isa(system, SystemSpecification)
            nsim, _ = simulateCSD(L, τ[val], T[val], Cf[val], system)
            lines!(ax, 1e6.*L, L.^3 ./1e6 .*nsim*1e4; color = colors[val], plot...)
        end
        c = @sprintf "Dataset %2.0f" val
        scatter!(ax, 1e6.*L, L.^3 ./1e6 .*n[:,val]*1e4; color = colors[val], 
            label = "$c, d₄₃ = $b μm", plot...) 
        xlims!(ax, 0, 2000)
        ylims!(ax, 0, 2.4)
    end

    if mode != :combined
        for i in 1:min(length(ids), 4)
            colgap!(fig.layout, i, Relative(0.03))
        end
    else
        axislegend(ax, framevisible = false, position = :rt)
    end
    display(fig)
    return fig
end

function plotFitQuality(system, dataset; mode = :combined, legend = true, ids = 1:length(dataset.T), kws...)
    colors = to_colormap(:Paired_12)
    colors[11] = RGBAf(0.0,0.0,0.0,1.0) #replace pale yellow color with black

    #quality of fit plot
    fig = Figure(; kws...)

    if mode == :combined
        ax = Axis(fig[1,1], xlabel = "particle size L [μm]", ylabel = "ln(n(L))", xticks = 0:400:2000, yticks = 0:5:40,
            aspect = 1)
    end
    (; L, τ, T, Cf, n) = dataset

    for (i, val) in enumerate(ids)
        τi = τ[val]
        Ti = T[val]
        Cfi = Cf[val]
        ni = n[:,val]
        if mode != :combined
            if i <= 5
                c = i
                r = 1
                xticklabelsvisible = false
                xlabelvisible = false
            else
                c = i - 5
                r = 2
                xticklabelsvisible = true
                xlabelvisible = true
            end
            if i == 1 || i == 6
                    yticklabelsvisible = true
                    ylabelvisible = true
            else
                yticklabelsvisible = false
                ylabelvisible = false
            end
        end

        nfiti, Cssi, d43i = simulateCSD(L, τi, Ti, Cfi, system)

        if mode != :combined
            ax = Axis(fig[r,c], xlabel = "particle size L [μm]", ylabel = "ln(n(L))", 
                title = "Dataset $val", xticks = 0:400:2000, yticks = 0:5:40,
                ylabelvisible = ylabelvisible, xlabelvisible = xlabelvisible,
                yticklabelsvisible = yticklabelsvisible, xticklabelsvisible = xticklabelsvisible)
        end
        c = @sprintf "Dataset %2.0f" val
        scatter!(ax, 1e6.*L[1:2:end], log.(ni[1:2:end]), color = colors[val])
        lines!(ax, 1e6.*L, log.(nfiti), color = colors[val], linestyle = :dash, label = "$c, d₄₃ = $(round(d43i*1e6)) μm")
        xlims!(ax, 0, 2000)
        ylims!(ax, 15, 35)        
    end
    if mode == :combined && legend == true
        axislegend(ax, position = :rt, framevisible = false)
    elseif mode != :combined
        for i in 1:min(length(ids), 4)
            colgap!(fig.layout, i, Relative(0.03))
        end    
    end
    return fig
end

function plotAttainableRegion_dynamic(system, conditions, d43, P; kws...)
    fig = Figure(; kws...)
    ax1 = Axis(fig[1,1], xlabel = "productivity, P [kg m⁻³ h⁻¹]", ylabel = "mean particle size, d₄₃ [μm]", 
            aspect = 1, xticks = [10, 250:250:1500...], yticks = 0:150:1200, title = "Attainable region")
    ax2 = Axis(fig[1,2], aspect = 1, xlabel = "particle size L [μm]", 
            ylabel = "L³n(L) × 10⁻⁴ [# μm⁻¹]", title = "Crystal size distribution")
    ax3 = Axis(fig[1,3], aspect = 1, xlabel = "temperature T [°C]", 
            ylabel = "concentration [kg m⁻³]", title = "Operating policy")

    plt = scatter!(ax1, P*3600, d43*1e6, markersize = 3, color = :gray70)

    (; solubility) = system
    ps, _, _ = __indexp(system)

    #make a reset button
    button = Button(fig[3,2], label = "Reset figure", tellwidth = false) 
    on(button.clicks) do n
        for i in 1:5
            nobs[i][] = NaN*zeros(size(L))
            ARobs[i][] = Point2(NaN, NaN)
            OPobs[i][] = NaN*zeros(2,3)
            notify(nobs[i])
            notify(ARobs[i])
            notify(OPobs[i])
            global counter = 0
            foreach(delete!, contents(fig[2,1:3]))
        end
    end

    rowsize!(fig.layout, 1,  Fixed(450))
    rowsize!(fig.layout, 2,  Fixed(80))
    rowsize!(fig.layout, 3,  Fixed(30))

    Text = extrema(conditions[3, :])
    Trange = range(Text[1]-5, Text[2]+5, length = 201)
    lines!(ax3, Trange, solubility.(Trange, Ref(ps)), color = :black, label = "solubility")
    counter = 0

    L = collect(range(0.0, 3000e-6, length = 201))

    nobs = [Observable(NaN*zeros(size(L))) for i in 1:5]
    ARobs = [Observable(Point2(NaN, NaN)) for i in 1:5]
    OPobs = [Observable(NaN*zeros(2,3)) for i in 1:5]
    labels = [0 for i in 1:5]

    for i in 1:5
        lines!(ax2, 1e6 .* L, nobs[i], color = Cycled(i))
        scatter!(ax1, ARobs[i], color = Cycled(i), markersize = 10)
        lines!(ax3, OPobs[i], color = Cycled(i))
    end
    xlims!(ax2, 0, 3000)

    on(events(ax1).mousebutton) do ev
        if ev.button == Mouse.left && ev.action == Mouse.press
            p, i = pick(ax1)        
            if p === plt && i > 0
                counter = rem(counter, 5) + 1
                Cf_i = conditions[1,i]
                tau_i = conditions[2,i]
                T_i = conditions[3,i]
                n_i, Css_i, d43_i = simulateCSD(L, tau_i, T_i, Cf_i, system)
                nobs[counter][] = L.^3 .* n_i ./1e6 .* 1e4
                ARobs[counter][] = Point2(P[i]*3600, d43[i]*1e6)
                prob = IntervalNonlinearProblem((T, p) -> system.solubility(T, ps) - Cf_i, extrema(Trange)) 
                sol = solve(prob)
                Tsat = sol.u
                OPobs[counter][] = [Tsat T_i T_i; Cf_i Cf_i Css_i]
                reset_limits!(ax2)        # recompute limits from all axis content
                notify(nobs[counter])
                labels[counter] = i

                elem = [LegendElement[]]
                labels2 = String[]
                for i in 1:5                
                    if !isnan(ARobs[i][][1])
                        push!(elem, [LineElement(color = Cycled(i), linestyle = nothing),
                        MarkerElement(color = Cycled(i), marker = :circle, markersize = 15,
                        strokecolor = :black)])
                        Cf = conditions[1,labels[i]]
                        tau = conditions[2,labels[i]]
                        T = conditions[3,labels[i]]
                        legendtext = "feed conc: $(round(Cf)) kg/m³\nres. time: $(round(tau/60)) min\ntemp.: $(round(T))°C"
                        push!(labels2, legendtext)
                    end
                end            
                foreach(delete!, contents(fig[2,1:3]))
                Legend(fig[2, 1:3],
                    elem[2:end],
                    labels2,
                    patchsize = (35, 35), rowgap = 10, orientation = :horizontal, labelsize = 14,
                    tellwidth = false)
            end
        end
        Consume(true)
    end
    
    xlims!(ax1, 10, 1500)
    ylims!(ax1, 0, 1200)

    display(fig)
    return fig
end
