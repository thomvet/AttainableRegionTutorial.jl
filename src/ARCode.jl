#This code should be read in conjunction with the tutorial chapter in the book.

using NonlinearSolve, GLMakie, Printf, Statistics

#Material constants
rhoc = 1200 #crystal density in kg/m^3
kv = pi/6

#Kinetics and solubility
B(S, M, p) = p[1]*S^p[2]*M^p[3] #nucleation rate, 
G(S, T, p) = p[1]*S^p[2]*exp(-p[3]/8.31441/(T+273.15)) #crystal growth rate [m/s]
Cstar(T) = 3.79e-2*T^2 + 3.77e-1*T + 2.07e1 #solubility [kg/m^3] against temperature [C]

#Define Dataset structure
struct Dataset
    τ::Vector{Float64} #mean residence time [s]
    T::Vector{Float64} #temperature [C]
    Cf::Vector{Float64} #feed concentration [kg/m^3]
    Css::Vector{Float64} #steady state concentration 
    L::Vector{Float64} #crystal length [m]
    n::Matrix{Float64} #noise free crystal size distribution [m^-4]
    nnoise::Matrix{Float64} #crystal size distribution with noise added [m^-4]
end

#Initialize Dataset, Css and n to be filled during generation of synthetic data
τm = [15.0, 20.0, 30.0, 60.0, 30.0, 50.0, 60.0, 80.0, 40.0, 40.0].*60 #converted to seconds
Tm = [20.0, 20.0, 20.0, 20.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
Cfm = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 80.0, 60.0]
Css = zeros(length(Tm))
L = collect(range(0, 2000e-6, length = 101))
n = zeros(length(L), length(Tm))
data = Dataset(τm, Tm, Cfm, Css, L, n, copy(n))

function findCss(Css, parameters, τ, T, Cf)
    #magma density calculated from mass balance
    M = Cf - Css
    
    #Kinetic expressions
    nucp = @view parameters[1:3]
    growthp = @view parameters[4:6]
    cstar = Cstar(T)
    S = log(Css/cstar)
    b = B(S, M, nucp)
    g = G(S, T, growthp)

    #magma density calculated from analytical solution of CSD
    M2 = 6 * rhoc * kv * b * τ^4 * g^3 

    return M - M2 #difference between magma density calculated through mass balance and magma 
                  #density calculated via CSD must be zero
end

#define error function
function errorFun(parameters, data)
    ssq = zero(eltype(parameters))
    L = data.L
    for i in 1:length(data.T)
        n = @view data.nnoise[:,i]
        T = data.T[i]
        τ = data.τ[i]
        Cf = data.Cf[i]
        nsim, Css = simulateCSD(L, τ, T, Cf, parameters)
        ssq = ssq + sum((log(n[j]+1) - log(nsim[j]+1))^2 for j in eachindex(n))
    end
    return ssq
end

#Kinetic Parameters
#p = [kb, b, j, kg, g, E]
parameters = [3e5, 2.0, 1.6, 3.34e-4, 1.1, 1.44e4]

function simulateCSD(L, τ, T, Cf, parameters)
    cstar = Cstar(T)
    midPoint = (cstar + Cf) / 2
    #we first attempt solving the nonlinear equation using a bracketing method focussing on 
    #the l.h.s interval (cstar, midPoint). This is to avoid the trivial steady state 
    #(M = 0), which always exists with the kinetics specified in the tutorial.
    f = (Css, parameters) -> findCss(Css, parameters, τ, T, Cf)
    prob = IntervalNonlinearProblem(f, (cstar, midPoint), parameters) 
    sol = solve(prob)
    if sol.retcode == SciMLBase.ReturnCode.InitialFailure #In case we did not find Css, consider other half of interval.
        prob = IntervalNonlinearProblem(f, (midPoint, Cf), parameters)
        sol = solve(prob)
    end
    Css = sol.u

    #Calculate CSD
    M = Cf - Css
    nucp = @view parameters[1:3]
    growthp = @view parameters[4:6]
    cstar = Cstar(T)
    S = log(Css/cstar)
    b = B(S, M, nucp)
    g = G(S, T, growthp)
    n = b./g.*exp.(-L./g./τ)
    d43 = 4*g*τ

    return n, Css, d43
end

m43 = zeros(length(data.T))
M = zeros(length(data.T))
S = zeros(length(data.T))
for i in 1:length(data.T)
    τ = data.τ[i]
    T = data.T[i]
    Cf = data.Cf[i]
    L = data.L
    #calculate PSD
    n, Css = simulateCSD(L, τ, T, Cf, parameters)
    M[i] = Cf - Css
    S[i] = log(Css/Cstar(T))

    #add into data structure
    data.Css[i] = Css
    data.n[:,i] = n  
end

fig = Figure(fontsize = 24, figure_padding = 30)
fig2 = Figure(fontsize = 24, figure_padding = 30)

ax2 = Axis(fig2[1,1], xlabel = "L [μm]", ylabel = "L³n(L) × 10⁻⁴ [#/μm]", xticks = 0:400:2000, 
            yticks = 0:0.4:2.4, aspect = 1)

colors = to_colormap(:Paired_12)
colors[11] = RGBAf(0.0,0.0,0.0,1.0) #replace pale yellow color with black
for i in 1:length(data.T)
    if i <= 4
        c = i
        r = 1
    elseif i <= 8
        c = i - 4
        r = 2
    else
        c = i - 8
        r = 3
    end

    if i == 1 || i == 5 || i == 9
        yticklabelsvisible = true
        ylabelvisible = true
    else
        yticklabelsvisible = false
        ylabelvisible = false
    end
    n = data.n[:,i]
    meanlogn = mean(log.(data.n))
    nnoise = exp.(log.(n) .+ 0.02.*meanlogn.*(0.5 .- rand(length(n))))
    data.nnoise[:,i] = nnoise

    m43[i] = sum(L.^4 .*nnoise) / sum(L.^3 .*nnoise)
    b = @sprintf "%.0f" m43[i]*1e6
    ax = Axis(fig[r,c], xlabel = "particle size L [μm]", ylabel = "L³n(L) × 10⁻⁴ [#/μm]", 
            title = "Dataset $i, d₄₃ = $b μm", xticks = 0:400:2000, 
            yticks = 0:0.4:2.4, yticklabelsvisible = yticklabelsvisible, 
            ylabelvisible = ylabelvisible)
    #note the units: n has [m^-3], L^3 has [m^3], therefore expressing L^3n(L) in [um^-1] 
    #means dividing by the numbers by 1e6, then we bring them on a nicer intervall by the 1e4.
    lines!(ax, 1e6.*L, L.^3 ./1e6 .*n*1e4, color = colors[i])
    scatter!(ax, 1e6.*L, L.^3 ./1e6 .*nnoise*1e4, color = colors[i]) 
    xlims!(ax, 0, 2000)
    ylims!(ax, 0, 2.4)

    if in(i, [1,4,6,10])
        lines!(ax2, 1e6.*L, L.^3 ./1e6 .*n*1e4, color = colors[i])
        c = @sprintf "Dataset %2.0f" i
        scatter!(ax2, 1e6.*L, L.^3 ./1e6 .*nnoise*1e4, color = colors[i], label = "$c, d₄₃ = $b μm",
            markersize = 10) 
        xlims!(ax2, 0, 2000)
        ylims!(ax2, 0, 2.4)
    end
end

for i in 1:3
    colgap!(fig.layout, i, Relative(0.02))
end

using Optimization, OptimizationPRIMA, OptimizationBBO, ForwardDiff, Optim, OptimizationOptimJL
using OptimizationPolyalgorithms, LikelihoodProfiler, OrdinaryDiffEq
using LinearAlgebra
optfun = OptimizationFunction(errorFun, AutoForwardDiff())
p0 = parameters .* (1.1 .- 0.2*rand(length(parameters)))
lb = [1e2, 0, 0, 1e-10, 0, 0,]
ub = [1e10, 10, 10, 1e-2, 10, 1e5]
prob1 = OptimizationProblem(optfun, p0, data, lb = lb, ub = ub)
# sol1 = solve(prob1, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 1e6)
# prob2 = OptimizationProblem(optfun, p0, data, lb = lb, ub = ub)
# sol2 = solve(prob2, LBFGS())
prob3 = OptimizationProblem(optfun, p0, data)
sol3 = solve(prob3, PolyOpt())
optparameters = sol3.minimizer
#optparameters = sol2.minimizer

plprob = ProfileLikelihoodProblem(prob1, optparameters)
profiler = IntegrationProfiler(integrator = Tsit5(), matrix_type = :hessian)
sol = solve(plprob, profiler)
#Plots.plot(sol)
#calculate confidence intervals
grad = ForwardDiff.gradient(p -> errorFun(p, data), optparameters)
hess = ForwardDiff.hessian((p) -> errorFun(p, data), optparameters)
se = sqrt.(diag(inv(hess)))
rel_se = se./parameters

#plot quality of fit
fig3 = Figure(fontsize = 24)
ax4 = Axis(fig2[1,2], xlabel = "L [μm]", ylabel = "ln(n(L))", xticks = 0:200:2000, yticks = 0:5:40,
    aspect = 1)
for i in 1:length(data.T)
    τ = data.τ[i]
    T = data.T[i]
    Cf = data.Cf[i]
    ndata = data.n[:,i]
    L = data.L
    if i <= 5
        c = i
        r = 1
    else
        c = i - 5
        r = 2
    end
    ax = Axis(fig3[r,c], xlabel = "particle size L [μm]", ylabel = "ln(n(L))", title = "Dataset $i", xticks = 0:200:2000, yticks = 0:5:40)
    scatter!(ax, 1e6.*L[1:2:end], log.(ndata[1:2:end]), color = colors[i])
    nfit, Css, d43 = simulateCSD(L, τ, T, Cf, optparameters)
    lines!(ax, 1e6.*L, log.(nfit), color = colors[i], linestyle = :dash)
    xlims!(ax, 0, 2000)
    ylims!(ax, 10, 40)

    if in(i, [1,4,6,10])
        scatter!(ax4, 1e6.*L[1:2:end], log.(ndata[1:2:end]), color = colors[i])
        lines!(ax4, 1e6.*L, log.(nfit), color = colors[i], linestyle = :dash)
        xlims!(ax4, 0, 2000)
        ylims!(ax4, 15, 35)
    end
end

feed_conc = range(Cstar(5), Cstar(75), length = 200)
residence_time = range(600, 36000, length = 200)
temperatures = range(5, 70, length = 200)
Nconditions = length(feed_conc)*length(residence_time)*length(temperatures)
Pall = NaN*ones(Nconditions)
d43all = NaN*ones(Nconditions)
conditions = zeros(3,Nconditions)

# Generate all possible combinations
i = 1
for fc in feed_conc
    for tau in residence_time
        for temp in temperatures
            conditions[:,i] = [fc, tau, temp]
            i = i + 1
        end
    end
end

L = collect(range(0, 6000e-6, length = 1001))
for i in 1:size(conditions, 2)
    @show i
    Cf_i = conditions[1,i]
    tau_i = conditions[2,i]
    T_i = conditions[3,i]
    if Cf_i >= Cstar(T_i) #run calc. only if feed conc smaller than solubility at temperature.
        n_i, Css_i, d43_i = simulateCSD(L, tau_i, T_i, Cf_i, optparameters)
        d43all[i] = d43_i
        Pall[i] = (conditions[1,i] - Css_i)/tau_i #kg/m3/s
    end       
end

fig4 = Figure(fontsize = 24, size = (1500, 900), figure_padding = 30)
ax5 = Axis(fig4[1,1], xlabel = "productivity, P [kg m⁻³ h⁻¹]", ylabel = "mean particle size, d₄₃ [μm]", 
        aspect = 1, xticks = [10, 250:250:1500...], yticks = 0:150:1200, title = "Attainable region")
ax6 = Axis(fig4[1,2], aspect = 1, xlabel = "particle size L [μm]", 
        ylabel = "L³n(L) × 10⁻⁴ [# μm⁻¹]", title = "Crystal size distribution")
ax7 = Axis(fig4[1,3], aspect = 1, xlabel = "temperature T [°C]", 
        ylabel = "concentration [kg m⁻³]", title = "Operating policy")

plt = scatter!(ax5, Pall*3600, d43all*1e6, markersize = 3, color = :gray48)

#make a reset button
button = Button(fig4[3,2], label = "Reset figure", tellwidth = false) 
on(button.clicks) do n
    for i in 1:5
        nobs[i][] = NaN*zeros(size(L))
        ARobs[i][] = Point2(NaN, NaN)
        OPobs[i][] = NaN*zeros(2,3)
        notify(nobs[i])
        notify(ARobs[i])
        notify(OPobs[i])
        global counter = 0
        foreach(delete!, contents(fig4[2,1:3]))
    end
end

rowsize!(fig4.layout, 1,  Fixed(450))
rowsize!(fig4.layout, 2,  Fixed(80))
rowsize!(fig4.layout, 3,  Fixed(30))


Trange = 0.0:1.0:80.0
lines!(ax7, Trange, Cstar.(Trange), color = :black, label = "solubility")
global counter = 0

nobs = [Observable(NaN*zeros(size(L))) for i in 1:5]
ARobs = [Observable(Point2(NaN, NaN)) for i in 1:5]
OPobs = [Observable(NaN*zeros(2,3)) for i in 1:5]
labels = [0 for i in 1:5]

for i in 1:5
    lines!(ax6, 1e6 .* L, nobs[i], color = Cycled(i))
    scatter!(ax5, ARobs[i], color = Cycled(i), markersize = 10)
    lines!(ax7, OPobs[i], color = Cycled(i))
end
xlims!(ax6, 0, 3000)

on(events(ax5).mousebutton) do ev
    if ev.button == Mouse.left && ev.action == Mouse.press
        p, i = pick(ax5)        
        if p === plt && i > 0
            global counter = rem(counter, 5) + 1
            Cf_i = conditions[1,i]
            tau_i = conditions[2,i]
            T_i = conditions[3,i]
            n_i, Css_i, d43_i = simulateCSD(L, tau_i, T_i, Cf_i, optparameters)
            nobs[counter][] = L.^3 .* n_i ./1e6 .* 1e4
            ARobs[counter][] = Point2(Pall[i]*3600, d43all[i]*1e6)
            prob = IntervalNonlinearProblem((T, p) -> Cstar(T) - Cf_i, extrema(Trange)) 
            sol = solve(prob)
            Tsat = sol.u
            OPobs[counter][] = [Tsat T_i T_i; Cf_i Cf_i Css_i]
            reset_limits!(ax6)        # recompute limits from all axis content
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
            foreach(delete!, contents(fig4[2,1:3]))
            Legend(fig4[2, 1:3],
                elem[2:end],
                labels2,
                patchsize = (35, 35), rowgap = 10, orientation = :horizontal, labelsize = 14,
                tellwidth = false)
        end
    end
    Consume(true)
end
   
xlims!(ax5, 10, 1500)
ylims!(ax5, 0, 1200)
