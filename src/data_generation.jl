
#Predefined kinetics and solubility as per tutorial
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

#Define System specification; contains all kinetic parameters, material constants and 
#solubility information
@kwdef struct SystemSpecification{F1,F2,F3}
    crystaldensity::Float64 = 1200.0
    shapefactor::Float64 = pi/6
    solubility::F1 = Cstar
    growthrate::F2 = G
    nucleationrate::F3 = B
    parameters::Vector{Float64} = [3e5, 2.0, 1.6, 3.34e-4, 1.1, 1.44e4]
    estimatedparameters::Vector{Float64} = NaN*ones(6)
end

function generateDataset(; system = SystemSpecification(), nBins = 101, maxSize = 2000e-6, noiselevel = 0.02)
    #Initialize Dataset, Css and n to be filled during generation of synthetic data
    τ = [15.0, 20.0, 30.0, 60.0, 30.0, 50.0, 60.0, 80.0, 40.0, 40.0].*60 #converted to seconds
    T = [20.0, 20.0, 20.0, 20.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    Cf = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 80.0, 60.0]
    Css = zeros(length(T))
    L = collect(range(0, maxSize, length = nBins))
    n = zeros(length(L), length(T))
    nnoise = zeros(length(L), length(T))

    for i in eachindex(T)
        τi = τ[i]
        Ti = T[i]
        Cfi = Cf[i]
        
        #calculate PSD
        ni, Cssi = simulateCSD(L, τi, Ti, Cfi, system)

        #add into data structure
        Css[i] = Cssi
        n[:,i] = ni 
    end

    meanlogn = mean(log.(n))
    for i in eachindex(T)        
        nnoise[:,i] = exp.(log.(n[:,i]) .+ noiselevel.*meanlogn.*(0.5 .- rand(rng, length(L))))
    end

    data = Dataset(τ, T, Cf, Css, L, n, nnoise)
    return data
end
