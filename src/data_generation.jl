
#Predefined kinetics and solubility as per tutorial
B(S, M, p) = p[1]*S^p[2]*M^p[3] #nucleation rate, 
G(S, T, p) = p[1]*S^p[2]*exp(-p[3]/8.31441/(T+273.15)) #crystal growth rate [m/s]
Cstar(T) = 3.79e-2*T^2 + 3.77e-1*T + 2.07e1 #solubility [kg/m^3] against temperature [C]

#Define Dataset structure
struct Dataset
    "Vector of residence times [s]."
    τ::Vector{Float64} 
    "Vector of temperatures [°C]."
    T::Vector{Float64} 
    "Vector of feed concentrations [kg m⁻³]."
    Cf::Vector{Float64} 
    "Vector of crystal sizes at which CSD data is available [m]."
    L::Vector{Float64} 
    "Matrix of CSDs [m⁻⁴]; each column is a CSD."
    n::Matrix{Float64} 
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

subscript(x) = join(["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"][digit+1] for digit in digits(x))

function Base.show(io::IO, ::MIME"text/plain", z::SystemSpecification)
    Np = length(z.parameters)

    column_labels_matconst = ["shapefactor", "crystaldensity"]
    column_labels_kintherm = ["solubility", "growthrate", "nucleationrate"]
    column_labels_parameters = ["p"*subscript(i) for i in 1:Np]
    row_labels_parameters = ["True", "Estimated", "Deviation (%)"]

    pt_matconst = pretty_table([z.shapefactor z.crystaldensity];
        column_labels = column_labels_matconst, title = "Material Constants")    
    pt_kintherm = pretty_table([z.solubility z.growthrate z.nucleationrate];
        column_labels = column_labels_kintherm, title = "Kin./Therm. Functions")
    
    dev = (z.estimatedparameters .- z.parameters) ./ z.parameters .* 100
    data = [z.parameters'; z.estimatedparameters'; dev']
    pt_parameters = pretty_table(data;
        column_labels = column_labels_parameters, row_labels = row_labels_parameters, 
            title = "Kinetic Parameters")
end

function Dataset(residencetimes, temperatures, feedconcentrations, 
        system::SystemSpecification = SystemSpecification(); 
        nBins = 101, maxSize = 2000e-6, noiselevel = 0.02)
    #Initialize Dataset, n will be filled during generation of synthetic data, hence 
    #requiring the system variable (where kinetics and thermodynamics are stored)
    L = collect(range(0, maxSize, length = nBins))
    n = zeros(length(L), length(temperatures))
    nnoise = zeros(length(L), length(temperatures))

    #simulate CSD for each operating condition given the kinetics provided in system
    for i in eachindex(temperatures)
        τi = residencetimes[i]
        Ti = temperatures[i]
        Cfi = feedconcentrations[i]
        ni, Cssi = simulateCSD(L, τi, Ti, Cfi, system)
        n[:,i] = ni 
    end

    #add noise to the CSDs
    meanlogn = mean(log.(n))
    for i in eachindex(temperatures)        
        nnoise[:,i] = exp.(log.(n[:,i]) .+ noiselevel.*meanlogn.*(0.5 .- rand(rng, length(L))))
    end

    #generate dataset structure
    data = Dataset(residencetimes, temperatures, feedconcentrations, L, nnoise)
    return data
end

function Base.show(io::IO, ::MIME"text/plain", z::Dataset)
    Nsets = length(z.T)

    # Create a Unicode plot (as a string)
    plotmat = Matrix{String}(undef, 1, Nsets)
    for i in 1:Nsets
        plt = lineplot(z.L*1e6, z.L.^3 .* z.n[:,i] ./ 1e6, width = 10, height = 3,
            xticks = false, yticks = false)
        plot_str = sprint(show, MIME("text/plain"), plt)
        plot_str = join([lstrip(line) for line in split(plot_str, '\n')], '\n')
        plotmat[i] = plot_str
    end
    data = [z.τ'; z.T'; z.Cf'; plotmat]

    footnotes = [
        (:row_label, 1, 2) => "Footnote in column label",
    ]

    source_notes = "τ: residence time [s], T: temperature [°C], Cf: feed concentration [kg m⁻³], L³n(L): volume-weighted CSD [μm⁻¹]"

    column_labels = ["Dataset $i" for i in 1:Nsets]
    row_labels = ["τ", "T", "Cf", "L³n(L)"]
    
    pt_parameters = pretty_table(data;
        column_labels = column_labels, row_labels = row_labels, 
            title = "Datasets", fixed_data_column_widths = 13, line_breaks = true,
            row_label_column_alignment = :l, source_notes = source_notes)
end

# ylabel = "L³n(L) [μm⁻¹]", xlabel = "L [μm]"
#[s] [°C] [kg m⁻³]