"""

Contains a function describing kinetics of a crystallization mechanism and information on 
how many optimizable parameters occur in the function.

## Fields
$FIELDS

## Constructors

```
FunctionalExpression(f::Function, Nparameters::Integer)
```

## Examples
`S` is supersaturation ln(c/c*) where c* is solubility.
`T` is temperature in °C
`M` is normalized suspension density

# Example 1: Growth rate definition
```
g = (S, T, p) -> p[1] * S^p[2] * exp(-p[3] / (T + 273.15)) 
growthrate = FunctionalExpression(g, 3) #3 parameters occur in g
```

# Example 2: Nucleation rate definitions

```
b1 = (S, M, p) -> p[1] * S^p[2] * M^p[3] 
nucleationrate1 = FunctionalExpression(b1, 3) #3 parameters occur in b
```
```
b2 = (S, M, p) -> p[1] * S^p[2]
nucleationrate2 = FunctionalExpression(b2, 2) 
```
Note that function signature must be maintained despite `M` not appearing on r.h.s. of b2.

# Example 3: Solubility line
```
cstar = (T, p) -> 3.79e-2*T^2 + 3.77e-1*T + 2.07e1
solubility = FunctionalExpression(cstar, 0) 
```

"""
struct FunctionalExpression{F}
    "Function describing kinetics"
    f::F
    "Number of optimizable parameters in function"
    Nparameters::Int64
end

function (x::FunctionalExpression)(y...)
    x.f(y[1:end-1]..., y[end])
end

#Predefined kinetics and solubility as per tutorial
B(S, M, p) = p[1] * S^p[2] * M^p[3] #nucleation rate [m⁻³ s ⁻¹]
G(S, T, p) = p[1] * S^p[2] * exp(- p[3] / 8.31441 / (T+273.15)) #crystal growth rate [m/s]
Cstar(T, p) = 3.79e-2*T^2 + 3.77e-1*T + 2.07e1 #solubility [kg/m^3] against temperature [C]

"""

A system specification contains material constants, a solubility expression, crystallization 
kinetics and parameters involved in those expressions. 

## Fields
$FIELDS

## Constructors

```
b = (S, M, p) -> p[1] * S^p[2] * M^p[3]
g = (S, T, p) -> p[1] * S^p[2] * exp(- p[3] / 8.31441 / (T+273.15))
cstar = (T,p) -> 3.79e-2*T^2 + 3.77e-1*T + 2.07e1
SystemSpecification(; crystaldensity = 1200.0,
         crystalshapefactor = pi / 6,
         solubility = FunctionalExpression(cstar, 0),
         nucleationrate = FunctionalExpression(b, 3),
         growthrate = FunctionalExpression(g, 3),
         parameters = [3.34e-4, 1.1, 1.44e4, 3e5, 2.0, 1.6])
```

"""
@kwdef struct SystemSpecification{F1,F2,F3}
    "Density of the crystalline material [kg m⁻³]"
    crystaldensity::Float64 = 1200.0
    "Volume shape factor k, so that volume of a crystal is V = kL³"
    shapefactor::Float64 = pi/6
    "Solubility function [kg m⁻³]"
    solubility::F1 = FunctionalExpression(Cstar, 0)
    "Growth rate  [m s ⁻¹]"
    growthrate::F2 = FunctionalExpression(G, 3)
    "Nucleation rate [kg m⁻³]"
    nucleationrate::F3 = FunctionalExpression(B, 3)
    "Vector of parameters occuring in the growth and nucleation rates"
    parameters::Vector{Float64} = [3.34e-4, 1.1, 1.44e4, 3e5, 2.0, 1.6]
end

#Define Dataset structure
"""

A dataset consists of MSMPR crystallizer operating conditions and associated measurement 
results. For the purpose of this tutorial data is assumed to consist of steady state crystal
size distributions that all possess the same discretization.

## Fields
$FIELDS

## Constructors

To generate artificial data from a model:
```
Dataset(residencetimes, temperatures, feedconcentrations, 
        system::SystemSpecification = SystemSpecification(); 
        nBins = 101, maxSize = 2000e-6, noiselevel = 0.02)
```

To directly supply measured data:
```
Dataset(residencetimes, temperatures, feedconcentrations, particlesizes, CSDs)
```

"""
struct Dataset
    "Vector of residence times [s]"
    τ::Vector{Float64} 
    "Vector of temperatures [°C]"
    T::Vector{Float64} 
    "Vector of feed concentrations [kg m⁻³]"
    Cf::Vector{Float64} 
    "Vector of crystal sizes at which CSD data is available [m]."
    L::Vector{Float64} 
    "Matrix of CSDs [m⁻⁴]; each column is a CSD."
    n::Matrix{Float64} 
end

function Base.show(io::IO, ::MIME"text/plain", z::SystemSpecification)
    Np = length(z.parameters)

    column_labels_matconst = ["shapefactor", "crystaldensity"]
    column_labels_kintherm = ["solubility", "growthrate", "nucleationrate"]
    row_labels_kintherm = ["Function", ["p"*subscript(i) for i in 1:Np]...]

    pt_matconst = pretty_table([z.shapefactor z.crystaldensity];
        column_labels = column_labels_matconst, title = "Material Constants")

    s = z.solubility
    g = z.growthrate
    n = z.nucleationrate
    p = z.parameters

    one = [p[1:s.Nparameters]..., ["" for _ in s.Nparameters+1:Np]...]
    two = [["" for _ in 1:s.Nparameters]..., p[s.Nparameters+1:s.Nparameters+g.Nparameters]..., ["" for _ in s.Nparameters+g.Nparameters+1:Np]...]
    three =  [["" for _ in 1:s.Nparameters+g.Nparameters]..., [p[s.Nparameters+g.Nparameters+1:Np]...]...]
    data_kintherm = [s.f g.f n.f; one two three]
    pt_kintherm = pretty_table(data_kintherm;
        column_labels = column_labels_kintherm, title = "Kin./Therm. Functions", 
        row_labels = row_labels_kintherm, row_label_column_alignment = :l)
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

    source_notes = "τ: residence time [s], T: temperature [°C], Cf: feed concentration [kg m⁻³], L³n(L): volume-weighted CSD [μm⁻¹]"
    column_labels = ["Dataset $i" for i in 1:Nsets]
    row_labels = ["τ", "T", "Cf", "L³n(L)"]
    
    pt_parameters = pretty_table(data;
        column_labels = column_labels, row_labels = row_labels, 
            title = "Datasets", fixed_data_column_widths = 13, line_breaks = true,
            row_label_column_alignment = :l, source_notes = source_notes)
end

function Base.copy(x::SystemSpecification)
    F = propertynames(x)    
    y = SystemSpecification([getproperty(x, f) for f in F]...)
end

#provides multi-digit subscripts programmatically
subscript(x) = join(["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"][digit+1] for digit in reverse(digits(x)))
