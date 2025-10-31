#This code should be read in conjunction with the tutorial chapter in the book.
module AttainableRegionTutorial

using StableRNGs
using NonlinearSolve, GLMakie, Printf
using Statistics 

using Optimization, OptimizationBBO, OptimizationPolyalgorithms, OptimizationOptimJL
using LinearAlgebra

using DocStringExtensions
using PrettyTables
using UnicodePlots

const rng = StableRNG(1234) #for reproducibility of random numbers

#TODO: would be nice to display the system and the dataset objects in a nice format.

function initializeRNG(seed)
    rng = StableRNG(seed)
    return nothing
end

include("CSD_calculation.jl")
include("data_generation.jl")
include("parameter_estimation.jl")
include("attainable_region_generation.jl")
include("plotting.jl")

export SystemSpecification, Dataset, FunctionalExpression, estimate_parameters
export generate_attainable_region, plot_CSDs, plot_fit_quality, plot_attainable_region_dynamic

end
