#This code should be read in conjunction with the tutorial chapter in the book.
module AttainableRegionTutorial

using StableRNGs
using NonlinearSolve, GLMakie, Printf
using Statistics 

using Optimization, OptimizationBBO, OptimizationPolyalgorithms
using LinearAlgebra

const rng = StableRNG(1234) #for reproducibility of random numbers

function initializeRNG(seed)
    rng = StableRNG(seed)
    return nothing
end

include("CSD_calculation.jl")
include("data_generation.jl")
include("parameter_estimation.jl")
include("attainable_region_generation.jl")
include("plotting.jl")

export SystemSpecification, generateDataset, estimateParameters, generateAttainableRegion
export plotCSDs, plotFitQuality, plotAttainableRegion_dynamic

end
