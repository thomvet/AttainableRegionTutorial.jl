using AttainableRegionTutorial
using Test
using Aqua

@testset "AttainableRegionTutorial.jl" begin
    @testset "Code quality (Aqua.jl)" begin #basic code quality checks
        Aqua.test_all(AttainableRegionTutorial, persistent_tasks = false)
    end
    # @testset "Constructors" begin #Checks that constructors work as intended.
    #     include("constructors.jl")
    # end
    @testset "Tutorial" begin #check of tutorial as it appears in the documentation
        include("tutorial.jl")
    end
    # @testset "Custom kinetics" begin
    #     include("custom_kinetics.jl")
    # end
    # @testset "Solubility estimation" begin
    #     include("solubility_estimation.jl")
    # end
end
