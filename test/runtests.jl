using AttainableRegionTutorial
using Test
using Aqua

@testset "AttainableRegionTutorial.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(AttainableRegionTutorial)
    end
end
