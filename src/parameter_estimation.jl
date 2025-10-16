#define error function
function errorFun(parameters, systemdataset)
    system = systemdataset[1]
    dataset = systemdataset[2]
    system.parameters .= parameters
    ssq = zero(eltype(parameters))
    (; L, T, τ, nnoise, Cf)  = dataset
    for i in 1:length(T)
        nnoisei = @view nnoise[:,i]
        Ti = T[i]
        τi = τ[i]
        Cfi = Cf[i]
        nsimi, Css = simulateCSD(L, τi, Ti, Cfi, system)
        ssq = ssq + sum((log(nnoisei[j]+1) - log(nsimi[j]+1))^2 for j in eachindex(L))
    end
    return ssq
end

function estimateParameters(system, dataset; algorithm = LBFGS(), errorFun = errorFun, initialguess = nothing)
    if isnothing(initialguess)
        kp = system.parameters
        initialguess = kp .* (1.1 .- 0.2*rand(rng, length(kp)))
    end
    ogparameters = copy(system.parameters)
    optfun = OptimizationFunction(errorFun, AutoFiniteDiff())
    lb = [1e2, 0, 0, 1e-10, 0, 0,]
    ub = [1e10, 10, 10, 1e-2, 10, 1e5]

    if typeof(algorithm) == typeof(LBFGS()) || typeof(algorithm) == typeof(BBO_adaptive_de_rand_1_bin_radiuslimited())
        prob = OptimizationProblem(optfun, initialguess, (system, dataset), lb = lb, ub = ub)
    elseif algorithm == PolyOpt()
        prob = OptimizationProblem(optfun, initialguess, (system, dataset))
    else
        error("Unknown optimization algorithm specified, valid options are: LBFGS(), PolyOpt(), BBO_adaptive_de_rand_1_bin_radiuslimited().")
    end    

    sol = solve(prob, algorithm)
    optparameters = sol.minimizer
    system.parameters .= ogparameters
    system.estimatedparameters .= optparameters
    return optparameters
end
