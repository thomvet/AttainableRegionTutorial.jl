#define error function
function errorFun(parameters, systemdataset)
    system = systemdataset[1]
    dataset = systemdataset[2]
    system.kineticparameters .= parameters
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

function estimateParameters(system, dataset, errorFun = errorFun, initialguess = nothing)
    if isnothing(initialguess)
        kp = system.kineticparameters
        initialguess = kp .* (1.1 .- 0.2*rand(rng,length(kp)))
    end
    ogparameters = copy(system.kineticparameters)
    optfun = OptimizationFunction(errorFun, AutoFiniteDiff())
    lb = [1e2, 0, 0, 1e-10, 0, 0,]
    ub = [1e10, 10, 10, 1e-2, 10, 1e5]
    # prob1 = OptimizationProblem(optfun, p0, dataset, lb = lb, ub = ub)
    # sol1 = solve(prob1, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 1e6)
    # prob2 = OptimizationProblem(optfun, p0, data, lb = lb, ub = ub)
    # sol2 = solve(prob2, LBFGS())
    prob3 = OptimizationProblem(optfun, initialguess, (system, dataset))
    @show initialguess
    sol3 = solve(prob3, PolyOpt())
    optparameters = sol3.minimizer
    #optparameters = sol2.minimizer
    system.kineticparameters .= ogparameters
    system.estimatedparameters .= optparameters
    return optparameters
end

function confidenceIntervalsViaProfiling()
    plprob = ProfileLikelihoodProblem(prob1, optparameters)
    profiler = IntegrationProfiler(integrator = Tsit5(), matrix_type = :hessian)
    sol = solve(plprob, profiler)
    #Plots.plot(sol)
    #calculate confidence intervals
    grad = ForwardDiff.gradient(p -> errorFun(p, data), optparameters)
    hess = ForwardDiff.hessian((p) -> errorFun(p, data), optparameters)
    se = sqrt.(diag(inv(hess)))
    rel_se = se./parameters
end
