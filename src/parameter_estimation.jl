@doc raw"""
```
f = errorFun(parameters, (system::SystemSpecification, dataset::Dataset))
```

Objective function that is minimized during parameter estimation. The implementation here 
minimizes the sum of squared differences between simulated and experimental CSDs, `n`:
``f = \sum_i^p\sum_j^q \left(\ln(n_{\mathrm{sim},i,j}+1) - \ln(n_{\mathrm{exp},i,j}+1)\right)^2``
where `p` is the number of datasets and `q` is the number of particle size bins.

"""
function errorFun(parameters, systemdataset)
    system = systemdataset[1]
    dataset = systemdataset[2]
    system.parameters .= parameters
    ssq = zero(eltype(parameters))
    (; L, T, τ, n, Cf)  = dataset
    for i in 1:length(T)
        nnoisei = @view n[:,i]
        Ti = T[i]
        τi = τ[i]
        Cfi = Cf[i]
        nsimi, Css = simulateCSD(L, τi, Ti, Cfi, system)
        ssq = ssq + sum((log(nnoisei[j]+1) - log(nsimi[j]+1))^2 for j in eachindex(L))
    end
    return ssq
end

"""
```
estimatedSystem::SystemSpecification = estimateParameters(system::SystemSpecification, 
    dataset::Dataset; initialguess = nothing, errorFun = errorFun, algorithm = PolyOpt())
```

`system` contains the model whose parameters should be estimated, `dataset` is the data that 
should be matched by the model, `initialguess` is the initial parameter guess, `errorFun` 
is a function returning the objective function value (a metric quantifying the difference 
between the data and the model) and `algorithm` is the optimization algorithm used to 
minimize the objective function.

`initialguess` defaults to `nothing` (initial guess is taken within +/- 10% of the supplied
parameters in `system`) or can be a user-supplied guess (must match length of parameters in 
`system`).

`algorithm` must be one of `PolyOpt()`, `BBO_adaptive_de_rand_1_bin_radiuslimited()`, or 
`LBFGS()`, 

`PolyOpt()` is a Poly-algorithm that explores dafdafa using the [ADAM optimizer](https://arxiv.org/abs/1412.6980) as built into [Optimisers.jl](https://github.com/FluxML/Optimisers.jl) 

`BBO_adaptive_de_rand_1_bin_radiuslimited()` is an evolutionary algorithm from the package [BlackBoxOptim.jl]().

`LBFGS()` is a quasi-Newton (Broyden-Fletcher-Goldfarb-Shanno) algorithm from the package [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/algo/lbfgs/)

"""
function estimateParameters(system, dataset; initialguess = nothing, algorithm = PolyOpt(), 
        errorFun = errorFun)
    if isnothing(initialguess)
        kp = system.parameters
        initialguess = kp .* (1.1 .- 0.2*rand(rng, length(kp)))
    end
    systemcopy = deepcopy(system)
    optfun = OptimizationFunction(errorFun, AutoFiniteDiff())
    lb = [1e2, 0, 0, 1e-10, 0, 0,]
    ub = [1e10, 10, 10, 1e-2, 10, 1e5]

    if typeof(algorithm) == typeof(LBFGS()) || typeof(algorithm) == typeof(BBO_adaptive_de_rand_1_bin_radiuslimited())
        prob = OptimizationProblem(optfun, initialguess, (systemcopy, dataset), lb = lb, ub = ub)
    elseif algorithm == PolyOpt()
        prob = OptimizationProblem(optfun, initialguess, (systemcopy, dataset))
    else
        error("Unknown optimization algorithm specified, valid options are: LBFGS(), PolyOpt(), BBO_adaptive_de_rand_1_bin_radiuslimited().")
    end    

    sol = solve(prob, algorithm)
    systemcopy.parameters .= sol.minimizer
    return systemcopy
end
