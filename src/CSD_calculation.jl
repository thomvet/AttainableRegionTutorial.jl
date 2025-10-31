@doc raw"""
```
diffMT = find_Css(Css, system, τ, T, Cf)
```

calculates the magma/suspension density using the overall mass balance of the MSMPR 
crystallizer (`M`) and the analytical solution of the PBE (at steady state) for the MSMPR 
crystallizer (``M_2``) and returns the difference between the two.

In mathematical form:
``M = C_f - C_{ss}``

``M2 = k_v \rho_c \int\limits_0^\infty L^3n(L) dL = 6 k_v \rho_c B  τ^4 G^3``

where ``k_v`` is the shape factor, ``\rho_c`` is the crystal density, ``B`` is the nucleation
rate, ``G`` is the growth rate, ``C_f`` is the feed concentration,  `τ` is the residence time, 
and ``C_{ss}`` is the steady state concentration.

"""
function find_Css(Css, system, τ, T, Cf)
    #magma density calculated from mass balance
    M = Cf - Css
    
    #Evaluate kinetic expressions and solubility
    parameters = system.parameters
    ps, pg, pn = __indexp(system)
    cstar = system.solubility(T, ps)
    S = log(Css/cstar)
    b = system.nucleation_rate(S, T, M, pn)
    g = system.growth_rate(S, T, pg)

    #magma density calculated from analytical solution of CSD
    (; density, shape_factor) = system 
    M2 = 6 * density * shape_factor * b * τ^4 * g^3 

    return M - M2 #difference between magma density calculated through mass balance and magma 
                  #density calculated via CSD must be zero
end

"""
```
n, Css, d43 = simulate_CSD(L, τ, T, Cf, system::SystemSpecification)
```

Calculates the CSD `n` (discretized at sizes `L`), the steady state concentration `Css`, and 
the volum-weighted mean particle size `d43` occuring in an MSMPR crystallizer at steady 
state when it is operated at temperature `T`, residence time `τ`, and feed concentration 
`Cf`. Kinetics and solubility are specified in `system`.

"""
function simulate_CSD(L, τ, T, Cf, system)
    (; solubility, nucleation_rate, growth_rate) = system
    ps, pg, pn = __indexp(system) 
    cstar = solubility(T, ps)
    midpoint = (cstar + Cf) / 2
    #we first attempt solving the nonlinear equation using a bracketing method focussing on 
    #the l.h.s interval (cstar, midpoint). This is to avoid the trivial steady state 
    #(M = 0), which always exists with the kinetics specified in the tutorial.
    f = (Css, p) -> find_Css(Css, system, τ, T, Cf)
    prob = IntervalNonlinearProblem(f, (cstar, midpoint), system) 
    sol = solve(prob)
    if sol.retcode == SciMLBase.ReturnCode.InitialFailure #In case we did not find Css, consider other half of interval.
        prob = IntervalNonlinearProblem(f, (midpoint, Cf), system)
        sol = solve(prob)
    end
    Css = sol.u

    #Calculate CSD
    M = Cf - Css
    S = log(Css/cstar)
    b = nucleation_rate(S, T, M, pn)
    g = growth_rate(S, T, pg)
    n = b./g.*exp.(-L./g./τ)
    d43 = 4*g*τ #comes from analytical solution of PBE
    return n, Css, d43
end

"""
```
ps, pg, pn = __indexp(system::SystemSpecification)
```

Internal implementation detail. Convenience function that returns the optimizable parameters
in the solubility, the growth rate and the nucleation rate (`ps`, `pg` and `pn`, 
respectively) occurring in `system`.

"""
function __indexp(system)
    parameters = system.parameters
    Ns = system.solubility.Nparameters
    Ng = system.growth_rate.Nparameters
    Nn = system.nucleation_rate.Nparameters
    ps = @view parameters[1:Ns]
    pg = @view parameters[Ns+1:Ns+Ng]
    pn = @view parameters[Ns+Ng+1:Ns+Ng+Nn]
    return ps, pg, pn
end
