function findCss(Css, system, τ, T, Cf)
    #magma density calculated from mass balance
    M = Cf - Css
    
    #Evaluate kinetic expressions and solubility
    parameters = system.kineticparameters
    nucp = @view parameters[1:3]
    growthp = @view parameters[4:6]
    cstar = system.solubility(T)
    S = log(Css/cstar)
    b = system.nucleationrate(S, M, nucp)
    g = system.growthrate(S, T, growthp)

    #magma density calculated from analytical solution of CSD
    (; crystaldensity, shapefactor) = system 
    M2 = 6 * crystaldensity * shapefactor * b * τ^4 * g^3 

    return M - M2 #difference between magma density calculated through mass balance and magma 
                  #density calculated via CSD must be zero
end

function simulateCSD(L, τ, T, Cf, system)
    (; kineticparameters, solubility, nucleationrate, growthrate) = system
    cstar = solubility(T)
    midPoint = (cstar + Cf) / 2
    #we first attempt solving the nonlinear equation using a bracketing method focussing on 
    #the l.h.s interval (cstar, midPoint). This is to avoid the trivial steady state 
    #(M = 0), which always exists with the kinetics specified in the tutorial.
    f = (Css, parameters) -> findCss(Css, system, τ, T, Cf)
    prob = IntervalNonlinearProblem(f, (cstar, midPoint), system) 
    sol = solve(prob)
    if sol.retcode == SciMLBase.ReturnCode.InitialFailure #In case we did not find Css, consider other half of interval.
        prob = IntervalNonlinearProblem(f, (midPoint, Cf), system)
        sol = solve(prob)
    end
    Css = sol.u

    #Calculate CSD
    M = Cf - Css
    nucp = @view kineticparameters[1:3]
    growthp = @view kineticparameters[4:6]
    S = log(Css/cstar)
    b = nucleationrate(S, M, nucp)
    g = growthrate(S, T, growthp)
    n = b./g.*exp.(-L./g./τ)
    d43 = 4*g*τ #comes from analytical solution of PBE
    return n, Css, d43
end
