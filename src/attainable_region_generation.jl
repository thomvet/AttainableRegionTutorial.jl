"""
```
generateAttainableRegion(system; temperaturerange = (5.0, 75.0), 
        feedconcrange = nothing,
        residencetimerange = (600.0, 36000.0), npointstemperature = 200, npointsfeedconc = 200,
        npointsresidencetime = 200)
```

TODO: complete docstring.

"""
function generateAttainableRegion(system; temperaturerange = (5.0, 75.0), 
        feedconcrange = nothing,
        residencetimerange = (600.0, 36000.0), npointstemperature = 200, npointsfeedconc = 200,
        npointsresidencetime = 200)
    #if no feedconcrange given, then use temperature range and solubility to set range
    if isnothing(feedconcrange)
        ps, _, _ = __indexp(system)
        feedconcrange = (system.solubility(temperaturerange[1], ps), system.solubility(temperaturerange[2], ps))
    end

    feed_conc = range(feedconcrange..., length = npointsfeedconc)
    residence_time = range(residencetimerange..., length = npointsresidencetime)
    temperatures = range(temperaturerange..., length = npointstemperature)
    Nconditions = length(feed_conc)*length(residence_time)*length(temperatures)
    Pall = NaN*ones(Nconditions)
    d43all = NaN*ones(Nconditions)
    conditions = zeros(3,Nconditions)

    # Generate all possible combinations
    i = 1
    for fc in feed_conc
        for tau in residence_time
            for temp in temperatures
                conditions[:,i] = [fc, tau, temp]
                i = i + 1
            end
        end
    end

    L = collect(range(0, 6000e-6, length = 1001))
    for i in 1:size(conditions, 2)
        @show i
        Cf_i = conditions[1,i]
        tau_i = conditions[2,i]
        T_i = conditions[3,i]
        if Cf_i >= system.solubility(T_i, ps) #run calc. only if feed conc smaller than solubility at temperature.
            n_i, Css_i, d43_i = simulateCSD(L, tau_i, T_i, Cf_i, system)
            d43all[i] = d43_i
            Pall[i] = (conditions[1,i] - Css_i)/tau_i #kg/m3/s
        end       
    end
    return conditions, d43all, Pall
end
