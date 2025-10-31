"""
```
generate_attainable_region(system::SystemSpecification; temperature_range = (5.0, 75.0), 
        feed_conc_range = nothing,
        residence_time_range = (600.0, 36000.0), npoints_temperature = 200, npoints_feed_conc = 200,
        npoints_residence_time = 200)
```

TODO: complete docstring.

"""
function generate_attainable_region(system; temperature_range = (5.0, 75.0), 
        feed_conc_range = nothing,
        residence_time_range = (600.0, 36000.0), npoints_temperature = 200, npoints_feed_conc = 200,
        npoints_residence_time = 200)
    #if no feed_conc_range given, then use temperature range and solubility to set range
    if isnothing(feed_conc_range)
        ps, _, _ = __indexp(system)
        feed_conc_range = (system.solubility(temperature_range[1], ps), system.solubility(temperature_range[2], ps))
    end

    feed_conc = range(feed_conc_range..., length = npoints_feed_conc)
    residence_time = range(residence_time_range..., length = npoints_residence_time)
    temperatures = range(temperature_range..., length = npoints_temperature)
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
        #TODO: maybe implement progress bar? This could also be done in parallel?
        Cf_i = conditions[1,i]
        tau_i = conditions[2,i]
        T_i = conditions[3,i]
        if Cf_i >= system.solubility(T_i, ps) #run calc. only if feed conc smaller than solubility at temperature.
            n_i, Css_i, d43_i = simulate_CSD(L, tau_i, T_i, Cf_i, system)
            d43all[i] = d43_i
            Pall[i] = (conditions[1,i] - Css_i)/tau_i #kg/m3/s
        end       
    end
    return conditions, d43all, Pall
end
