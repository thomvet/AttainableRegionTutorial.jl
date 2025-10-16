using AttainableRegionTutorial

#Initialize System
system = SystemSpecification()

#generate dataset
residencetimes = [15.0, 20.0, 30.0, 60.0, 30.0, 50.0, 60.0, 80.0, 40.0, 40.0].*60 #converted to seconds
temperatures = [20.0, 20.0, 20.0, 20.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
feedconcentrations = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 80.0, 60.0]
dataset = Dataset(residencetimes, temperatures, feedconcentrations, system)

#make plot of data
#Re-create tutorial figure 
fig = plotCSDs(dataset, ids = [1,4,6,10], fontsize = 20, 
    figure_padding = 30)
#Or plot all data
fig2 = plotCSDs(dataset, mode = :separate, fontsize = 20, 
    figure_padding = 30)

#Do parameter estimation
estimatedParameters = estimateParameters(system, dataset, algorithm)

#Plot quality of fit
#Re-create tutorial figure 
fig3 = plotFitQuality(system, dataset, ids = [1,4,6,10], 
    mode = :combined, fontsize = 20, figure_padding = 30, legend = false)

#Or plot all data
fig4 = plotFitQuality(system, dataset, mode = :separate, fontsize = 20, figure_padding = 30)

#Generate data for the attainable region
conditions, d43, P = generateAttainableRegion(system)

#Plot attainable region
#static plot for the tutorial.

#dynamic plot for exploring
fig5 = plotAttainableRegion_dynamic(system, conditions, d43, P, fontsize = 24, size = (1500, 900), figure_padding = 30)
