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
fig = plotCSDs(dataset, system, ids = [1,4,6,10], figure = (fontsize = 20, 
    figure_padding = 30))
#Or plot all data
fig2 = plotCSDs(dataset, system, mode = :separate, figure = (fontsize = 20, 
    figure_padding = 30))

#Do parameter estimation
estsystem = estimateParameters(system, dataset)

#Plot quality of fit
#Re-create tutorial figure 
fig3 = plotFitQuality(estsystem, dataset, ids = [1,4,6,10], mode = :combined, 
    axislegend = false)

#Or plot all data
fig4 = plotFitQuality(estsystem, dataset, mode = :separate)

#Generate data for the attainable region
conditions, d43, P = generateAttainableRegion(estsystem)

#Plot attainable region
#static plot for the tutorial.

#dynamic plot for exploring
fig5 = plotAttainableRegion_dynamic(system, conditions, d43, P)
