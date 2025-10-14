include("src/AttainableRegionTutorial.jl")

using .AttainableRegionTutorial

#Initialize System
system = AttainableRegionTutorial.SystemSpecification()

#generate dataset
dataset = AttainableRegionTutorial.generateDataset(system = system)

#make plot of data
#Re-create tutorial figure 
fig = AttainableRegionTutorial.plotCSDs(dataset, ids = [1,4,6,10], fontsize = 20, 
    figure_padding = 30)
#Or plot all data
fig2 = AttainableRegionTutorial.plotCSDs(dataset, mode = :separate, fontsize = 20, 
    figure_padding = 30)

#Do parameter estimation
estimatedParameters = AttainableRegionTutorial.estimateParameters(system, dataset)

#Plot quality of fit
#Re-create tutorial figure 
fig3 = AttainableRegionTutorial.plotFitQuality(system, dataset, ids = [1,4,6,10], 
    mode = :combined, fontsize = 20, figure_padding = 30, legend = false)

#Or plot all data
fig4 = AttainableRegionTutorial.plotFitQuality(system, dataset, mode = :separate, fontsize = 20, figure_padding = 30)

#Generate data for the attainable region
conditions, d43all, Pall = AttainableRegionTutorial.generateAttainableRegion(system)

#Plot attainable region
#static plot for the tutorial.

#dynamic plot for exploring
fig5 = AttainableRegionTutorial.plotAttainableRegion_dynamic(system, conditions, d43all, Pall, fontsize = 24, size = (1500, 900), figure_padding = 30)
