using AttainableRegionTutorial
using Test

#Initialize System
system = SystemSpecification()

#generate dataset
residencetimes = [15.0, 20.0, 30.0, 60.0, 30.0, 50.0, 60.0, 80.0, 40.0, 40.0].*60 #converted to seconds
temperatures = [20.0, 20.0, 20.0, 20.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
feedconcentrations = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 80.0, 60.0]
dataset = Dataset(residencetimes, temperatures, feedconcentrations, system)

@test dataset isa Dataset

#make plot of data
#Re-create tutorial figure 
fig = plot_CSDs(dataset, system, ids = [1,4,6,10], figure = (fontsize = 20, 
    figure_padding = 30))

#Or plot all data
fig2 = plot_CSDs(dataset, system, mode = :separate, figure = (fontsize = 20, 
    figure_padding = 30))

#Do parameter estimation
estsystem = estimate_parameters(system, dataset)

@test estsystem isa SystemSpecification

#Plot quality of fit
#Re-create tutorial figure 
fig3 = plot_fit_quality(estsystem, dataset, ids = [1,4,6,10], mode = :combined, 
    axislegend = false)

#Or plot all data
fig4 = plot_fit_quality(estsystem, dataset, mode = :separate)

#Generate data for the attainable region
conditions, d43, P = generate_attainable_region(estsystem, npoints_temperature = 50, 
                                                npoints_feed_conc = 50, 
                                                npoints_residence_time = 50)

#Plot attainable region
#static plot for the tutorial.

#dynamic plot for exploring
fig5 = plot_attainable_region_dynamic(system, conditions, d43, P)
