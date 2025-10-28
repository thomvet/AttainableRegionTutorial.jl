```@meta
CurrentModule = AttainableRegionTutorial
```

# Book tutorial

```jldoctest tutorial; output = false
using AttainableRegionTutorial

#Initialize System
system = SystemSpecification()

# output

```

```@raw html
<img src="../assets/Tutorial_SystemSpecification.png" alt="REPL Output showing the default system specification used in the tutorial" width="310"/>
```

```jldoctest tutorial; output = false

#generate dataset
residencetimes = [15.0, 20.0, 30.0, 60.0, 30.0, 50.0, 60.0, 80.0, 40.0, 40.0].*60 #converted to seconds
temperatures = [20.0, 20.0, 20.0, 20.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
feedconcentrations = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 80.0, 60.0]
dataset = Dataset(residencetimes, temperatures, feedconcentrations, system)

# output

```
After executing the last command, you will see a summary of the dataset, as shown below. 
Note that the Unicode ("dot") plots show volume-weighted crystal size distributions (CSDs) 
and are meant to provide a quick overview of the CSDs contained within the dataset. They are 
depicted to scale, i.e., broader distributions have lower peaks and lower areas below the 
curves represent lower suspension densities.

![REPL Output showing the datasets](assets/Tutorial_Dataset.png)

While nice, we may want to have some more professional looking plots of the CSDs. We first 
recreate the figure depicting the CSDs as shown in the book chapter.

```jldoctest tutorial; output = false
#make plot of data
#Re-create tutorial figure 
fig = plotCSDs(dataset, system, ids = [1,4,6,10], figure = (fontsize = 20, 
    figure_padding = 30))

# output
Figure()
```

```@raw html
<img src="../assets/Tutorial_CSDs_combined.png" alt="REPL Output showing the CSDs in a combined plot" width="466"/>
```

```jldoctest tutorial; output = false
#Now plot all the data in separate plots
fig2 = plotCSDs(dataset, system, mode = :separate, figure = (fontsize = 20, 
    figure_padding = 30))

# output
Figure()
```

![REPL Output showing the datasets](assets/Tutorial_CSDs_separate.png)

```jldoctest tutorial; output = false
#Do parameter estimation
estsystem = estimateParameters(system, dataset)

#Plot quality of fit
#Re-create tutorial figure 
fig3 = plotFitQuality(estsystem, dataset, ids = [1,4,6,10], 
    mode = :combined, fontsize = 20, figure_padding = 30, legend = false)

#Or plot all data
fig4 = plotFitQuality(estsystem, dataset, mode = :separate, fontsize = 20, figure_padding = 30)

# output

```

![Image showing data fit quality in combined plot](assets/Tutorial_FitQuality.png)


```jldoctest tutorial; output = false
#Generate data for the attainable region
conditions, d43, P = generateAttainableRegion(estsystem, npointstemperature = 50, npointsfeedconc = 50,
        npointsresidencetime = 50)

#Plot attainable region
#static plot for the tutorial.

#dynamic plot for exploring
fig5 = plotAttainableRegion_dynamic(system, conditions, d43, P)

# output
Figure()
```
