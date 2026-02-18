```@meta
CurrentModule = AttainableRegionTutorial
```

# Tutorial

## Overview
In this tutorial we will detail how attainable regions can be computed for a single mixed 
suspension mixed product removal (MSMPR) crystallizer at steady state. For the sake of 
simplicity, we will do this based on synthetic data, though originating from case study 
published by [Power et al. (2015)](https://www.sciencedirect.com/science/article/abs/pii/S0009250915001207)
where experimental data was used to establish kinetics for the paracetamol crystallizing from 
isopropanol/water mixtures in an MSMPR.

We proceed by the following steps:
* Background on the MSMPR crystallizer model
* Specification of model solubility, kinetics and material constants and generation of synthetic data
* Parameter estimation
* Generating the attainable region and exploring it

## Background: Mixed Suspension Mixed Product Removal Crystallizer Model
 The population balance equation describing the crystal size distribution in an MSMPR 
crystallizer governed by nucleation and growth can be written as:
```math
\frac{\partial n(L,t)}{\partial t} = -\frac{\partial G(S,T)n(L,t)}{\partial L} - \frac{n(L,t)}{\tau}
```
where ``n(L,t)`` is the crystal size distribution, such that ``n(L,t)\mathrm{d}L`` is the 
number of crystals with sizes between ``L`` and ``L+\mathrm{d}L`` per volume of solution 
(it therefore has units of [m``^{-4}``]). ``t`` is time, ``L`` is crystal size, ``G(S,T)`` 
is the supersaturation- (``S``) and temperature-dependent (``T``) growth rate 
(in [m s``^{-1}``]) and ``\tau`` is the mean residence time of the crystallizer (determined 
by the ratio of the volume of solution inside the crystallizer, ``V``, and the flow rate of 
the particle-free inlet, ``Q`` [m``^3`` s``^{-1}``]). With the typical assumption that nucleation 
occurs at negligibly small particle size (treated as ``L=0``), a boundary condition 
``n(L=0,t) = B(S,T,M_\mathrm{T})/G(S,T)`` completes the description of the crystal phase in 
the crystallizer. ``B(S,T,M_\mathrm{T})`` is the nucleation rate in [m``^{-3}`` s``^{-1}``] 
and ``M_\mathrm{T}`` is the magma/suspension density in [kg m``^{-3}``].

At steady state and with a growth rate that is not size-dependent, this becomes:
```math
0 = -G(S,T)\frac{\mathrm{d} n(L)}{\mathrm{d} L} - \frac{n(L)}{\tau}
```
This is a separable first order differential equation and thus can be directly solved to 
obtain:
```math
n(L) = \frac{B(S,T,M_\mathrm{T})}{G(S,T)}\exp\left(-\frac{L}{G\tau}\right)
```
We define the supersaturation, ``S`` as:
```math
S = \ln \left(\frac{C_\mathrm{ss}}{C_\star(T)}\right)
```
where ``C_\star(T)`` is the solubility. Therefore, the nucleation and growth rates need to 
be calculated at the steady state concentration in the crystallizer. This couples the liquid 
phase and the crystalline phase to each other via the overall mass balance:
```math
C_\mathrm{f} = M_\mathrm{T} + C_\mathrm{ss} 
```
where ``C_\mathrm{f}`` is the feed concentration and ``C_\mathrm{ss}`` is the steady state 
concentration, both in [kg m``^{-3}``]. The suspension density can be calculated from the 
steady state crystal size distribution:
```math
\begin{align*}
M_\mathrm{T} &=&& k_\mathrm{v}\rho_\mathrm{c}\int\limits_0^\infty L^3 n(L) \mathrm{d} L \\
             &=&& k_\mathrm{v}\rho_\mathrm{c}\int\limits_0^\infty L^3 \frac{B(S,T,M_\mathrm{T})}{G(S,T)}\exp\left(-\frac{L}{G\tau}\right) \mathrm{d} L \\
             &=&& 6k_\mathrm{v}\rho_\mathrm{c} τ^4 B G^3 
\end{align*} 
```
where ``\rho_\mathrm{c}`` is the density of the crystals in [kg m``^{-3}``] and ``k_\mathrm{v}`` 
is the volume shape factor of the crystals, defined so that ``V_\mathrm{c} = k_\mathrm{v}L^3`` 
where ``V_\mathrm{c}`` is the volume of a crystal with characteristic size ``L``. This means, 
for example, that ``k_\mathrm{v} = 1`` for cubes with side length ``L``, 
``k_\mathrm{v} = \pi/6`` for spheres of diameter ``L``, 
``k_\mathrm{v} = \left(6\sqrt{2}\right)^{-1}`` for regular tetrahedra with edge length 
``L``, etc.

See [Solving root finding problem](@ref) for details on how AttainableRegionTutorial.jl 
solves the mass balance and population balance equation numerically.

## Generation of synthetic data using a specified model

We start out by defining a description of the "system" we are crystallizing. This includes 
defining material constants, such as the crystal density and the volumetric shape factor, as 
well as defining an expression of the solubility against temperature and the kinetics of 
nucleation and crystal growth. As outlined in the book chapter, we will use the kinetics 
and solubility established by [Power et al. (2015)](https://www.sciencedirect.com/science/article/abs/pii/S0009250915001207) 
for paracetamol crystallizing from isopropanol/water mixtures.

Solubility in [kg m``^3``] against temperature in [°C]:
```math
C_\star(T) = 3.79 \times 10^{-2} T^2 + 3.77 \times 10^{-1} T + 2.07 \times 10^{1} 
```

#nucleation rate [m⁻³ s ⁻¹]
```math
B(S, T, M_\mathrm{T}) = p[1] * S^p[2] * \left(\frac{M_\mathrm{T}}\right)^p[3] 
```
#crystal growth rate [m/s]
```math
G(S, T) = p[1] * S^p[2] * exp(- p[3] / 8.31441 / (T + 273.15)) 
```

These kinetics and the solubility line are built into AttainableRegionTutorial.jl by default 
and specifying them has therefore been simplified as much as possible. We can do it by 
running the following two lines of Julia code:

```julia
using AttainableRegionTutorial

#Initialize System
system = SystemSpecification()
```
You should see the following output in your Julia REPL. Note that the (tunable) parameters 
in the kinetic expressions have been numbered consecutively and that it is clear which 
parameter occurs in which functional expression. 

```@raw html
<img src="../assets/Tutorial_SystemSpecification.png" alt="REPL Output showing the default system specification used in the tutorial" width="310"/>
```

We next generat

```julia
#generate dataset
residence_times = [15.0, 20.0, 30.0, 60.0, 30.0, 50.0, 60.0, 80.0, 40.0, 40.0].*60 #converted to seconds
temperatures = [20.0, 20.0, 20.0, 20.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
feed_concentrations = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 80.0, 60.0]
dataset = Dataset(residence_times, temperatures, feed_concentrations, system)
```

After executing the last command, you will see a summary of the dataset, as shown below. 
Note that the Unicode ("dot") plots show volume-weighted crystal size distributions (CSDs) 
and are meant to provide a quick overview of the CSDs contained within the dataset. They are 
depicted to scale, i.e., broader distributions have lower peaks and lower areas below the 
curves represent lower suspension densities.

![REPL Output showing the datasets](./assets/Tutorial_Dataset.png)

While nice, we may want to have some more professional looking plots of the CSDs. We first 
recreate the figure depicting the CSDs as shown in the book chapter.

```julia
#make plot of data
#Re-create tutorial figure 
fig = plot_CSDs(dataset, system, ids = [1,4,6,10])
```

```@raw html
<img src="../assets/Tutorial_CSDs_combined.png" alt="REPL Output showing the CSDs in a combined plot" width="466"/>
```

```julia
#Now plot all the data in separate plots
fig2 = plot_CSDs(dataset, system, mode = :separate, figure = (fontsize = 20, 
    figure_padding = 30))
```

![REPL Output showing the datasets](./assets/Tutorial_CSDs_separate.png)

## Estimating parameters: fitting a model to the synthetic data
We will now estimate the kinetic parameters from the dataset. For this purpose, we are using 
the function `estimate_parameters()`. We will run this with the default values in this case,
but see the page [Tuning parameter estimation](@ref) for ways to tune the parameter 
estimation procedure. Also consult the docstring of the `estimate_parameters` function by 
writing `?estimate_parameters` in your Julia session.

```julia
estsystem = estimate_parameters(system, dataset)
```
`estsystem` should be a new system specification that contains the estimated parameters, which 
will be shown in the Julia session. It should look like this:
```@raw html
<img src="../assets/Tutorial_EstimatedSystem.png" alt="REPL Output showing the system 
specification after parameter estimation" width="310"/>
```

Note that the precise values may be slightly different based on the optimizer used when 
estimating the parameters and also on the precise noise added when generating the dataset 
above.

During the estimation procedure you will likely encounter warnings like this:
```julia
┌ Warning: The interval is not an enclosing interval, opposite signs at the boundaries are required.
└ @ BracketingNonlinearSolve C:\Users\thvt\.julia\packages\BracketingNonlinearSolve\Ucmxz\src\itp.jl:87
```
These are nothing to be concerned about. They stem for the particular way the root finding 
problem required to solve the PBE and mass balance together is implemented. For a deeper 
explanation, see [Solving the root finding problem](@ref)

```julia
#Plot quality of fit
#Re-create tutorial figure 
fig3 = plot_fit_quality(estsystem, dataset, ids = [1,4,6,10], mode = :combined, 
    axislegend = false)
```

![Image showing data fit quality in combined plot](./assets/Tutorial_FitQuality.png)

```julia
#Or plot all data
fig4 = plot_fit_quality(estsystem, dataset, mode = :separate)
```

TODO: insert plot with the separate fit qualities.

## Generating an attainable region
We can now proceed to populate the attainable region by varying process conditions. The 
`generate_attainable_region` function allows to do this in an easy manner. Specifically, we 
vary the temperature in the crystallizer, the feed concentration, as well as the residence 
time. By default, the function selects reasonable ranges for these values for this case 
study. More information on this can be obtained by typing `?generate_attainable_region` in your 
Julia session or directly here: [`generate_attainable_region`](@ref). 
The function will return volume-weighted mean particle sizes at each operating point, `d43` 
in [m], as well as productivity values `P` in [``\mathrm{kg} \mathrm{m}^{-3} \mathrm{s}^{-1}``].

```julia
#Generate data for the attainable region
conditions, d43, P = generate_attainable_region(estsystem, npoints_temperature = 50, npoints_feed_conc = 50,
        npoints_residence_time = 50)
```
Note that we have intentionally reduced the number of points investigated, so that the code 
finishes running quicker. The output in the book chapter was generated using ``200 \times 200 \times 200`` 
points. 

We can now plot all attainable combinations of mean particle sizes against productivity 
values in a static plot.

```julia
#Plot attainable region
#static plot for the tutorial.
#TODO: add code for fig5
```

In order to explore which operating condition generates which point in the attainable region, 
there is a dynamic analysis tool available:
```julia
#dynamic plot for exploring
fig6 = plot_attainable_region_dynamic(system, conditions, d43, P)
```

`fig6` is an interactive figure. By left-clicking on specific points in the attainable 
(left) region plot, the volume-weighted CSD and the operating policy generating it are 
displayed in the middle and right plot, respetively. Multiple points can be compared using 
multiple left-clicks (up to a maximum of 5 - more clicks result in previously selecting point 
getting deselected). The functionality is showcased in the video below:

TODO: insert video showing the dynamic functionality