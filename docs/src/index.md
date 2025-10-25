```@meta
CurrentModule = AttainableRegionTutorial
```

# AttainableRegionTutorial.jl

Welcome to AttainableRegionTutorial.jl! This website hosts documentation and further 
material for the Attainable Region Tutorial appearing in "Chemical Engineering in the 
Pharmaceutical Industry, Third Edition", Chapter 28: "Population Balance Modelling for Batch 
and Continuous Crystallization Processes". Due to appear in 2026.

[Link to second edition](https://www.wiley.com/en-it/Chemical+Engineering+in+the+Pharmaceutical+Industry%3A+Active+Pharmaceutical+Ingredients%2C+2nd+Edition-p-9781119285878).

# Installation

## Julia
To run the tutorial, the [Julia programming language](https://julialang.org/) needs to be 
installed. It is recommended to install ``juliaup``. Alternatively, pre-made installers are 
also available. Details are provided by the [Julia Installation Instructions](https://julialang.org/install/).

Once installed, Julia can be used by running the app directly, by using it from within a 
code editor like [Microsoft Visual Studio Code](https://code.visualstudio.com/), directly 
from the Command Prompt (Windows, cmd.exe) or Terminal (Mac OS), etc.

## Package installation
Package installation is simple. In a Julia session, simply execute:
```
using Pkg
Pkg.add(url = "https://github.com/thomvet/AttainableRegionTutorial.jl")
```
This will add the AttainableRegionTutorial.jl package to the already active environment. 
Alternatively, you could navigate to a desired folder using `cd(path/to/folder)` and use 
`Pkg.activate()` to create a new environment in that folder (or activate the existing 
environment).

More information on environments and package management in Julia can be found [here](https://docs.julialang.org/en/v1/stdlib/Pkg/).
