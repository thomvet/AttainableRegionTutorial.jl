```@meta
CurrentModule = AttainableRegionTutorial
```

# How it works
This page will contain more in depth explanations on how algorithms are implemented "under 
the hood".

## Solving the root finding problem
As detailed in [Background: Mixed Suspension Mixed Product Removal Crystallizer Model](@ref) 
the MSMPR crystallizer model at steady state contains a population balance equation and the 
mass balance that are coupled to each other and must be solved at the same time. Since 
analytical solutions are not generally available (depends on the kinetic expressions), this 
can be solved numerically.

In AttainableRegionTutorial.jl this is implemented as a root finding problem.

## Generating attainable regions
There is at least two ways of generating attainable regions once a suitable process model is 
established. The first one was showcased in [Tutorial](@ref)

