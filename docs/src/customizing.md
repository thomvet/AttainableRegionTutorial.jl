```@meta
CurrentModule = AttainableRegionTutorial
```

# Using custom kinetics and thermodynamics
The package accepts user-supplied kinetic and solubility functions via `FunctionalExpression`.
Below is a minimal example that builds a system with custom kinetics and thermodynamics:

```julia
using AttainableRegionTutorial

# Define custom growth, nucleation and solubility functions
g = (S, T, p) -> p[1] * S^p[2] * exp(-p[3] / (T + 273.15))            # growth: G(S,T; p)
b = (S, T, M, p) -> p[1] * S^p[2] * M^p[3]                           # nucleation: B(S,T,M; p)
cstar = (T, p) -> 3.79e-2 * T^2 + 3.77e-1 * T + 2.07e1               # solubility: c*(T)

# Wrap functions in FunctionalExpression with the number of parameters
growth = FunctionalExpression(g, 3)
nucleation = FunctionalExpression(b, 3)
solubility = FunctionalExpression(cstar, 0)

# Parameter vector: order matches solubility, growth, nucleation Nparameters
params = [3.34e-4, 1.1, 1.44e4, 3e5, 2.0, 1.6]

# Build a SystemSpecification with the custom kinetics
system = SystemSpecification(solubility = solubility,
							 growth_rate = growth,
							 nucleation_rate = nucleation,
							 parameters = params)
```

# Tuning parameter estimation
TODO Lorem ipsum
