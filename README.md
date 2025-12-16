# MultivariateCreativeTelescoping.jl

**MultivariateCreativeTelescoping.jl** is a Julia package for computing linear differential equations satisfied by parametric integrals using *creative telescoping* (see, e.g., [this reference](https://theses.hal.science/tel-01069831) for an introduction).

The package implements:
- the integration algorithm introduced in [this paper](https://arxiv.org/abs/2504.12724), and  
- an algorithm for approximating the Weyl closure, presented in my PhD thesis ([link](https://hbrochet.github.io/articles/thesis.pdf)).

Combined, these algorithms make it possible to compute integrals with natural boundaries of the form
```math
I(t) = \int f(t,\boldsymbol{x})\,\mathrm{d}\boldsymbol{x},
```

where $\boldsymbol{x} = (x_1,\dots,x_n)$ and $f$ is a D-finite/holonomic function satisfying a system of partial differential equations with coefficients in $\mathbb{Q}(t,\boldsymbol{x})$.

The package is still under development, but an early version is already available.  
The documentation can be found [here](https://hbrochet.github.io/MultivariateCreativeTelescoping.jl/dev/).