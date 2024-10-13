# Introduction

MultivariateCreativeTelescoping is a package to compute linear differential equations satisfied by parametric integrals using a method called creative telescoping (see e.g. [here](https://theses.hal.science/tel-01069831) for an introduction to creative telescoping). This package implements an algorithm recently presented at the JNCF conference [[slides]](https://www.cirm-math.fr/RepOrga/3047/Slides/Brochet.pdf) and an article is being written.

This algorithm is dedicated to integrals with natural boundaries of the form 
```math
    I(t) = \int f(t, \boldsymbol{x}) d \boldsymbol{x}
```
where ``\boldsymbol{x}=(x_1,\dots,x_n)`` and ``f`` is a D-finite/holonomic function that satisfies PDEs with coefficients in ``\mathbb{Q}(t,\boldsymbol{x})``.


## Installation

Open Julia in a terminal and execute the following line to install the package.
```
] add MultivariateCreativeTelescoping
```

Then you should be able to load the package in Julia's REPL with the command
```
using MultivariateCreativeTelescoping
```


