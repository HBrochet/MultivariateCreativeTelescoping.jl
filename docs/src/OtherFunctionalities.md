## Gröbner bases 
This package implements two algorithms to compute Gröbner bases, namely f4 and f5. The f5 implementation is based on [this](https://arxiv.org/abs/2210.13788) paper.
It is recommended to use the f5 implementation in priority as it is much faster.

```@docs
f4
```
```@docs
f5
```

#### Example

```jldoctest OtherFunc
julia> using MultivariateCreativeTelescoping

julia>  A = OreAlg(order="grevlex x y > grevlex dx dy",poldiffvars=(["x","y"],["dx","dy"]))
Ore algebra

julia> g = parse_vector_OrePoly("[dx - 2*x, dy + 3*y^2]",A)
Vector of Ore polynomials
    
julia> gb1 = f4(g,A)
Vector of Ore polynomials

julia> gb2 = f5(g,A)
Vector of Ore polynomials

julia> gb1 == gb2
true 

julia> prettyprint(gb1,A)
vector of 2 OrePoly
(1)y^2 + (1//3)dy
(1)x + (-1//2)dx
```


## Weyl closure

Let ``D`` be the following algebra 
```math
 D = \mathbb{Q}(\boldsymbol{t})[\boldsymbol{x}]\langle \partial_{\boldsymbol{t}}, \partial_{\boldsymbol{x}}\rangle
``` 
where ``\boldsymbol{t} =(t_1,\dots,t_m)`` and ``\boldsymbol{x} = (x_1,\dots,x_n)``.

Let ``R`` be the same algebra as with rational coefficients
```math
R = \mathbb{Q}(\boldsymbol{t},\boldsymbol{x})\langle \partial_{\boldsymbol{t}}, \partial_{\boldsymbol{x}}\rangle
``` 
The Weyl closure of a submodule ``S`` of ``D^r`` for ``r\in\mathbb{N}`` is defined to be the ``D``-module ``R\cdot S\cap D^r``.
This definition generalizes the definition of Tsai given in his PhD [thesis](https://dl.acm.org/doi/10.5555/931963).

Assuming that ``D^r/S`` has finite holonomic rank and that generators gens of ``S`` are given as input of the function weyl_closure, it returns a Gröbner basis of a ``D``-module ``S'`` such that ``D^r/S'`` is holonomic and ``S' \subset R\cdot S\cap D^r ``.

```@docs
weyl_closure_init
```
```@docs
weyl_closure
```

#### Example

```jldoctest
julia> using MultivariateCreativeTelescoping

julia> A = OreAlg(order = "grevlex x y > grevlex dx dy",poldiffvars=(["x","y"],["dx","dy"]))
Ore algebra

julia> p = parse_OrePoly("x^2-y^3",A)
Ore polynomial

julia> gens = [parse_OrePoly("dx*(x^2-y^3)",A),parse_OrePoly("dy*(x^2-y^3)",A)]
Vector of Ore polynomials 

julia> init = weyl_closure_init(A)
WeylClosureInit

julia> gb = weyl_closure(gens,A,init)
Vector of Ore polynomials 

julia> prettyprint(f5(gens,A),A)
vector of 2 OrePoly
(1)dxy^3 + (-1)dxx^2
(1)dyy^3 + (-1)dyx^2

julia> prettyprint(gb,A)
vector of 3 OrePoly
(1)dyy^3 + (-1)dyx^2
(3)dxy^2 + (2)dyx
(3)dxx + (2)dyy + (1)
```

## Integral of a module 

Keeping the same notation as in the previous subsection, let ``S`` be a submodule of ``D^r`` and let ``M \stackrel{def}{=} D^r/S``. The integral of the module ``M`` 
is defined as 
```math
\int M = M/\partial M \simeq D^r/ (S + \sum_{i=1}^n \partial_i D^r).
```

This package provides tools to compute in this quotient. 
During the precomputation step the user has to choose a value for the parameter ``sigma``.
The larger sigma is, the more relations between elements of  ``\int M`` the algorithm will know.
After this precomputation step is done, it is possible to try to find a smaller representative for an element of the quotient w.r.t. the order of A
and to get every element up to a certain degree that can do not have smaller representative (at least for the current sigma). 


```@docs
representative_in_integral_module_precomp
``` 
```@docs
representative_in_integral_module
``` 
```@docs
irreducible_monomials
``` 

#### Example
```jldoctest
julia> using MultivariateCreativeTelescoping

julia> s = "(x^3 + y^3 + z^3)"
"(x^3 + y^3 + z^3)"

julia> A = OreAlg(order = "grevlex x y z > grevlex dx dy dz",poldiffvars=(["x","y","z"],["dx","dy","dz"]))
Ore algebra

julia> p = parse_OrePoly(s,A)
Ore polynomial

julia> ann = [parse_OrePoly("dx*"*s,A),parse_OrePoly("dy*"*s,A),parse_OrePoly("dz*"*s,A)]
Vector of Ore polynomials 

julia> init = weyl_closure_init(A)
WeylClosureInit

julia> gb = weyl_closure(ann,A,init)
Vector of Ore polynomials 

julia> precomp = representative_in_integral_module_precomp(gb,5,A)
RepInIntMod

julia> irr = irreducible_monomials(precomp,3,A)
Vector of Ore polynomials 

julia> prettyprint(irr,A)
vector of 1 OrePoly
(1)

julia> red = representative_in_integral_module(precomp,parse_OrePoly("x^2",A),A)
Ore polynomial

julia> prettyprint(red,A)

julia> red2 = representative_in_integral_module(precomp,parse_OrePoly("1",A),A)
Ore polynomial

julia> prettyprint(red2,A)
(1)
```

