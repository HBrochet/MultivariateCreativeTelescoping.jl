## Gröbner bases 
This package implements two algorithms to compute Gröbner bases, namely F4 and F5. The F5 implementation is based on [Lairez2024](https://arxiv.org/abs/2210.13788).
It is recommended to use the F5 implementation as it is faster.

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

Let ``D`` be the algebra 
```math
 D = \mathbb{Q}(\boldsymbol{t})[\boldsymbol{x}]\langle \partial_{\boldsymbol{t}}, \partial_{\boldsymbol{x}}\rangle
``` 
where ``\boldsymbol{t} =(t_1,\dots,t_m)`` and ``\boldsymbol{x} = (x_1,\dots,x_n)``.

Let ``R`` be the same algebra as with rational coefficients
```math
R = \mathbb{Q}(\boldsymbol{t},\boldsymbol{x})\langle \partial_{\boldsymbol{t}}, \partial_{\boldsymbol{x}}\rangle.
``` 
The Weyl closure of a submodule ``S`` of ``D^r`` for ``r\in\mathbb{N}`` is defined to be the ``D``-module ``R\cdot S\cap D^r``.
This definition generalizes the definition of Tsai given in his PhD [thesis](https://dl.acm.org/doi/10.5555/931963) which corresponds to the case ``m=0``.

Assuming that ``D^r/S`` has finite holonomic rank and that generators gens of ``S`` are given as input of the function weyl\_closure, it returns a Gröbner basis of a ``D``-module ``S'`` such that ``D^r/S'`` is holonomic and ``S' \subset R\cdot S\cap D^r ``. It is not guaranteed that the algorithm returns the full weyl closure.


```@docs
weyl_closure_init
```
```@docs
weyl_closure
```

#### Example
Let ``W_2 = \mathbb{Q}[x,y]\langle \partial_t, \partial_x\rangle`` be the second Weyl algebra and let ``f = \frac{1}{x^2 -y^3}``. The function ``f`` is annihilated by the operators ``\partial_x  (x^2-y^3)`` and ``\partial_y (x^2-y^3)``. But these two equations do not generate the full annihilator of ``f`` in ``W_2``. A third equation can be found using the Weyl closure. 

First we load the package, define the algebra and the two obvious equations annihilating ``f``.
```jldoctest wc
julia> using MultivariateCreativeTelescoping

julia> A = OreAlg(order = "grevlex x y > grevlex dx dy",poldiffvars=(["x","y"],["dx","dy"]))
Ore algebra

julia> gens = [parse_OrePoly("dx*(x^2-y^3)",A),parse_OrePoly("dy*(x^2-y^3)",A)]
Vector of Ore polynomials 
```

Then we can compute its weyl closure
```jldoctest wc
julia> init = weyl_closure_init(A)
WeylClosureInit

julia> gb = weyl_closure(gens,A,init)
Vector of Ore polynomials 
```
And indeed we found a third equation.
```jldoctest wc
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
First a precomputation has to be done for which the user has to choose a value for the parameter sigma.
The larger sigma is, the more relations between elements of  ``\int M`` the algorithm will know.
After this precomputation step is done, it is possible to try to find a smaller representative for an element of the quotient w.r.t. the order of A
and to get every element up to a certain degree that can do not have smaller representatives (at least for the current sigma). 

It is known that for sigma sufficiently large these functions returns respectively the smallest representative and all the irreducible monomials.
The user is encouraged to experiment with this parameter.
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
In this example we take ``S`` equal to a holonomic ideal annihilating the rational function ``\frac{1}{x^3 + y^3 + z^3}``.

First we load the package and define the algebra.
```jldoctest intmod
julia> using MultivariateCreativeTelescoping

julia> s = "(x^3 + y^3 + z^3)"
"(x^3 + y^3 + z^3)"

julia> A = OreAlg(order = "grevlex x y z > grevlex dx dy dz",poldiffvars=(["x","y","z"],["dx","dy","dz"]))
Ore algebra
```

Then we compute a Gröbner basis of ``S``.
```jldoctest intmod
julia> p = parse_OrePoly(s,A)
Ore polynomial

julia> ann = [parse_OrePoly("dx*"*s,A),parse_OrePoly("dy*"*s,A),parse_OrePoly("dz*"*s,A)]
Vector of Ore polynomials 

julia> init = weyl_closure_init(A)
WeylClosureInit

julia> gb = weyl_closure(ann,A,init)
Vector of Ore polynomials 
```

Then comes the precomputation step with sigma equals to 5
```jldoctest intmod
julia> precomp = representative_in_integral_module_precomp(gb,5,A)
RepInIntMod
```

Now we can compute every monomials of degree less than 4 that do not have smaller representatives in the quotient (here sigma is large enough).
```jldoctest intmod
julia> irr = irreducible_monomials(precomp,3,A)
Vector of Ore polynomials 

julia> prettyprint(irr,A)
vector of 1 OrePoly
(1)
```

Indeed we can check that ``x^2`` reduces to ``0`` but ``1`` doesn't. 
```jldoctest intmod
julia> red = representative_in_integral_module(precomp,parse_OrePoly("x^2",A),A)
Ore polynomial

julia> prettyprint(red,A)

julia> red2 = representative_in_integral_module(precomp,parse_OrePoly("1",A),A)
Ore polynomial

julia> prettyprint(red2,A)
(1)
```

