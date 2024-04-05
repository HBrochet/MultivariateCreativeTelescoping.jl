```@meta
CurrentModule = MultivariateCreativeTelescoping
DocTestSetup  = quote
    using MultivariateCreativeTelescoping
end
```
Let us demonstrate how this package work on a very simple example based on Cauchy's integral formula. We will check that the integral 
```math
I(t) = \int_\gamma \frac{x}{x-t}dx
```
where ``\gamma`` is a loop satisfies the diffential equation ``t\partial_t -1``. 



 We load the package

```jldoctest Quickstart
julia> using MultivariateCreativeTelescoping
```
 We define an algebra of differential operators. The parameter of the integral must be a rational variable and the integration variables must be polynomial variables.
 The order must eliminate ``dt`` and in general it is advised to take an order of the form ``dt > [\text{polynomial variables}] > [\text{differential variables}]``

```jldoctest Quickstart
julia> A = OreAlg(order = "lex dt x dx",ratdiffvars=(["t"],["dt"]),poldiffvars=(["x"],["dx"]))
Ore algebra
```

We remove the largest polynomial factor in the integrand (here ``x``). Then we define a set of differential equations satisfied by what remains, that is ``1/(x-t)``, and take its weyl closure.
```jldoctest Quickstart
julia> ann = [parse_OrePoly("dt*(x-t)",A), parse_OrePoly("dx*(x-t)",A)]
Vector of Ore polynomials

julia> init = weyl_closure_init(A)
WeylClosureInit

julia> gb = weyl_closure(ann,A,init)
Vector of Ore polynomials 
```
We can now apply the integration algorithm
```jldoctest Quickstart
julia> LDE = MCT(parse_OrePoly("x",A), gb, A)
Ore polynomial 

julia> prettyprint(LDE,A)
(t)dt + (-1)
```
