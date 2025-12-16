```@meta
CurrentModule = MultivariateCreativeTelescoping
DocTestSetup  = quote
    using MultivariateCreativeTelescoping
end
```

## Quickstart 


Let us demonstrate how this package works on a very simple example based on Cauchy's integral formula. We will check that the integral 
```math
I(t) = \int_\gamma \frac{x}{x-t}dx
```
where ``\gamma`` is a loop satisfies the differential equation ``t\partial_t -1`` for any ``t`` inside that loop. This reflects the well-known fact ``I(t) = 2i\pi t``.




 We load the package

```jldoctest Quickstart
julia> using MultivariateCreativeTelescoping
```

We define the algebra 
```math
 \mathbb{Q}(t)[x]\langle \partial_t, \partial_x\rangle
```
with the order lex ``dt`` > grevlex ``x`` ``dx``. 

```jldoctest Quickstart
julia> A = OreAlg(order = "lex dt > grevlex x dx",ratdiffvars=(["t"],["dt"]),poldiffvars=(["x"],["dx"]))
Ore algebra
```

We remove the largest polynomial factor in the integrand (here ``x``). This is done for efficiency reasons.
 Then we define a set of differential equations satisfied by what remains, that is ``1/(x-t)``. 
They must generate a D-finite ideal.
```jldoctest Quickstart
julia> ann = [parse_OrePoly("dt*(x-t)",A), parse_OrePoly("dx*(x-t)",A)]
Vector of Ore polynomials
```
These equations do not necessarily define a holonomic ideal hence we have to take the Weyl closure.
Note that the WeylClosure command returns only an holonomic approximation.
```jldoctest Quickstart
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
