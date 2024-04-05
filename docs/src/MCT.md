### Multivariate Creative Telescoping

Let us define the algebra ``D`` 
```math
D = \mathbb{Q}(t)[\boldsymbol{x}]\langle \partial_{t}, \partial_{\boldsymbol{x}}\rangle.
```
where ``\boldsymbol{x} = (x_1,\dots,x_n)``.

Let `` I`` be a holonomic ideal of ``D`` and ``p\in D``.

The function MCT returns an operator ``L\in \mathbb{Q}(t)\langle \partial_t\rangle`` such that ``Lp\in I + \sum_{i=1}^n \partial_i D``.
In particular if ``I`` is the annihilator of a holonomic function ``f`` then ``L`` annihilates the integral 

```math
\int pf \boldsymbol{dx}
```
provided it has natural boundaries.


The signature of the function the following.
```@docs
MCT
```

Assumption: 
 - spol corresponds to the operator ``p``
 - gb is a Gröbner basis of a holonomic subideal of ``I``
 - A corresponds to the algebra ``D`` and its order eliminates ``dt``

It is also advised to use an order of the form ``\text{lex} > \text{grevlex} [\text{polynomial variables}] > \text{grevlex} [\text{differential variables}]``


#### Example

In this example we prove that the integral
```math
I(t) = \int_\gamma \frac{x}{x-t}dx
```
where ``\gamma`` is a loop satisfies the diffential equation ``t\partial_t -1``. 

```jldoctest Quickstart
julia> using MultivariateCreativeTelescoping
```
 We define an algebra of differential operators. The parameter of the integral must be a rational variable and the integration variables must be polynomial variables.
 The order must eliminate ``dt`` and in general it is advised to take an order of the form ``dt > [\text{polynomial variables}] > [\text{differential variables}]``

First we compute an holonomic ideal included in the annihilator of ``1/(x-t)``. For more details see the weyl closure subsection of the other functionnalities section.
```jldoctest Quickstart
julia> A = OreAlg(order = "lex dt x dx",ratdiffvars=(["t"],["dt"]),poldiffvars=(["x"],["dx"]))
Ore algebra

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
