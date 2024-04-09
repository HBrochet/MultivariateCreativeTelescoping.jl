The first step is to define an algebra using the OreAlg function.
## Definition of an Ore algebra 
The following command creates an object of type OreAlg.
```@docs
OreAlg
``` 
- char is the characteristic of the base field and is ``0`` by default .
- ratvars is a vector of string representing rational variables in the base field.
- ratdiffvars is a pair of vectors of same length representing pairs of variables ``(t_i,d_{t_i})`` 
  satisfying the relation ``d_{t_i}t_i = t_i d_{t_i} +1``. The variables ``t_i`` are part of the base field 
  and the variables ``d_{t_i}`` are polynomial variables. 
- poldiffvars is a pair of vectors of same length representing pairs of variables ``(x_i,d_{x_i})`` 
  satisfying the relation ``d_{x_i}x_i = x_i d_{x_i} +1``. The variables ``x_i`` and ``d_{x_i}`` are polynomial variables 
- polvars is a vector of string representing polynomial variables.
- locvars is a pair of vectors of same length representing pairs ``(T_i,p_i)`` where ``T_i`` are variable names 
  and ``p_i`` are parsable polynomials in all the previous variables except for ``d_{t_i}`` and ``d_{x_i}``. 
  The variables ``T_i`` corresponds to the inverse of ``p_i`` and satisfy the commutation rule 
  ``d_u T = Td_u - \partial_u(p_i)T^2`` for any ``u\in \{ x_j,t_j\}``.
- order is a parsable order of the form "ord var1 ... varn > ord var1 ... varm > ..." 
  where ord is currently either lex or grevlex and vari names of the previous variables. 
- nomul is a vector of strings. If a name of a polynomial variable is in this vector then no multiplication by this  
  variable will be performed during Gröbner basis computations.


#### Example
The algebra ``\mathbb{Q}(t)[x]\langle \partial_t, \partial_x\rangle`` with the order lex ``dt`` > grevlex ``x`` ``dx`` can be defined with the command
```jldoctest DefOA
julia> using MultivariateCreativeTelescoping

julia> A = OreAlg(order = "lex dt > grevlex x dx",ratdiffvars=(["t"],["dt"]),poldiffvars=(["x"],["dx"]))
Ore algebra
```

## Creating Ore polynomials
The most convenient method to create elements of an algebra A is to use the parse\_OrePoly and parse\_vector\_OrePoly functions. 
```@docs
parse_OrePoly
``` 

```@docs
parse_vector_OrePoly
``` 
#### Example
Continuing the previous example:
```jldoctest DefOA
julia> p = parse_OrePoly("t*dt + x*dx",A)
Ore polynomial

julia> vec = parse_vector_OrePoly("[x*dx, dt*t, dx*(t-x), x^4*dt^5]", A)
Vector of Ore polynomials
```

## Printing Ore polynomials
Printing is done with the prettyprint function based on information stored in the algebra A.
Remark that the package uses a standard monomial basis with the derivatives w.r.t. polynomial variables on the left. This is unusual but
actually more efficient for the computations done in this package.

```@docs
prettyprint
```

#### Example
Continuing the previous example:
```jldoctest DefOA
julia> prettyprint(p,A)
(t)dt + (1)dxx + (-1)

julia> prettyprint(vec,A)
vector of 4 OrePoly
(1)dxx + (-1)
(t)dt + (1)
(-1)dxx + (t)dx
(1)dt^5x^4
```

## Operations on Ore polynomials
Here is a list of basic operations that are supported by this package. Feel free to contact me if you wish to have access to other functions. 


```@docs
add(::OrePoly, ::OrePoly, ::OreAlg)
```
```@docs
sub(::OrePoly, ::OrePoly, ::OreAlg)
```
```@docs
mul(::OrePoly, ::OrePoly, ::OreAlg)
```
```@docs
length(::OrePoly)
```




## Exporting Ore polynomials to other systems
It is possible to convert an OrePoly of a vector of OrePoly to a parsable string which can be parsed by other computer algebra systems.

```@docs
mystring
```

#### Example
Continuing the previous example:
```jldoctest DefOA
julia> mystring(p,A)
"(t)*dt*1 + (1)*dx*x*1 + (-1)*1"

julia> mystring(vec,A)
"[(1)*dx*x*1 + (-1)*1, (t)*dt*1 + (1)*1, (-1)*dx*x*1 + (t)*dx*1, (1)*dt^5*x^4*1]"
```