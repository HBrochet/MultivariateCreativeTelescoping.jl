```@meta
CurrentModule = MultivariateCreativeTelescoping
DocTestSetup  = quote
    using MultivariateCreativeTelescoping
end
```

The first step is to define an algebra using the OreAlg function.
## Definition of an Algebra 
The following command create an object of type OreAlg.
```@docs
OreAlg
``` 
- char is the characteristic of the base field and is by default $0$.
- ratvars is a vector of string representing rational variables in the base field.
- ratdiffvars is a pair of vectors of same length representing pairs of variables $(t_i,d_{t_i})$ 
  that satisfy the relation $d_{t_i}t_i = t_i d_{t_i} +1$. The variables $t_i$ are part of the base field 
  and the variables $d_{t_i}$ are polynomial variables. 
- poldiffvars is a pair of vectors of same length representing pairs of variables $(x_i,d_{x_i})$ 
  that satisfy the relation $d_{x_i}x_i = x_i d_{x_i} +1$. The variables $x_i$ and $d_{x_i}$ are polynomial variables 
- polvars is a vector of string representing polynomial variables.
- locvars is a pair of vectors of same length representing pairs $(T_i,p_i)$ where $T_i$ are variable names 
  and $p_i$ are parsable polynomials in all the previous variables except for $d_{t_i}$ and $d_{x_i}$. 
  The variables $T_i$ corresponds to the inverse of $p_i$ and satisfy the commutation rule 
  $d_u T = Td_u - \partial_u(p_i)T^2$ for any $u\in \{ x_i,t_i\}$.
- order is a parsable order of the form "ord var1 ... varn > ord var1 ... varm > ..." 
  where ord is currently either lex or grevlex and vari names of the previous variables. 
- nomul is a vector of strings. If a name of a polynomial variable is in this vector then no multiplication by this  
  variable will be performed during Gröbner basis computation.

Here is a simple example
```jldoctest
julia> using MultivariateCreativeTelescoping

julia> A = OreAlg(order = "lex dt > grevlex x dx",ratdiffvars=(["t"],["dt"]),poldiffvars=(["x"],["dx"]))
OreAlgebra
```