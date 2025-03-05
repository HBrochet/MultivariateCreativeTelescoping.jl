```@meta
CurrentModule = MultivariateCreativeTelescoping
DocTestSetup  = quote
    using MultivariateCreativeTelescoping
end
```

This section describes concrete problems that have been solved using this package. 

## Counting k-regular graphs 

A graph is $k$ regular if all its vertices have degree exactly $k$. 
Let $u(n,k)$ be the number of $k$-regular graphs on $n$ vertices. 
We are interested in finding a recurrence relation on  $u(n,k)$ for fixed values of $k$. 
Chyzak and Mishna were able to compute such a recurrence for $k$ up to 7 [[link]](https://arxiv.org/abs/2406.04753). The computation for k=7 took 24 weeks of CPU time. 
With this package, we were able to obtain the same result in 7 minutes and to compute the recurrence for k=8. 

#### Example 
Below is the script for k=6. It returns a linear differential equation annihilating the generating series 
of the $u(n,6)$. 
Excluding compilation time this computation takes approximately 3 seconds and returns an operator 
of order 6 and degree 213.

```julia
using MultivariateCreativeTelescoping
Rpt = OreAlg(order = "lex dt > grevlex p1 p2 p3 p4 p5 p6 > grevlex d1 d2 d3 d4 d5 d6",
                   ratdiffvars = (["t"], ["dt"]),
                   poldiffvars = (["p1", "p2", "p3", "p4", "p5", "p6"], ["d1", "d2", "d3", "d4", "d5", "d6"])
                  )
Gstring = "[0+(1/6*t*p1^3-1/120*t*p1^5-t*p1+1/3*t*p3-1/5*t*p5-1/6*t*p1^2*p3+1/4*t*p1*p4-1/2*t*p1*p2+1/6*t*p2*p3+1/12*t*p1^3*p2-1/8*t*p1*p2^2+p1)*1+(-1/6*t*p3-1/12*t*p1^3+1/2*t*p1+1/4*t*p1*p2)*d1^2+(-t)*d2*d3+(1/2*t)*d1^2*d3+(-1/2*t*p1)*d2^2+(1/2*t)*d1*d2^2+(-1/6*t)*d1^3*d2+(t+1/2*t*p2-1/2*t*p1^2)*d1*d2+(1/120*t)*d1^5+(-1/6*t-1/12*t*p2+1/12*t*p1^2)*d1^3+(t-1/4*t*p4+1/8*t*p2^2+1/2*t*p2+1/24*t*p1^4-1/2*t*p1^2+1/3*t*p1*p3-1/4*t*p1^2*p2)*d1+(1/3*t*p3+1/6*t*p1^3-t*p1-1/2*t*p1*p2)*d2+(-t-1/2*t*p2+1/2*t*p1^2)*d3+(t*p1)*d4+(t)*d5+(-1/24*t*p1)*d1^4+(-t*p1)*d1*d3+(1/2*t*p1)*d1^2*d2+(-t)*d1*d4, 0+(-1/2*t*p2+1/4*t*p4-1/8*t*p2^2+1/2*t*p1^2-1/24*t*p1^4-1/3*t*p1*p3+1/4*t*p1^2*p2+p2-t)*1+(1/2*t+1/4*t*p2-1/4*t*p1^2)*d1^2+(-1/2*t)*d2^2+(-t*p1)*d1*d2+(1/6*t*p1)*d1^3+(1/3*t*p3+1/6*t*p1^3-t*p1-1/2*t*p1*p2)*d1+(-t-1/2*t*p2+1/2*t*p1^2)*d2+(t*p1)*d3+(t)*d4+(-1/24*t)*d1^4+(-t)*d1*d3+(1/2*t)*d1^2*d2, 0+(p3-1/6*t*p1^3+t*p1+1/2*t*p1*p2-1/3*t*p3)*1+(-1/2*t*p1)*d1^2+(-t)*d1*d2+(1/6*t)*d1^3+(-t-1/2*t*p2+1/2*t*p1^2)*d1+(t*p1)*d2+(t)*d3, 0+(p4-1/2*t*p1^2+t+1/2*t*p2)*1+(-1/2*t)*d1^2+(t*p1)*d1+(t)*d2, 0+(t)*d1+(-t*p1+p5)*1, 0+(p6-t)*1, 0+(1+1/3*p1*p3-1/18*p3^2+1/2*p2-1/4*p4+1/6*p6+1/8*p2^2-1/2*p1^2-1/8*p2*p4-1/5*p1*p5-1/720*p1^6+1/48*p1^4*p2-1/16*p1^2*p2^2-1/18*p1^3*p3+1/24*p1^4-1/4*p1^2*p2+1/48*p2^3+1/6*p1*p2*p3+1/8*p1^2*p4)*1+(-1)*d1*d5+(-1)*d2*d4+(1/2)*d1^2*d4+(-1/2)*d3^2+(1/6)*d2^3+(-1/720)*d1^6+(1/8*p4-1/2-1/16*p2^2-1/48*p1^4-1/6*p1*p3+1/8*p1^2*p2-1/4*p2+1/4*p1^2)*d1^2+(-p1)*d2*d3+(1/2*p1)*d1^2*d3+(1/2+1/4*p2-1/4*p1^2)*d2^2+(1/2*p1)*d1*d2^2+(-1/6*p1)*d1^3*d2+(1/2*p1*p2+p1-1/3*p3-1/6*p1^3)*d1*d2+(-1/6)*d1^3*d3+(1/120*p1)*d1^5+(1/18*p3+1/36*p1^3-1/12*p1*p2-1/6*p1)*d1^3+(1/5*p5+p1+1/2*p1*p2-1/4*p1*p4+1/120*p1^5-1/6*p2*p3+1/6*p1^2*p3+1/8*p1*p2^2-1/12*p1^3*p2-1/3*p3-1/6*p1^3)*d1+(-1/4*p4+1/2*p2+1+1/8*p2^2+1/24*p1^4+1/3*p1*p3-1/4*p1^2*p2-1/2*p1^2)*d2+(1/3*p3+1/6*p1^3-1/2*p1*p2-p1)*d3+(-1/2*p2+1/2*p1^2-1)*d4+(p1)*d5+(1)*d6+(-1/4)*d1^2*d2^2+(1/24)*d1^4*d2+(1)*d1*d2*d3+(1/48*p2-1/48*p1^2+1/24)*d1^4+(1/2*p2-1/2*p1^2+1)*d1*d3+(1)*dt+(-1/4*p2+1/4*p1^2-1/2)*d1^2*d2+(-p1)*d1*d4]"

G = parse_vector_OrePoly(Gstring, Rpt)
B = f5(G, Rpt)
res = MCT(one(Rpt), B, Rpt)
```

