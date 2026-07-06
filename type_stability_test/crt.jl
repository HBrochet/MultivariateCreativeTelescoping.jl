using JET
using MultivariateCreativeTelescoping
using Cthulhu 

A = OreAlg(order = "lex x dx", poldiffvars=(["x"],["dx"]))

function foo(p :: OrePoly, q :: OrePoly, A :: OreAlg)
    return add(p,q, A) 
end 

p = parse_OrePoly("x", A)
q = parse_OrePoly("dx", A)

res = compute_with_CRT(foo, A, p, q, A)


@code_warntype compute_with_CRT(foo, A, p, q, A)
@report_opt compute_with_CRT(foo, A, p, q, A)


@descend compute_with_CRT(foo, A, p, q, A)
@descend_code_warntype compute_with_CRT(foo, A, p, q, A)

