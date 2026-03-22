using JET
using MultivariateCreativeTelescoping
using Cthulhu 

A = OreAlg(order = "lex dt x dx",ratdiffvars=(["t"],["dt"]), poldiffvars=(["x"],["dx"]),char=primes[1])

p = parse_OrePoly("t",A)

function foo(A::OreAlg,p :: OrePoly)
    cs = [mul(c,c,ctx(A)) for c in coeffs(p)]
    return OrePoly(cs,deepcopy(mons(p)))
end


@code_warntype compute_with_cauchy_interpolation(foo,A,p)
@report_opt compute_with_cauchy_interpolation(foo,A,p)


@descend compute_with_cauchy_interpolation(foo,A,p)
@descend_code_warntype compute_with_cauchy_interpolation(foo,A,p)
