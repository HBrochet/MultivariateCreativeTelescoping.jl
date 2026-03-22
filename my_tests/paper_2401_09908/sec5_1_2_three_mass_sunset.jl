include(joinpath(@__DIR__, "_common.jl"))

# Section 5.1.2 of arXiv:2401.09908
# Two-loop sunset, three-mass case, in D = 2 - 2*eps dimensions.
# Projective chart x3 = 1.

A = OreAlg(
    order = "lex dt > grevlex x y > grevlex dx dy",
    ratvars = ["eps", "m1", "m2", "m3"],
    ratdiffvars = (["t"], ["dt"]),
    poldiffvars = (["x", "y"], ["dx", "dy"]),
)

U = "(x*y + x + y)"
F = "(((m1^2)*x + (m2^2)*y + m3^2)*($U) - t*x*y)"
integrand = "1/($F)*(($U)^3/(($F)^2))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
