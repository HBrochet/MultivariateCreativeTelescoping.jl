include(joinpath(@__DIR__, "_common.jl"))

# Section 5.2 of arXiv:2401.09908
# Two-point one-mass kite in D = 4 - 2*eps dimensions.
# Projective chart x5 = 1.

A = OreAlg(
    order = "lex dX > grevlex x1 x2 x3 x4 > grevlex dx1 dx2 dx3 dx4",
    ratvars = ["eps"],
    ratdiffvars = (["X"], ["dX"]),
    poldiffvars = (["x1", "x2", "x3", "x4"], ["dx1", "dx2", "dx3", "dx4"]),
)

U = "((x1 + x2)*(x3 + x4) + x1 + x2 + x3 + x4)"
F = "(X*(x2*x3*x4 + x1*x3*x4 + x1*x2*x4 + x1*x2*x3 + (x1 + x4)*(x2 + x3)) - (x1 + x3 + 1)*($U))"
integrand = "1/(($U)*($F))*((($U)^3)/(($F)^2))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
