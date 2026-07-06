include(joinpath(@__DIR__, "_common.jl"))

# Section 4.1 of arXiv:2401.09908
# Massless box graph in D = 4 - 2*eps dimensions.
# Projective chart x4 = 1 and single-scale variable X = t / s.

A = OreAlg(
    order = "lex dX > grevlex x1 x2 x3 > grevlex dx1 dx2 dx3",
    ratvars = ["eps"],
    ratdiffvars = (["X"], ["dX"]),
    poldiffvars = (["x1", "x2", "x3"], ["dx1", "dx2", "dx3"]),
)

U = "(x1 + x2 + x3 + 1)"
F = "(x1*x3 + X*x2)"
integrand = "1/(($F)^2)*((($U)^2)/($F))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
