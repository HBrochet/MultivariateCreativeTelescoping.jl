include(joinpath(@__DIR__, "_common.jl"))

# Section 5.6 of arXiv:2401.09908
# Witten ice-cream cone diagram in analytic regularisation.
# Projective chart x4 = 1.

A = OreAlg(
    order = "lex du dv > grevlex x1 x2 x3 > grevlex dx1 dx2 dx3",
    ratvars = ["kappa"],
    ratdiffvars = (["u", "v"], ["du", "dv"]),
    poldiffvars = (["x1", "x2", "x3"], ["dx1", "dx2", "dx3"]),
)

U = "(x1*x2 + (x1 + x2)*(x3 + 1))"
F = "(u*x1*x2*x3 + v*x1*x2 + (x1 + x2)*x3)"
integrand = "1/(($U)^2)*((x1*x2*x3*($U))/($F))^kappa"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
