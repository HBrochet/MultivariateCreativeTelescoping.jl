include(joinpath(@__DIR__, "_common.jl"))

# Section 6.1.2 of arXiv:2401.09908
# Three-loop sunset, mass configuration [211].
# Projective chart x4 = 1.

A = OreAlg(
    order = "lex dt > grevlex x1 x2 x3 > grevlex dx1 dx2 dx3",
    ratvars = ["eps", "ma", "mb", "mc"],
    ratdiffvars = (["t"], ["dt"]),
    poldiffvars = (["x1", "x2", "x3"], ["dx1", "dx2", "dx3"]),
)

U = "(x1*x2*x3 + x1*x2 + x1*x3 + x2*x3)"
F = "(($U)*((ma^2)*(x1 + x2) + (mb^2)*x3 + mc^2) - t*x1*x2*x3)"
integrand = "1/($F)*(($U)^4/(($F)^3))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
