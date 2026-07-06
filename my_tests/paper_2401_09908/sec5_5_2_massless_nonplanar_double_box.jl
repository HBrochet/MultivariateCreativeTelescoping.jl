include(joinpath(@__DIR__, "_common.jl"))

# Section 5.5.2 of arXiv:2401.09908
# Massless non-planar double-box graph in D = 4 - 2*eps dimensions.
# Projective chart x7 = 1 and X = t / s.

A = OreAlg(
    order = "lex dX > grevlex x1 x2 x3 x4 x5 x6 > grevlex dx1 dx2 dx3 dx4 dx5 dx6",
    ratvars = ["eps"],
    ratdiffvars = (["X"], ["dX"]),
    poldiffvars = (["x1", "x2", "x3", "x4", "x5", "x6"], ["dx1", "dx2", "dx3", "dx4", "dx5", "dx6"]),
)

U = "((x1 + x3 + x4)*(x2 + x5 + x6 + 1) + (x2 + 1)*(x5 + x6))"
F = "(x1*x3*(x2 + x5 + x6 + 1) + x1*x6 + x2*(x3*x5 - x4*x6) + X*x4*(x5 - x2*x6))"
integrand = "($U)/(($F)^3)*((($U)^3)/(($F)^2))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
