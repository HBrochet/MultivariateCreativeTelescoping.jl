include(joinpath(@__DIR__, "_common.jl"))

# Section 5.4 of arXiv:2401.09908
# Three-point non-planar triangle-box graph in D = 4 - 2*eps dimensions.
# Projective chart x6 = 1.

A = OreAlg(
    order = "lex dX > grevlex x1 x2 x3 x4 x5 > grevlex dx1 dx2 dx3 dx4 dx5",
    ratvars = ["eps"],
    ratdiffvars = (["X"], ["dX"]),
    poldiffvars = (["x1", "x2", "x3", "x4", "x5"], ["dx1", "dx2", "dx3", "dx4", "dx5"]),
)

U = "((x1 + x2)*(x3 + x4 + x5 + 1) + (x3 + x4)*(x5 + 1))"
F = "(-(((x3 + x4 + x5 + 1)*x1*x2 + x1*x3*x5 + x2*x4)*X) + (x3 + x4 + x5 + 1)*($U))"
integrand = "1/(($F)^2)*((($U)^3)/(($F)^2))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
