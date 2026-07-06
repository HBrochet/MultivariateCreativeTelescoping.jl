include(joinpath(@__DIR__, "_common.jl"))

# Section 4.2 of arXiv:2401.09908
# Cross Witten diagram in AdS4.
# Projective chart x3 = 1.

A = OreAlg(
    order = "lex dz dzb > grevlex x1 x2 > grevlex dx1 dx2",
    ratvars = ["eps"],
    ratdiffvars = (["z", "zb"], ["dz", "dzb"]),
    poldiffvars = (["x1", "x2"], ["dx1", "dx2"]),
)

U = "(x1 + x2 + 1)"
F = "(x1*x2 + z*zb*x1 + (1-z)*(1-zb)*x2)"
integrand = "1/(($U)*($F))*((($F)^2)/(x2^4))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
