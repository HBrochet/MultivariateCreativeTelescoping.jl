include(joinpath(@__DIR__, "_common.jl"))

# Section 5.3 of arXiv:2401.09908
# Two-loop ice-cream cone, equal-kinematics case.
# Projective chart z = 1.

A = OreAlg(
    order = "lex dt > grevlex y1 y2 x1 > grevlex dy1 dy2 dx1",
    ratvars = ["eps"],
    ratdiffvars = (["t"], ["dt"]),
    poldiffvars = (["y1", "y2", "x1"], ["dy1", "dy2", "dx1"]),
)

U = "((y1 + y2)*(x1 + 1) + x1)"
V = "(y1*y2*(1 + x1) + x1*(y1 + y2))"
F = "((y1 + y2 + x1 + 1)*($U) - t*($V))"
integrand = "($U)/(($F)^2)*((($U)^3)/(($F)^2))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
