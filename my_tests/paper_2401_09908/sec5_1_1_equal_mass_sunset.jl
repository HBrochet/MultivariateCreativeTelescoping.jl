include(joinpath(@__DIR__, "_common.jl"))

# Section 5.1.1 of arXiv:2401.09908
# Two-loop equal-mass sunset in D = 2 - 2*eps dimensions.
# We deprojectivize the projective form on the affine chart x3 = 1.
#
# U = x*y + x + y
# F = (x + y + 1)*U - t*x*y
# integrand = 1/F * (U^3/F^2)^eps

A = OreAlg(
    order = "lex dt > grevlex x y > grevlex dx dy",
    ratvars = ["eps"],
    ratdiffvars = (["t"], ["dt"]),
    poldiffvars = (["x", "y"], ["dx", "dy"]),
)

U = "(x*y + x + y)"
F = "((x + y + 1)*($U) - t*x*y)"
integrand = "1/($F)*(($U)^3/(($F)^2))^eps"

ann = annihilator_from_dfinite_parser(integrand, A)

init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)

spol = one(A)
sys = MCTMany(spol, gb, A)
