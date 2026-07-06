include(joinpath(@__DIR__, "_common.jl"))

# Section 6.1.1 of arXiv:2401.09908
# Higher-loop equal-mass sunset family in D = 2 - 2*eps dimensions.
# The default below is n = 5, corresponding to a 4-loop sunset.

function _monomial_str(vars::Vector{String})
    isempty(vars) && return "1"
    return join(vars, "*")
end

function _equal_mass_higher_loop_sunset(n::Int)
    vars = ["x$i" for i in 1:n-1]
    dvars = ["dx$i" for i in 1:n-1]
    allvars = vcat(vars, ["1"])
    uterms = String[]
    for i in eachindex(allvars)
        others = [allvars[j] for j in eachindex(allvars) if j != i]
        push!(uterms, _monomial_str(others))
    end
    U = "(" * join(uterms, " + ") * ")"
    S = "(" * join(allvars, " + ") * ")"
    prodx = _monomial_str(vars)
    F = "(($U)*($S) - t*($prodx))"
    order = "lex dt > grevlex " * join(vars, " ") * " > grevlex " * join(dvars, " ")
    A = OreAlg(
        order = order,
        ratvars = ["eps"],
        ratdiffvars = (["t"], ["dt"]),
        poldiffvars = (vars, dvars),
    )
    integrand = "1/($F)*(($U)^$n/(($F)^$(n-1)))^eps"
    return A, integrand
end

n = 5
A, integrand = _equal_mass_higher_loop_sunset(n)
ann = annihilator_from_dfinite_parser(integrand, A)
init = weyl_closure_init(A)
gb = weyl_closure(ann, A, init)
sys = MCTMany(one(A), gb, A)
