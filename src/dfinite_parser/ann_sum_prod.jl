### The functions in this file assume that the ring contains only ratdiffvars and ratvars variables, with one polynomial variable used for closure.
### the ideals considered shouldn't involve this polynomial variable and should be D-finite

@inline function _shift_diff_block_d1_to_d2(m::OreMonVE{N,E}, n::Int) where {N,E}
    return OreMonVE(SVector{N,E}(i <= n ? E(0) : (i <= 2*n ? m[i] + m[i-n] : m[i]) for i in 1:N))
end

function _coeff_substitution_x1_to_x2(A::OreAlg, n::Int)
    ratdiffvars = A.inp.ratdiffvars[1]
    ratvars = A.inp.ratvars
    sub = Vector{eltype_co(A)}(undef, length(ratdiffvars) + length(ratvars))

    @inbounds for i in eachindex(ratdiffvars)
        name = i <= n ? ratdiffvars[i+n] : ratdiffvars[i]
        sub[i] = A.ratvars[name]
    end
    off = length(ratdiffvars)
    @inbounds for i in eachindex(ratvars)
        sub[off + i] = A.ratvars[ratvars[i]]
    end

    return sub
end

function _rename_x1_d1_to_x2_d2(
    p::OrePoly{T,M},
    A::OreAlg,
    n::Int,
    coeff_sub::Vector{K}
) where {T<:RatFun,K,M<:OreMonVE}
    len = length(p)
    q = undefOrePoly(len, A)
    c = coeffs(q)
    m = mons(q)

    @inbounds for j in 1:len
        c[j] = convertn(evaluate(coeff(p, j), coeff_sub), ctx(A))
        m[j] = _shift_diff_block_d1_to_d2(mon(p, j), n)
    end
    normalize!(q, A)
    return q
end

function ann_sum(g1 :: Vector{OrePoly{T,M}}, g2 :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T,M}
    _check_ann_sum_product_assumptions(A)

    # form the new ideal ug1 + (1-u)g2 with the pol. variable
    u = makemon(A.nrdv + 2*A.npdv + 1, A)
    u_poly = makepoly(one(ctx(A)), u)
    omu = sub(one(A), u_poly, A) # equals 1-u
    gens = OrePoly{T,M}[]
    sizehint!(gens, length(g1) + length(g2))

    for p in g1
        push!(gens, mul(u_poly, p, A))
    end
    for p in g2
        push!(gens, mul(omu, p, A))
    end

    # now compute a gb
    gb = f4(gens, A)

    # keep only operators without u
    u_idx = A.nrdv + 1
    filter!(p -> lm(p)[u_idx] == 0, gb)

    return gb
end

function ann_product(g1 :: Vector{OrePoly{T,M}}, g2 :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T,M}
    _check_ann_sum_product_assumptions(A)

    n = div(A.nrdv, 2)
    gens = OrePoly{T,M}[]
    coeff_sub = _coeff_substitution_x1_to_x2(A, n)

    append!(gens, g1)
    for p in g2
        push!(gens, _rename_x1_d1_to_x2_d2(p, A, n, coeff_sub))
    end

    gb = fglm(gens, A, A, true)

    return gb
end

function _check_ann_sum_product_assumptions(A :: OreAlg)
    if !isempty(A.inp.poldiffvars[1]) || !(length(A.inp.polvars) == 1) || !isempty(A.inp.locvars[1])
        error("ann_sum and ann_product currently require algebras with only ratdiffvars and ratvars")
    end
end
