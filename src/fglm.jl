
# assume that g is a reduced GB for ord(A_gl)
# implementation according to the article 
# Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering
@inline function _find_predecessor(
    m::M,
    normal_forms::Dict{M,OrePoly{T,M}},
    A::OreAlg
) where {T,M}
    for i in 1:A.nrdv
        if m[i] > 0
            d = makemon(i, A)
            mprime = m / d
            if haskey(normal_forms, mprime)
                return true, d, mprime
            end
        end
    end
    one_mon = makemon(-1, A)
    return false, one_mon, one_mon
end

@inline function _predecessor_action(
    normal_form::OrePoly{T,M},
    di::M,
    A_gl::OreAlg,
    ::Val{false}
) where {T,M}
    return mul(makepoly(one(ctx(A_gl)), di), normal_form, A_gl)
end

@inline function _predecessor_action(
    normal_form::OrePoly{T,M},
    di::M,
    A_gl::OreAlg,
    ::Val{true}
) where {T,M}
    idx = findfirst(!iszero, exp(di))
    di2 = makemon(idx + div(A_gl.nrdv, 2), A_gl)
    return add(
        mul(makepoly(one(ctx(A_gl)), di), normal_form, A_gl),
        mul(makepoly(one(ctx(A_gl)), di2), normal_form, A_gl),
        A_gl)

end

@inline function _normal_form_image(
    m::M,
    g::Vector{OrePoly{T,M}},
    normal_forms::Dict{M,OrePoly{T,M}},
    geob::GeoBucket{T,M},
    tmp_poly::ReuseOrePoly{T,M},
    one_mon::M,
    A_gl::OreAlg,
    ann_prod::Val
) where {T,M}
    if m == one_mon
        return div!(geob, tmp_poly, makepoly(one(ctx(A_gl)), m), g, A_gl)
    end
    _, di, mprime = _find_predecessor(m, normal_forms, A_gl)
    mult = _predecessor_action(normal_forms[mprime], di, A_gl, ann_prod)
    return div!(geob, tmp_poly, mult, g, A_gl)
end

function _ann_prod_coeff_diagonal_substitution(A_gl::OreAlg)
    n = div(A_gl.nrdv, 2)
    ratdiffvars = A_gl.inp.ratdiffvars[1]
    ratvars = A_gl.inp.ratvars
    sub = Vector{eltype_co(A_gl)}(undef, length(ratdiffvars) + length(ratvars))

    @inbounds for i in eachindex(ratdiffvars)
        name = i <= n ? ratdiffvars[i] : ratdiffvars[i - n]
        sub[i] = A_gl.ratvars[name]
    end
    off = length(ratdiffvars)
    @inbounds for i in eachindex(ratvars)
        sub[off + i] = A_gl.ratvars[ratvars[i]]
    end

    return sub
end

function _specialize_ann_prod_normal_form(
    p::OrePoly{T,M},
    A_gl::OreAlg,
    coeff_sub::Vector{K}
) where {T<:RatFun,K,M}
    len = length(p)
    q = undefOrePoly(len, A_gl)
    c = coeffs(q)
    m = mons(q)

    @inbounds for j in 1:len
        c[j] = convertn(evaluate(coeff(p, j), coeff_sub), ctx(A_gl))
        m[j] = mon(p, j)
    end
    normalize!(q, A_gl)
    return q
end

@inline _specialize_ann_prod_normal_form(p::OrePoly, ::OreAlg, ::Val{false}, ::Nothing) = p
@inline function _specialize_ann_prod_normal_form(p::OrePoly{T,M}, A_gl::OreAlg, ::Val{true}, coeff_sub::Vector{K}) where {T<:RatFun,K,M}
    return _specialize_ann_prod_normal_form(p, A_gl, coeff_sub)
end

function _reconstruct_relation(m :: M,
    sec_v::Vector{T},
    mbasis1::Vector{M},
    A_gl::OreAlg
) where {T,M}
    @assert length(sec_v) == length(mbasis1) + 1
    cpoly = [sec_v[end]]
    mpoly = [m]
    for i in length(mbasis1):-1:1
        mo = mbasis1[i]
        c = sec_v[i]
        if !iszero(c, ctx(A_gl))
            push!(cpoly, c)
            push!(mpoly, mo)
        end
    end
    res = OrePoly(cpoly, mpoly) 
    ctx(A_gl) isa RingCtx && primitive_part!(res, A_gl)
    return res
end

@inline function update_todo_fglm!(
    todo::SortedSet{M},
    m::M,
    A_gl::OreAlg,
    ::Val{false}
) where {M}
    for i in 1:A_gl.nrdv
        push!(todo, m * makemon(i, A_gl))
    end
    return nothing
end

@inline function update_todo_fglm!(
    todo::SortedSet{M},
    m::M,
    A_gl::OreAlg,
    ::Val{true}
) where {M}
    ub = div(A_gl.nrdv, 2)
    for i in 1:ub
        push!(todo, m * makemon(i, A_gl))
    end
    return nothing
end

function fglm(g :: Vector{OrePoly{T,M}}, 
              A_gl :: OreAlg, 
              A_elim :: OreAlg,
              ann_prod :: Bool = false) where {T,M}
    return fglm(g, A_gl, A_elim, Val(ann_prod))
end

function fglm(g :: Vector{OrePoly{T,M}}, 
              A_gl :: OreAlg, 
              A_elim :: OreAlg,
              ann_prod :: Val{B}) where {T,M,B}
    newbasis = OrePoly{T,M}[]
    one_mon = makemon(-1,A_gl)
    if any(divide(mon(p,1), one_mon, A_gl) for p in g)
        return g 
    end

    # Monomials to process in target order, it corresponds to ListOfNext
    todo = SortedSet{M}(order(A_elim), [one_mon])
    done = SortedSet{M}(order(A_elim))

    # Independent monomials in the quotient (staircase)
    staircase = M[]


    # Echelon form of independent normal forms wrt `g` (in A_gl)
    # it is possible to reconstruct MBasis[2] from it (cf the article)
    echelon = OrePoly{T,M}[]
    echelonvect = Vector{T}[]

    # mbasis1 corresponds to the monomials under the stair for the new order. 
    mbasis1 = M[]

    # # Leading monomials already found for the target basis, used to prune.
    # lm_res = SortedSet{M}(order(A_elim))

    geob = GeoBucket(one(A_gl))
    tmp_poly = ReuseOrePoly(1, A_gl)
    normal_forms = Dict{M,OrePoly{T,M}}()
    ring_ctx = ctx(A_gl) isa RingCtx
    unitary_echelon = !ring_ctx
    coeff_sub = B ? _ann_prod_coeff_diagonal_substitution(A_gl) : nothing

    while !isempty(todo) 
        m = first(todo)
        delete!(todo, m)
        m in done && continue
        push!(done,m)

        # If m is a multiple of a known LT for the new basis, it is already handled.
        if !any(divide(mm, m, A_elim) for mm in staircase)
            vec = _normal_form_image(m, g, normal_forms, geob, tmp_poly, one_mon, A_gl, ann_prod)
            vec = _specialize_ann_prod_normal_form(vec, A_gl, ann_prod, coeff_sub)

            if ring_ctx
                primitive_part!(vec, A_gl)
            end
            normal_forms[m] = vec
            tmp, sec_v = reduce_with_echelon_augmented!(echelon, vec, A_gl, echelonvect)
            if ring_ctx
                primitive_part!(rel, A_gl)
            end
            if iszero(tmp,A_gl) # there exists a linear relation
                rel = _reconstruct_relation(m, sec_v, mbasis1, A_gl)
                push!(newbasis, rel)
                push!(staircase,m)
            else 
                add_echelon!(echelon, tmp, A_gl, augmented = true, echelonvect = echelonvect, vect = sec_v, unitary = unitary_echelon)
                push!(mbasis1,m)
                update_todo_fglm!(todo, m, A_gl, ann_prod)
            end
        end
    end
    return newbasis
end
