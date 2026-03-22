function compute_next_rel_many(der_maps :: Vector{Dict{M,OrePoly{T,M}}}, pol :: OrePoly{T,M}, ind :: Int, A :: OreAlg{T,C,M,O}) where {T,C,M,O}
    res = diff(pol, ind, A)
    for (c,m) in pol
        res = add!(res, mul(c, der_maps[ind][m], A), A)
    end
    return res
end

function normalize_relation_many!(rel :: OrePoly{T,M}, A :: OreAlg{T,C,M,O}) where {T,C,M,O}
    normalize!(rel, A)
    clear_denominators!(rel, A)
    normalize!(rel, A)
    ctx(A) isa RingCtx && primitive_part!(rel, A)
    return rel
end

"""
    find_LDE_direct_many(der_maps :: Vector{Dict{M,OrePoly{T,M}}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}

Compute a mixed system of differential equations satisfied by the integral, by iterating over
monomials in `A.inp.ratdiffvars[2]` in increasing total degree and detecting linear dependencies
between their reduced images.
"""
function find_LDE_direct_many(der_maps :: Vector{Dict{M,OrePoly{T,M}}}, spol :: OrePoly{T,M}, A :: OreAlg{T,C,M,O}) where {T,C,M,O}
    A.nrdv == 0 && return OrePoly{T,M}[]

    one_mon = makemon(-1, A)
    todo = [SortedSet{M}(order(A), [one_mon]), SortedSet{M}(order(A))]
    done = SortedSet{M}(order(A))
    leading_mons = M[]
    normal_forms = Dict{M,OrePoly{T,M}}()
    newbasis = OrePoly{T,M}[]

    echelon = OrePoly{T,M}[]
    echelonvect = Vector{T}[]
    mbasis1 = M[]
    ring_ctx = ctx(A) isa RingCtx
    unitary_echelon = !ring_ctx
    current_deg = 0

    while !isempty(todo[current_deg+1])
        while !isempty(todo[current_deg+1])
            m = popfirst!(todo[current_deg+1])
            m in done && continue
            push!(done, m)

            any(mm -> divide(mm, m, A), leading_mons) && continue

            vec = if m == one_mon
                deepcopy(spol)
            else
                found, di, mprime = _find_predecessor(m, normal_forms, A)
                found || error("Unable to reconstruct the predecessor of a ratdiff monomial")
                compute_next_rel_many(der_maps, normal_forms[mprime], findfirst(!iszero, exp(di)), A)
            end
            normal_forms[m] = vec

            tmp, sec_v = reduce_with_echelon_augmented!(echelon, vec, A, echelonvect)
            if iszero(tmp, A)
                rel = _reconstruct_relation(m, sec_v, mbasis1, A)
                normalize_relation_many!(rel, A)
                push!(newbasis, rel)
                push!(leading_mons, m)
                continue
            end

            add_echelon!(echelon, tmp, A, augmented = true, echelonvect = echelonvect, vect = sec_v, unitary = unitary_echelon)
            push!(mbasis1, m)
            for i in 1:A.nrdv
                push!(todo[current_deg+2], m * makemon(i, A))
            end
        end

        isempty(todo[current_deg+2]) && return newbasis
        push!(todo, SortedSet{M}(order(A)))
        current_deg += 1
    end

    return newbasis
end

function find_LDE_direct_many(data :: NamedTuple, A :: OreAlg)
    return find_LDE_direct_many(data.der_maps, data.spol, A)
end

"""
    MCTMany(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg; param :: MCTParam = mct_param())

Compute a mixed system of differential equations directly over the current coefficient field.
"""
function MCTMany(spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}, A::OreAlg{T,C,M,O}; param :: MCTParam = mct_param()) where {T,C,M,O}
    push!(A.nomul, 1)
    try
        red = der_red_map_many(A, spol, gb, param)
        debug(param) && @debug "der_red_map_many computed directly, starting to compute the mixed LDE system"
        return find_LDE_direct_many(red, A)
    finally
        deleteat!(A.nomul, length(A.nomul))
    end
end

function MCTMany(A::OreAlg{T,C,M,O}, spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}; param :: MCTParam = mct_param()) where {T,C,M,O}
    return MCTMany(spol, gb, A; param = param)
end
