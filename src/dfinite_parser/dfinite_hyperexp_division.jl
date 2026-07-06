function _dfinite_divide_by_hyperexp_ann_main(num::Vector{OrePoly{T,M}}, den, A::OreAlg) where {T,M}
    den_ann = _dfinite_expr_to_ann(den, A)
    R = _dfinite_hyperexp_logder_from_ann(den_ann, A)

    nactive = div(A.nrdv, 2)
    shifts = Vector{OrePoly{T,M}}(undef, nactive)
    shift_pos = zeros(Int, nvars(A))

    @inbounds for j in 1:nactive
        dname = A.inp.ratdiffvars[2][j]
        idx = A.strvar_to_indexp[dname]
        shift_pos[idx] = j
        shifts[j] = add!(makepoly(one(ctx(A)), makemon(idx, A)), makepoly(R[j], makemon(-1, A)), A)
    end

    pow_caches = [OrePoly{T,M}[one(A), shifts[j]] for j in 1:nactive]
    out = Vector{OrePoly{T,M}}(undef, length(num))
    @inbounds for i in eachindex(num)
        out[i] = _dfinite_hyperexp_shift_operator(num[i], shifts, shift_pos, pow_caches, A)
    end
    return out
end

function _dfinite_hyperexp_logder_from_ann(gens::Vector{<:OrePoly{T,M}}, A::OreAlg) where {T,M}
    nactive = div(A.nrdv, 2)
    R = Vector{T}(undef, nactive)
    found = falses(nactive)

    @inbounds for j in 1:nactive
        idx = A.strvar_to_indexp[A.inp.ratdiffvars[2][j]]
        for g in gens
            ok, a, b = _dfinite_hyperexp_first_order_parts(g, idx, A)
            ok || continue
            if !iszero(a, ctx(A))
                R[j] = -b / a
                found[j] = true
                break
            end
        end
    end

    all(found) || error("Division is currently supported only by hyperexponential denominators")
    return R
end

@inline function _dfinite_hyperexp_shift_pow!(
    shift::OrePoly{T,M},
    e::Int,
    cache::Vector{OrePoly{T,M}},
    A::OreAlg
) where {T,M}
    while length(cache) <= e + 1
        push!(cache, mul(cache[end], shift, A))
    end
    return cache[e + 1]
end

function _dfinite_hyperexp_shift_operator(
    L::OrePoly{T,M},
    shifts::Vector{OrePoly{T,M}},
    shift_pos::Vector{Int},
    pow_caches::Vector{Vector{OrePoly{T,M}}},
    A::OreAlg
) where {T,M}
    res = zero(A)
    nv = nvars(A)
    @inbounds for k in 1:length(L)
        m = mon(L, k)
        term = makepoly(coeff(L, k), makemon(-1, A))
        for i in 1:nv
            ei = m[i]
            ei == 0 && continue
            p = shift_pos[i]
            if p == 0
                term = mul(term, makepoly(one(ctx(A)), makemon(i, A)^ei), A)
            else
                factor = _dfinite_hyperexp_shift_pow!(shifts[p], Int(ei), pow_caches[p], A)
                term = mul(term, factor, A)
            end
        end
        res = add!(res, term, A)
    end
    return res
end

function _dfinite_hyperexp_first_order_parts(
    p::OrePoly{T,M},
    idx::Int,
    A::OreAlg
) where {T,M}
    one_mon = makemon(-1, A)
    d_mon = makemon(idx, A)
    a = zero(ctx(A))
    b = zero(ctx(A))

    @inbounds for k in 1:length(p)
        m = mon(p, k)
        c = coeff(p, k)
        if m == one_mon
            b = add(b, c, ctx(A))
        elseif m == d_mon
            a = add(a, c, ctx(A))
        else
            return false, a, b
        end
    end

    return true, a, b
end
