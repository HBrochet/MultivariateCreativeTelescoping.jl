
function _ann_poly_power_data(f_in::Union{String,Expr,Symbol}, s_in::Union{String,Symbol}, A::OreAlg)
    f_expr = f_in isa String ? Meta.parse(f_in) : f_in
    f_rat = _dfinite_expr_to_ratfun(f_expr, A)
    den = Nemo.denominator(f_rat, false)
    isone(den) || error("ann_poly_power expects a polynomial expression")

    sname = s_in isa String ? s_in : String(s_in)
    sname in A.inp.ratvars || error("ann_poly_power expects $(sname) to be declared in ratvars")
    scoeff = A.ratvars[sname]

    return f_rat, scoeff
end

function _ann_ratfun_derivative(r, xname::String, A::OreAlg)
    if ctx(A) isa MRatFunCtx
        return derivative(r, ctx(A).vars[A.drvars_to_int[xname]])
    end
    return derivative(r)
end

function ann_multiply_by_poly_power(num::Vector{OrePoly{T,M}},
                                    f_in::Union{String,Expr,Symbol},
                                    s_in::Union{String,Symbol},
                                    A::OreAlg{T,C,MA,O}) where {T,C<:RatFunCtx,M,MA,O}
    f_rat, scoeff = _ann_poly_power_data(f_in, s_in, A)
    nactive = div(A.nrdv, 2)

    shifts = Vector{OrePoly{T,M}}(undef, nactive)
    shift_pos = zeros(Int, nvars(A))

    @inbounds for j in 1:nactive
        xname = A.inp.ratdiffvars[1][j]
        dname = A.inp.ratdiffvars[2][j]
        idx = A.strvar_to_indexp[dname]
        shift_pos[idx] = j
        Rj = scoeff * _ann_ratfun_derivative(f_rat, xname, A) / f_rat
        shifts[j] = sub!(makepoly(one(ctx(A)), makemon(idx, A)), makepoly(Rj, makemon(-1, A)), A)
    end

    pow_caches = [OrePoly{T,M}[one(A), shifts[j]] for j in 1:nactive]
    out = Vector{OrePoly{T,M}}(undef, length(num))
    @inbounds for i in eachindex(num)
        out[i] = _dfinite_hyperexp_shift_operator(num[i], shifts, shift_pos, pow_caches, A)
    end
    return out
end

function ann_poly_power(f_in::Union{String,Expr,Symbol},
                        s_in::Union{String,Symbol},
                        A::OreAlg)
    nactive = div(A.nrdv, 2)
    gens = [parse_OrePoly(A.inp.ratdiffvars[2][i], A) for i in 1:nactive]
    return ann_multiply_by_poly_power(gens, f_in, s_in, A)
end
