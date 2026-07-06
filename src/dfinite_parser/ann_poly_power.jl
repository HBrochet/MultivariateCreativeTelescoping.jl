
function _ann_poly_power_exponent_coeff(s_in, A::OreAlg)
    if s_in isa String || s_in isa Symbol
        sname = s_in isa String ? s_in : String(s_in)
        sname in A.inp.ratvars || error("ann_poly_power expects $(sname) to be declared in ratvars")
        return A.ratvars[sname]
    end

    _dfinite_expr_is_ratvar_rational(s_in, A) || error("ann_poly_power expects an exponent that is a rational expression in ratvars only")
    return _dfinite_expr_to_ratfun(s_in, A)
end

function _ann_poly_power_data(f_in::Union{String,Expr,Symbol}, s_in, A::OreAlg)
    f_expr = f_in isa String ? Meta.parse(f_in) : f_in
    f_rat = _dfinite_expr_to_ratfun(f_expr, A)
    scoeff = _ann_poly_power_exponent_coeff(s_in, A)
    return f_rat, scoeff
end

function _ann_ratfun_derivative(r, xname::String, A::OreAlg)
    return coeff_derivative(r, A.drvars_to_int[xname], A)
end

function ann_multiply_by_poly_power(num::Vector{OrePoly{T,M}},
                                    f_in::Union{String,Expr,Symbol},
                                    s_in,
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
                        s_in,
                        A::OreAlg)
    nactive = div(A.nrdv, 2)
    gens = [parse_OrePoly(A.inp.ratdiffvars[2][i], A) for i in 1:nactive]
    return ann_multiply_by_poly_power(gens, f_in, s_in, A)
end
