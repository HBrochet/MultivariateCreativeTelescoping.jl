"""
    ann_comp_right_rat(fname::Symbol, g, A::OreAlg; params = Dict{Symbol,Any}())

Closure by right composition for a database function symbol:

- reads the univariate annihilator template from `datab_LDE[fname]`;
- in closure algebra `A = dfinite_closure_A(A0)`, uses only the first half of
  ratdiff variables and ignores the duplicated second half;
- substitutes `x -> g` and `dx -> dx/g` (for the first active differential variable).
"""

function _dfinite_expr_to_ratfun(g, A::OreAlg)
    gp = parse_OrePoly(string(g), A)
    length(gp) == 1 || error("Expected a rational expression, got $(g)")
    mon(gp, 1) == makemon(-1, A) || error("Expected a rational expression, got $(g)")
    return coeff(gp, 1)
end

function ann_comp_right_rat(
    fname::Symbol,
    g :: Expr,
    A::OreAlg{T,C,MA,O};
    params::Dict{Symbol,Any} = Dict{Symbol,Any}()) where {T,C<:RatFunCtx,MA,O}

    lexpr = datab_LDE[fname]
    for (k, v) in params
        lexpr = _subs_expr(lexpr, Dict{Symbol,Any}(k => v))
    end
    nactive = div(A.nrdv, 2)
    g_rat = _dfinite_expr_to_ratfun(g, A)

    res = OrePoly{T,MA}[]

    for j in 1:nactive 
        xname = A.inp.ratdiffvars[1][j]
        dname = A.inp.ratdiffvars[2][j]
        idx = A.strvar_to_indexp[dname]

        # evaluate _x at g and _D at D[j]
        tmp = ann_comp_right_rat_ev_LDE(A, g, j, lexpr)

        # differentiate g 
        dg = ctx(A) isa MRatFunCtx ? derivative(g_rat, ctx(A).vars[A.drvars_to_int[xname]]) : derivative(g_rat)
        
        # evaluate D[j] at D[j]/dg or return D[j] is dg = 0 
        tmp = ann_comp_right_rat_ev_D(A, dg, idx, tmp)
        push!(res,tmp)
    end
    return res
end

function ann_comp_right_rat_ev_LDE(A::OreAlg, g::Expr, j::Int, L::Expr)
    xname = A.inp.ratdiffvars[1][j]
    dname = A.inp.ratdiffvars[2][j]
    idx = A.strvar_to_indexp[dname]
    subs = Dict{Symbol,Any}(:_x => g, :_D => Symbol(dname))
    lexpr = _subs_expr(L, subs)
    return parse_OrePoly(string(lexpr), A)
end
    
# evaluate D[j] at D[j]/g 
function ann_comp_right_rat_ev_D(A::OreAlg{T,C,MA,O}, dg::T, idx::Int, L::OrePoly{T,MA}) where {T,C,MA,O}
    if iszero(dg)
        return makepoly(one(ctx(A)), makemon(idx,A))
    end
    Dog = makepoly(inv(dg),makemon(idx,A)) # D over g
    powpoly = one(A)
    pow = 0
    res = zero(A)
    for i in length(L):-1:1
        targetpow = mon(L, i)[idx]
        while pow < targetpow
            powpoly = mul(powpoly, Dog, A)
            pow += 1
        end
        res = add!(res, mul(coeff(L, i), powpoly, A), A)
    end
    return res 
end
