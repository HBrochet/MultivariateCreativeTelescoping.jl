"""
    weyl_closure_internal_single_T_nomul_di(gens, A, init; param = f4_param())

Variant of `weyl_closure_internal` for experimentation:

- compute the singular locus polynomial `p`;
- add one localization variable `_T` for `1/p`;
- add the single localization relation `p*_T - 1`;
- add only one premultiple `_T * g` for each generator `g`;
- forbid Grobner reductions that change the exponents of differential variables
  by adding all differential-variable indices to `nomul`;
- run one `f4` computation only, without the iterative saturation loop.

Returns a named tuple with:

- `gb_full`: the full Grobner basis in the localized algebra;
- `gb_T0`: the subfamily whose leading monomial has `_T`-degree `0`;
- `localized_algebra`: the localized algebra used for the computation;
- `singular_locus`: the square-free singular-locus polynomial used for `_T`.
"""
function weyl_closure_internal_single_T_nomul_di(
    gens::Vector{OrePoly{T,M}},
    A::OreAlg,
    init::WeylClosureInit;
    param = f4_param(),
) where {T,M}
    sl = singular_locus(gens, A, init.SLI)

    if sl == one(A)
        gb = f4(gens, A; param = param)
        return (gb_full = gb,
                gb_T0 = gb,
                localized_algebra = A,
                singular_locus = sl)
    end

    nemo_p = prod(fact[1] for fact in factor(sl))
    p = OrePoly([nemo_p], [makemon(-1, A)])

    nA = add_loc_var("_T", "1", A)
    (polsloc, diff_pols_loc) = compute_locpol_and_derivs([mystring(p, A)], nA)
    empty!(nA.pols_loc)
    append!(nA.pols_loc, polsloc)
    empty!(nA.diff_pols_loc)
    append!(nA.diff_pols_loc, diff_pols_loc)

    ngens = map_algebras(gens, A, nA)

    Tidx = nA.strvar_to_indexp["_T"]
    Tpoly = makepoly(one(ctx(nA)), makemon(Tidx, nA))

    locrel = sub(mul(nA.pols_loc[1], Tpoly, nA), one(nA), nA)
    work = copy(ngens)
    append!(work, [mul(Tpoly, g, nA) for g in ngens])
    push!(work, locrel)

    dnames = String[]
    append!(dnames, nA.inp.ratdiffvars[2])
    append!(dnames, nA.inp.poldiffvars[2])
    didxs = [nA.strvar_to_indexp[dname] for dname in dnames]

    dvars = [makemon(didx, nA) for didx in didxs]
    dmons = typeof(makemon(-1, nA))[]
    expT = eltype(mon(locrel, 1).exp)
    seen_dexps = Set{typeof(ntuple(_ -> zero(expT), length(didxs)))}()
    for g in work
        for i in 1:length(g)
            m = mon(g, i)
            dexps = ntuple(j -> m[didxs[j]], length(didxs))
            all(iszero, dexps) && continue
            dexps in seen_dexps && continue
            push!(seen_dexps, dexps)
            dm = makemon(-1, nA)
            for j in 1:length(didxs)
                e = dexps[j]
                e == 0 && continue
                dm *= dvars[j]^e
            end
            push!(dmons, dm)
        end
    end
    append!(work, [mul(dm, locrel, nA) for dm in dmons])

    maxT = maxdeg(work, Tidx)
    for k in 1:maxT
        push!(work, mul(makemon(Tidx, nA)^k, locrel, nA))
    end

    for didx in didxs
        didx in nA.nomul || push!(nA.nomul, didx)
    end

    gb = f4(work, nA; param = param)
    gb_T0 = filter(g -> mon(g, 1)[Tidx] == 0, gb)

    return (gb_full = gb,
            gb_T0 = gb_T0,
            localized_algebra = nA,
            singular_locus = p)
end

"""
    weyl_closure_internal_single_T_nomul_di2(gens, A, init; param = f4_param())

Variant of `weyl_closure_internal_single_T_nomul_di` where, before adding the
single premultiples by `_T`, we enlarge the generator family by all left
products `d_i * g` with `g in gens` and `d_i` ranging over the differential
variables of the algebra.
"""
function weyl_closure_internal_single_T_nomul_di2(
    gens::Vector{OrePoly{T,M}},
    A::OreAlg,
    init::WeylClosureInit;
    param = f4_param(),
) where {T,M}
    sl = singular_locus(gens, A, init.SLI)

    if sl == one(A)
        gb = f4(gens, A; param = param)
        return (gb_full = gb,
                gb_T0 = gb,
                localized_algebra = A,
                singular_locus = sl)
    end

    nemo_p = prod(fact[1] for fact in factor(sl))
    p = OrePoly([nemo_p], [makemon(-1, A)])

    nA = add_loc_var("_T", "1", A)
    (polsloc, diff_pols_loc) = compute_locpol_and_derivs([mystring(p, A)], nA)
    empty!(nA.pols_loc)
    append!(nA.pols_loc, polsloc)
    empty!(nA.diff_pols_loc)
    append!(nA.diff_pols_loc, diff_pols_loc)

    ngens = map_algebras(gens, A, nA)

    Tidx = nA.strvar_to_indexp["_T"]
    Tpoly = makepoly(one(ctx(nA)), makemon(Tidx, nA))

    dnames = String[]
    append!(dnames, nA.inp.ratdiffvars[2])
    append!(dnames, nA.inp.poldiffvars[2])
    didxs = [nA.strvar_to_indexp[dname] for dname in dnames]

    gens_di = copy(ngens)
    for didx in didxs
        dpoly = makepoly(one(ctx(nA)), makemon(didx, nA))
        append!(gens_di, [mul(dpoly, g, nA) for g in ngens])
    end

    locrel = sub(mul(nA.pols_loc[1], Tpoly, nA), one(nA), nA)
    work = copy(gens_di)
    append!(work, [mul(Tpoly, g, nA) for g in gens_di])
    push!(work, locrel)

    dvars = [makemon(didx, nA) for didx in didxs]
    dmons = typeof(makemon(-1, nA))[]
    expT = eltype(mon(locrel, 1).exp)
    seen_dexps = Set{typeof(ntuple(_ -> zero(expT), length(didxs)))}()
    for g in gens_di
        for i in 1:length(g)
            m = mon(g, i)
            dexps = ntuple(j -> m[didxs[j]], length(didxs))
            all(iszero, dexps) && continue
            dexps in seen_dexps && continue
            push!(seen_dexps, dexps)
            dm = makemon(-1, nA)
            for j in 1:length(didxs)
                e = dexps[j]
                e == 0 && continue
                dm *= dvars[j]^e
            end
            push!(dmons, dm)
        end
    end
    append!(work, [mul(dm, locrel, nA) for dm in dmons])

    maxT = maxdeg(work, Tidx)
    Tvar = makemon(Tidx, nA)
    Tmon = Tvar
    sizehint!(work, length(work) + maxT * (length(dmons) + 1))
    for k in 1:maxT
        push!(work, mul(Tmon, locrel, nA))
        for dm in dmons
            push!(work, mul(dm * Tmon, locrel, nA))
        end
        Tmon *= Tvar
    end

    for didx in didxs
        didx in nA.nomul || push!(nA.nomul, didx)
    end

    gb = f4(work, nA; param = param)
    gb_T0 = filter(g -> mon(g, 1)[Tidx] == 0, gb)

    return (gb_full = gb,
            gb_T0 = gb_T0,
            localized_algebra = nA,
            singular_locus = p)
end
