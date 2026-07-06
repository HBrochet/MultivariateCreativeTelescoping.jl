function _ann_comp_right_alg_poly_xt_to_U(poly, xidx::Int, base_vars, K, U, t)
    out = zero(U)
    lp = length(poly)
    @inbounds for i in 1:lp
        ev = Nemo.exponent_vector(poly, i)
        td = ev[xidx]
        c = K(Nemo.coeff(poly, i))
        for j in 1:length(base_vars)
            j == xidx && continue
            ej = ev[j]
            ej == 0 && continue
            c *= K(base_vars[j]^ej)
        end
        out += U(c) * t^td
    end
    return out
end

function _ann_comp_right_alg_coeff_to_vec(c, xidx::Int, base_vars, p, n::Int, K, U, t)
    num = Nemo.numerator(c, false)
    den = Nemo.denominator(c, false)

    numu = _ann_comp_right_alg_poly_xt_to_U(num, xidx, base_vars, K, U, t)
    denu = _ann_comp_right_alg_poly_xt_to_U(den, xidx, base_vars, K, U, t)

    numv = _dfinite_poly_to_vec(rem(numu, p), n, K)
    denv = _dfinite_poly_to_vec(rem(denu, p), n, K)
    invdenv = _dfinite_inv_vec_mod_p(denv, p, n, U, t, K)
    return _dfinite_mul_vec_mod_p(numv, invdenv, p, n, U, t, K)
end

function _ann_comp_right_alg_apply_D!(out_blocks, blocks, dy, recs, p, xidx::Int, n::Int, r::Int, K, U, t)
    @inbounds for i in 1:r
        out_blocks[i] = _dfinite_derivative_vec(blocks[i], dy, p, xidx, n, K, U, t)
    end

    @inbounds for i in 1:(r - 1)
        v = _dfinite_mul_vec_mod_p(blocks[i], dy, p, n, U, t, K)
        out_blocks[i + 1] = map(+, out_blocks[i + 1], v)
    end

    brdy = _dfinite_mul_vec_mod_p(blocks[r], dy, p, n, U, t, K)
    @inbounds for j in 1:r
        v = _dfinite_mul_vec_mod_p(brdy, recs[j], p, n, U, t, K)
        out_blocks[j] = map(+, out_blocks[j], v)
    end

    return out_blocks
end

function _ann_comp_right_alg_blocks_to_col(blocks, n::Int, r::Int)
    col = Vector{eltype(blocks[1])}(undef, n * r)
    @inbounds for i in 1:r
        off = (i - 1) * n
        for j in 1:n
            col[off + j] = blocks[i][j]
        end
    end
    return col
end

function _ann_comp_right_alg_zero_blocks(n::Int, r::Int, K)
    return [[zero(K) for _ in 1:n] for _ in 1:r]
end

function ann_comp_right_alg(A::OreAlg{T,C,MA,O}, L::OrePoly{T,MA}, e::Expr) where {T,C<:RatFunCtx,MA,O}
    nactive = div(A.nrdv, 2)
    nactive == 0 && (nactive = A.nrdv)
    nactive > 0 || error("ann_comp_right_alg requires at least one differential variable")

    pmin = minimal_polynomial(e, A)
    pxy, base_syms = _dfinite_project_minpoly_to_base_tmp(pmin, A)

    dep_idxs = Int[]
    dep_dnames = String[]
    out_gens = OrePoly{T,MA}[]
    @inbounds for j in 1:nactive
        xnamej = Symbol(A.inp.ratdiffvars[1][j])
        dnamej = A.inp.ratdiffvars[2][j]
        if !_dfinite_expr_depends_on_symbol(e, xnamej)
            push!(out_gens, parse_OrePoly(dnamej, A))
            continue
        end
        push!(dep_idxs, j)
        push!(dep_dnames, dnamej)
    end

    ndep = length(dep_idxs)
    ndep == 0 && return out_gens

    # L is expected to be univariate in the first active differential variable.
    dmain = A.inp.ratdiffvars[2][1]
    didx = A.strvar_to_indexp[dmain]

    xname = Symbol(A.inp.ratdiffvars[1][1])
    xidx = findfirst(==(xname), base_syms)
    xidx === nothing && error("variable $(xname) not found in minimal polynomial base symbols")

    up, S0, K, U, t = _dfinite_minpoly_to_univariate_over_frac(pxy, base_syms)
    n = Nemo.degree(up)
    n > 0 || error("minimal polynomial has degree 0 in __tmp")

    lc = Nemo.coeff(up, n)
    p = lc == one(K) ? up : Nemo.divexact(up, lc)

    pt = derivative(p)
    ptv = _dfinite_poly_to_vec(rem(pt, p), n, K)
    invptv = _dfinite_inv_vec_mod_p(ptv, p, n, U, t, K)

    dy_all = Vector{Vector{typeof(one(K))}}(undef, ndep)
    @inbounds for jj in 1:ndep
        xj = dep_idxs[jj]
        px = zero(U)
        for i in 0:n
            ci = Nemo.coeff(p, i)
            iszero(ci) && continue
            px += U(derivative(ci, xj)) * t^i
        end
        pxv = _dfinite_poly_to_vec(rem(px, p), n, K)
        dy_all[jj] = _dfinite_mul_vec_mod_p(invptv, map(-, pxv), p, n, U, t, K)
    end

    maxpow = 0
    @inbounds for i in 1:length(L)
        m = mon(L, i)
        for j in 1:A.npdv
            j == didx && continue
            iszero(m[j]) || error("ann_comp_right_alg expects a univariate operator in $(dmain)")
        end
        pj = m[didx]
        pj > maxpow && (maxpow = pj)
    end
    r = Int(maxpow)
    r > 0 || error("ann_comp_right_alg expects a positive-order operator")

    az = zero(ctx(A))
    a = fill(az, r + 1)
    @inbounds for i in 1:length(L)
        k = mon(L, i)[didx]
        a[k + 1] = a[k + 1] + coeff(L, i)
    end

    base_vars = gens(S0)
    avec = Vector{Vector{typeof(one(K))}}(undef, r + 1)
    @inbounds for k in 0:r
        avec[k + 1] = _ann_comp_right_alg_coeff_to_vec(a[k + 1], xidx, base_vars, p, n, K, U, t)
    end

    ar = avec[r + 1]
    ar_inv = _dfinite_inv_vec_mod_p(ar, p, n, U, t, K)
    recs = Vector{Vector{typeof(one(K))}}(undef, r)
    @inbounds for j in 1:r
        recs[j] = _dfinite_mul_vec_mod_p(map(-, avec[j]), ar_inv, p, n, U, t, K)
    end

    Ndim = n * r

    blocks0 = _ann_comp_right_alg_zero_blocks(n, r, K)
    blocks0[1][1] = one(K)

    mon_basis = Vector{Vector{Int}}()
    vec_basis = Vector{Vector{typeof(one(K))}}()
    leading_mons = Vector{Vector{Int}}()
    normal_forms = Dict{Vector{Int},Vector{Vector{typeof(one(K))}}}()

    m0 = zeros(Int, ndep)
    normal_forms[m0] = blocks0
    todo = Vector{Vector{Int}}(undef, 1)
    todo[1] = m0
    todo_set = Set{Vector{Int}}()
    push!(todo_set, m0)
    head = 1

    while head <= length(todo)
        m = todo[head]
        head += 1
        delete!(todo_set, m)

        skip_m = false
        @inbounds for lm in leading_mons
            if _dfinite_is_multiple_of_monvec(m, lm)
                skip_m = true
                break
            end
        end
        skip_m && continue

        b = get(normal_forms, m, nothing)
        if b === nothing
            found = false
            @inbounds for i in 1:ndep
                m[i] == 0 && continue
                mp = copy(m)
                mp[i] -= 1
                bp = get(normal_forms, mp, nothing)
                bp === nothing && continue
                bo = _ann_comp_right_alg_zero_blocks(n, r, K)
                _ann_comp_right_alg_apply_D!(bo, bp, dy_all[i], recs, p, dep_idxs[i], n, r, K, U, t)
                normal_forms[m] = bo
                b = bo
                found = true
                break
            end
            found || error("failed to construct derivative image in ann_comp_right_alg")
        end

        v = _ann_comp_right_alg_blocks_to_col(b, n, r)
        k = length(vec_basis)
        if k == 0
            push!(vec_basis, v)
            push!(mon_basis, m)
            @inbounds for i in 1:ndep
                mc = copy(m)
                mc[i] += 1
                if !haskey(normal_forms, mc) && !(mc in todo_set)
                    push!(todo, mc)
                    push!(todo_set, mc)
                end
            end
            continue
        end

        M = matrix_space(K, Ndim, k + 1)()
        @inbounds for j in 1:k, i in 1:Ndim
            M[i, j] = vec_basis[j][i]
        end
        @inbounds for i in 1:Ndim
            M[i, k + 1] = v[i]
        end

        null_dim, N = nullspace(M)
        if null_dim == 0
            push!(vec_basis, v)
            push!(mon_basis, m)
            @inbounds for i in 1:ndep
                mc = copy(m)
                mc[i] += 1
                if !haskey(normal_forms, mc) && !(mc in todo_set)
                    push!(todo, mc)
                    push!(todo_set, mc)
                end
            end
            continue
        end

        rel_col = 0
        @inbounds for c in 1:null_dim
            if !iszero(N[k + 1, c])
                rel_col = c
                break
            end
        end
        if rel_col == 0
            push!(vec_basis, v)
            push!(mon_basis, m)
            @inbounds for i in 1:ndep
                mc = copy(m)
                mc[i] += 1
                if !haskey(normal_forms, mc) && !(mc in todo_set)
                    push!(todo, mc)
                    push!(todo_set, mc)
                end
            end
            continue
        end

        leadc = N[k + 1, rel_col]
        terms = String[]
        tm = _dfinite_term_string(N[k + 1, rel_col] / leadc, m, dep_dnames, ndep)
        !isempty(tm) && push!(terms, tm)
        @inbounds for i in 1:k
            ti = _dfinite_term_string(N[i, rel_col] / leadc, mon_basis[i], dep_dnames, ndep)
            !isempty(ti) && push!(terms, ti)
        end

        isempty(terms) && error("computed zero annihilator in ann_comp_right_alg")
        push!(out_gens, parse_OrePoly(join(terms, " + "), A))
        push!(leading_mons, m)
    end

    return out_gens
end
