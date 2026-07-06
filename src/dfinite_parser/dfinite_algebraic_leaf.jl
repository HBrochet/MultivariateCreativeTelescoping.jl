function _dfinite_expr_depends_on_symbol(x, s::Symbol)
    if x isa Symbol
        return x == s
    elseif x isa Number
        return false
    elseif x isa Expr
        @inbounds for a in x.args
            _dfinite_expr_depends_on_symbol(a, s) && return true
        end
    end
    return false
end

function _dfinite_project_minpoly_to_base_tmp(p, A::OreAlg)
    raw_syms = ctx(A).R.S
    base_syms = raw_syms isa Symbol ? Symbol[raw_syms] : Symbol.(raw_syms)
    nbase = length(base_syms)

    Rsrc = parent(p)
    src_vars = gens(Rsrc)
    nsrc = length(src_vars)
    nsrc >= nbase + 1 || error("minimal polynomial ring has unexpected number of variables")
    ysrc_idx = nsrc

    Rdst, dst_vars = polynomial_ring(base_ring(Rsrc), vcat(base_syms, [:__tmp]))

    vals0 = [zero(Rdst) for _ in 1:nsrc]
    vals1 = [zero(Rdst) for _ in 1:nsrc]
    @inbounds for j in 1:nbase
        vals0[j] = dst_vars[j]
        vals1[j] = dst_vars[j]
    end
    @inbounds for j in (nbase + 1):(ysrc_idx - 1)
        vals0[j] = zero(Rdst)
        vals1[j] = one(Rdst)
    end
    vals0[ysrc_idx] = dst_vars[end]
    vals1[ysrc_idx] = dst_vars[end]

    pdst = Nemo.evaluate(p, vals0)
    pdst_check = Nemo.evaluate(p, vals1)
    pdst == pdst_check || error("minimal polynomial still depends on auxiliary variables")

    return pdst, base_syms
end

function _dfinite_minpoly_to_univariate_over_frac(pxy, base_syms)
    Rxy = parent(pxy)
    vars = gens(Rxy)
    nbase = length(base_syms)
    nvars = length(vars)
    nvars == nbase + 1 || error("unexpected polynomial ring shape")

    S0, rem_vars = polynomial_ring(base_ring(Rxy), base_syms)
    K = fraction_field(S0)
    U, t = polynomial_ring(K, :__tmp)

    up = zero(U)
    @inbounds for i in 1:length(pxy)
        ev = Nemo.exponent_vector(pxy, i)
        d = ev[end]
        c = S0(Nemo.coeff(pxy, i))
        for j in 1:nbase
            ej = ev[j]
            ej == 0 && continue
            c *= rem_vars[j]^ej
        end
        up += K(c) * t^d
    end

    return up, S0, K, U, t
end

function _dfinite_vec_to_poly(v, U, t)
    out = zero(U)
    @inbounds for i in 1:length(v)
        iszero(v[i]) && continue
        out += U(v[i]) * t^(i - 1)
    end
    return out
end

function _dfinite_poly_to_vec(f, n::Int, K)
    v = [zero(K) for _ in 1:n]
    @inbounds for i in 0:(n - 1)
        v[i + 1] = Nemo.coeff(f, i)
    end
    return v
end

function _dfinite_mul_vec_mod_p(a, b, p, n::Int, U, t, K)
    ap = _dfinite_vec_to_poly(a, U, t)
    bp = _dfinite_vec_to_poly(b, U, t)
    rp = rem(ap * bp, p)
    return _dfinite_poly_to_vec(rp, n, K)
end

function _dfinite_inv_vec_mod_p(a, p, n::Int, U, t, K)
    M = matrix_space(K, n, n)()
    @inbounds for j in 1:n
        ej = [zero(K) for _ in 1:n]
        ej[j] = one(K)
        col = _dfinite_mul_vec_mod_p(a, ej, p, n, U, t, K)
        for i in 1:n
            M[i, j] = col[i]
        end
    end

    b = matrix_space(K, n, 1)()
    b[1, 1] = one(K)
    @inbounds for i in 2:n
        b[i, 1] = zero(K)
    end
    x = solve(M, b; side = :right)
    return [x[i, 1] for i in 1:n]
end

function _dfinite_derivative_vec(v, dy, p, xidx::Int, n::Int, K, U, t)
    vp = _dfinite_vec_to_poly(v, U, t)

    dcoeff = zero(U)
    @inbounds for i in 0:Nemo.degree(vp)
        ci = Nemo.coeff(vp, i)
        iszero(ci) && continue
        dcoeff += U(derivative(ci, xidx)) * t^i
    end

    dt = derivative(vp)
    dyp = _dfinite_vec_to_poly(dy, U, t)
    return _dfinite_poly_to_vec(rem(dcoeff + dt * dyp, p), n, K)
end

function _dfinite_nullspace_operator_from_minpoly(pxy, base_syms, xidx::Int, dname::String, A::OreAlg)
    up, _, K, U, t = _dfinite_minpoly_to_univariate_over_frac(pxy, base_syms)
    n = Nemo.degree(up)
    n > 0 || error("minimal polynomial has degree 0 in __tmp")

    lc = Nemo.coeff(up, n)
    p = lc == one(K) ? up : Nemo.divexact(up, lc)

    pt = derivative(p)
    px = zero(U)
    @inbounds for i in 0:n
        ci = Nemo.coeff(p, i)
        iszero(ci) && continue
        px += U(derivative(ci, xidx)) * t^i
    end

    yt = _dfinite_poly_to_vec(rem(t, p), n, K)
    ptv = _dfinite_poly_to_vec(rem(pt, p), n, K)
    pxv = _dfinite_poly_to_vec(rem(px, p), n, K)
    invptv = _dfinite_inv_vec_mod_p(ptv, p, n, U, t, K)
    dy = _dfinite_mul_vec_mod_p(invptv, map(-, pxv), p, n, U, t, K)

    cols = Vector{Vector{typeof(one(K))}}(undef, n + 1)
    cols[1] = yt
    @inbounds for k in 2:(n + 1)
        cols[k] = _dfinite_derivative_vec(cols[k - 1], dy, p, xidx, n, K, U, t)
    end

    M = matrix_space(K, n, n + 1)()
    @inbounds for j in 1:(n + 1), i in 1:n
        M[i, j] = cols[j][i]
    end

    null_dim, N = nullspace(M)
    null_dim > 0 || error("failed to find annihilator from minimal polynomial")

    coeffs = [N[i, 1] for i in 1:(n + 1)]
    terms = String[]
    @inbounds for i in 0:n
        ci = coeffs[i + 1]
        iszero(ci) && continue
        cs = "(" * string(ci) * ")"
        if i == 0
            push!(terms, cs)
        elseif i == 1
            push!(terms, cs * "*" * dname)
        else
            push!(terms, cs * "*" * dname * "^" * string(i))
        end
    end
    isempty(terms) && error("computed zero operator from minimal polynomial")
    return parse_OrePoly(join(terms, " + "), A)
end

@inline function _dfinite_is_multiple_of_monvec(a::Vector{Int}, b::Vector{Int})
    @inbounds for i in 1:length(a)
        a[i] < b[i] && return false
    end
    return true
end

function _dfinite_monomial_string(m::Vector{Int}, dep_dnames::Vector{String}, ndep::Int)
    parts = String[]
    @inbounds for i in 1:ndep
        ei = m[i]
        ei == 0 && continue
        di = dep_dnames[i]
        if ei == 1
            push!(parts, di)
        else
            push!(parts, di * "^" * string(ei))
        end
    end
    return isempty(parts) ? "1" : join(parts, "*")
end

function _dfinite_term_string(c, m::Vector{Int}, dep_dnames::Vector{String}, ndep::Int)
    iszero(c) && return ""
    cs = "(" * string(c) * ")"
    ms = _dfinite_monomial_string(m, dep_dnames, ndep)
    return ms == "1" ? cs : cs * "*" * ms
end

function _dfinite_leaf_algebraic_ann_from_minpoly(expr::Expr, A::OreAlg, params::Dict{Symbol,<:Integer})
    p = minimal_polynomial(expr, A)
    pxy, base_syms = _dfinite_project_minpoly_to_base_tmp(p, A)

    nactive = div(A.nrdv, 2)
    nactive == 0 && (nactive = A.nrdv)
    gens = OrePoly{eltype1_ctx(ctx(A)),eltype_mo(A)}[]
    dep_idxs = Int[]
    dep_dnames = String[]
    @inbounds for j in 1:nactive
        xname = Symbol(A.inp.ratdiffvars[1][j])
        dname = A.inp.ratdiffvars[2][j]
        if !_dfinite_expr_depends_on_symbol(expr, xname)
            push!(gens, parse_OrePoly(dname, A))
            continue
        end
        push!(dep_idxs, j)
        push!(dep_dnames, dname)
    end

    ndep = length(dep_idxs)
    ndep == 0 && return gens

    up, _, K, U, t = _dfinite_minpoly_to_univariate_over_frac(pxy, base_syms)
    n = Nemo.degree(up)
    n > 0 || error("minimal polynomial has degree 0 in __tmp")

    lc = Nemo.coeff(up, n)
    p = lc == one(K) ? up : Nemo.divexact(up, lc)

    pt = derivative(p)
    yt = _dfinite_poly_to_vec(rem(t, p), n, K)
    ptv = _dfinite_poly_to_vec(rem(pt, p), n, K)
    invptv = _dfinite_inv_vec_mod_p(ptv, p, n, U, t, K)

    KT = typeof(one(K))
    dy_all = Vector{Vector{KT}}(undef, ndep)
    @inbounds for jj in 1:ndep
        xidx = dep_idxs[jj]
        px = zero(U)
        for i in 0:n
            ci = Nemo.coeff(p, i)
            iszero(ci) && continue
            px += U(derivative(ci, xidx)) * t^i
        end
        pxv = _dfinite_poly_to_vec(rem(px, p), n, K)
        dy_all[jj] = _dfinite_mul_vec_mod_p(invptv, map(-, pxv), p, n, U, t, K)
    end

    mon_basis = Vector{Vector{Int}}()
    vec_basis = Vector{Vector{KT}}()
    leading_mons = Vector{Vector{Int}}()
    normal_forms = Dict{Vector{Int},Vector{KT}}()

    m0 = zeros(Int, ndep)
    normal_forms[m0] = yt
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

        v = get(normal_forms, m, nothing)
        if v === nothing
            found = false
            @inbounds for i in 1:ndep
                m[i] == 0 && continue
                mp = copy(m)
                mp[i] -= 1
                vp = get(normal_forms, mp, nothing)
                vp === nothing && continue
                vv = _dfinite_derivative_vec(vp, dy_all[i], p, dep_idxs[i], n, K, U, t)
                normal_forms[m] = vv
                v = vv
                found = true
                break
            end
            found || error("failed to construct derivative image in algebraic leaf ann")
        end

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

        M = matrix_space(K, n, k + 1)()
        @inbounds for j in 1:k, i in 1:n
            M[i, j] = vec_basis[j][i]
        end
        @inbounds for i in 1:n
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
        isempty(terms) && error("computed zero operator from minimal polynomial")
        push!(gens, parse_OrePoly(join(terms, " + "), A))
        push!(leading_mons, m)
    end

    return gens
end
