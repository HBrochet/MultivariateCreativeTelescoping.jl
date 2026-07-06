function _ann_parse_exponent(expo, ratvars::Vector{String} = String[])
    if expo isa Integer
        return Int(expo)
    elseif expo isa Rational{<:Integer}
        return Int(Base.numerator(expo)) // Int(Base.denominator(expo))
    elseif expo isa Symbol
        s = String(expo)
        !(s in ratvars) && error("Symbolic exponents of polynomials should be declared in ratvars")
        return expo
    elseif expo isa Expr && expo.head === :call && (expo.args[1] === :// || expo.args[1] === :/) && length(expo.args) == 3
        num = _ann_parse_exponent(expo.args[2], ratvars)
        den = _ann_parse_exponent(expo.args[3], ratvars)
        if (num isa Int || num isa Rational{Int}) && (den isa Int || den isa Rational{Int})
            return num // den
        end
    end
    return nothing
end

function _ann_alg_count_aux(e)
    if e isa Number || e isa Symbol
        return 0
    elseif !(e isa Expr)
        throw(ArgumentError("unsupported node type in algebraic expression"))
    end

    e.head === :call || throw(ArgumentError("only :call Expr nodes are supported"))
    args = e.args
    op = args[1]

    if op === :+ || op === :- || op === :* || op === :x
        c = 0
        @inbounds for i in 2:length(args)
            c += _ann_alg_count_aux(args[i])
        end
        return c
    elseif op === :/ || op === ://
        length(args) == 3 || throw(ArgumentError("/ expects exactly 2 operands"))
        return 1 + _ann_alg_count_aux(args[2]) + _ann_alg_count_aux(args[3])
    elseif op === :^
        length(args) == 3 || throw(ArgumentError("^ expects exactly 2 operands"))
        base = args[2]
        expo = _ann_parse_exponent(args[3])
        c = _ann_alg_count_aux(base)
        if expo isa Int
            return expo < 0 ? c + 1 : c
        elseif expo isa Rational{Int}
            return c + 1
        else
            throw(ArgumentError("power must be an integer or rational"))
        end
    end

    throw(ArgumentError("unsupported internal node $op"))
end

function _ann_alg_parse_expr!(e, R, sym_to_var, aux_vars, next_aux::Base.RefValue{Int}, aux_eqs, aux_cache::Dict{Expr,Any})
    if e isa Number
        return R(e)
    elseif e isa Symbol
        v = get(sym_to_var, e, nothing)
        v === nothing && throw(ArgumentError("unknown symbol in algebraic expression"))
        return v
    elseif !(e isa Expr)
        throw(ArgumentError("unsupported node type in algebraic expression"))
    end

    e.head === :call || throw(ArgumentError("only :call Expr nodes are supported"))
    args = e.args
    op = args[1]

    if op === :+
        length(args) >= 3 || throw(ArgumentError("+ expects at least 2 operands"))
        acc = _ann_alg_parse_expr!(args[2], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        @inbounds for i in 3:length(args)
            acc += _ann_alg_parse_expr!(args[i], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        end
        return acc
    elseif op === :-
        length(args) >= 2 || throw(ArgumentError("- expects at least 1 operand"))
        if length(args) == 2
            return -_ann_alg_parse_expr!(args[2], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        end
        acc = _ann_alg_parse_expr!(args[2], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        @inbounds for i in 3:length(args)
            acc -= _ann_alg_parse_expr!(args[i], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        end
        return acc
    elseif op === :* || op === :x
        length(args) >= 3 || throw(ArgumentError("* expects at least 2 operands"))
        acc = _ann_alg_parse_expr!(args[2], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        @inbounds for i in 3:length(args)
            acc *= _ann_alg_parse_expr!(args[i], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        end
        return acc
    elseif op === :/ || op === ://
        length(args) == 3 || throw(ArgumentError("/ expects exactly 2 operands"))
        v = get(aux_cache, e, nothing)
        v === nothing || return v
        num = _ann_alg_parse_expr!(args[2], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        den = _ann_alg_parse_expr!(args[3], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        next_aux[] += 1
        u = aux_vars[next_aux[]]
        push!(aux_eqs, u * den - num)
        aux_cache[e] = u
        return u
    elseif op === :^
        length(args) == 3 || throw(ArgumentError("^ expects exactly 2 operands"))
        base = _ann_alg_parse_expr!(args[2], R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)
        pow = _ann_parse_exponent(args[3])

        if pow isa Int
            if pow >= 0
                return base^pow
            end
            v = get(aux_cache, e, nothing)
            v === nothing || return v
            next_aux[] += 1
            u = aux_vars[next_aux[]]
            push!(aux_eqs, u * (base^(-pow)) - one(R))
            aux_cache[e] = u
            return u
        elseif pow isa Rational{Int}
            v = get(aux_cache, e, nothing)
            v === nothing || return v
            p = Base.numerator(pow)
            q = Base.denominator(pow)
            next_aux[] += 1
            u = aux_vars[next_aux[]]
            if p >= 0
                push!(aux_eqs, u^q - base^p)
            else
                push!(aux_eqs, (u^q) * (base^(-p)) - one(R))
            end
            aux_cache[e] = u
            return u
        end

        throw(ArgumentError("power must be an integer or rational"))
    end

    throw(ArgumentError("unsupported internal node"))
end

function _ann_alg_embed_small_poly_to_big(s, Rbig, big_vars, rem_to_big)
    out = zero(Rbig)
    ls = length(s)
    @inbounds for i in 1:ls
        ev = Nemo.exponent_vector(s, i)
        term = Rbig(Nemo.coeff(s, i))
        for j in 1:length(ev)
            ej = ev[j]
            ej == 0 && continue
            term *= big_vars[rem_to_big[j]]^ej
        end
        out += term
    end
    return out
end

function _ann_alg_to_univariate_over_poly(p, var_idx::Int, U, t, S0, rem_vars, rem_to_big)
    up = zero(U)
    lp = length(p)
    @inbounds for i in 1:lp
        ev = Nemo.exponent_vector(p, i)
        d = ev[var_idx]
        c = S0(Nemo.coeff(p, i))
        for j in 1:length(rem_vars)
            ej = ev[rem_to_big[j]]
            ej == 0 && continue
            c *= rem_vars[j]^ej
        end
        up += U(c) * t^d
    end
    return up
end

function _ann_alg_eliminate_var_by_resultant(p, q, var_idx::Int, all_syms)
    Rbig = parent(p)
    big_vars = gens(Rbig)
    rem_to_big = [i for i in 1:length(big_vars) if i != var_idx]
    rem_syms = all_syms[rem_to_big]

    S0, rem_vars = polynomial_ring(base_ring(Rbig), rem_syms)
    U, t = polynomial_ring(S0, Symbol(all_syms[var_idx]))

    up = _ann_alg_to_univariate_over_poly(p, var_idx, U, t, S0, rem_vars, rem_to_big)
    uq = _ann_alg_to_univariate_over_poly(q, var_idx, U, t, S0, rem_vars, rem_to_big)
    ur = resultant(up, uq)
    iszero(ur) && return zero(Rbig)
    return _ann_alg_embed_small_poly_to_big(ur, Rbig, big_vars, rem_to_big)
end

function minimal_polynomial(e :: Expr, A :: OreAlg)
    raw_syms = ctx(A).R.S
    base_syms = raw_syms isa Symbol ? Symbol[raw_syms] : Symbol.(raw_syms)
    nbase = length(base_syms)
    naux = _ann_alg_count_aux(e)
    aux_syms = [Symbol("__aux", i) for i in 1:naux]
    all_syms = vcat(base_syms, aux_syms, [:__tmp])

    R, vars = polynomial_ring(base_ring(ctx(A).R), all_syms)

    sym_to_var = Dict{Symbol,Any}(base_syms[i] => vars[i] for i in 1:nbase)
    aux_vars = vars[nbase + 1 : nbase + naux]
    y = vars[end]
    #todo: what you do with next_aux is very strange. It does not look like julia programming. Can you explain ?
    aux_eqs = Vector{typeof(vars[1])}()
    next_aux = Ref(0)
    aux_cache = Dict{Expr,Any}()
    rhs = _ann_alg_parse_expr!(e, R, sym_to_var, aux_vars, next_aux, aux_eqs, aux_cache)

    p = y - rhs
    @inbounds for i in next_aux[]:-1:1
        p = _ann_alg_eliminate_var_by_resultant(p, aux_eqs[i], nbase + i, all_syms)
    end

    pclean, _ = _dfinite_project_minpoly_to_base_tmp(p, A)
    return pclean
end
