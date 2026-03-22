# Main function for dfinite_expr_to_ann 


function dfinite_expr_to_ann(expr_in::Union{String,Expr}, A::OreAlg; ratvars::Union{Vector{Symbol},Vector{String}} = Symbol[])
    isempty(ratvars) || error("`ratvars` is only used when the algebra is inferred automatically")
    expr = expr_in isa String ? Meta.parse(expr_in, raise=true) : expr_in
    if expr isa Expr && expr.head == :incomplete
        error("Input expression is incomplete or malformed")
    end
    Ac = _dfinite_temp_closure_alg(A)
    gens = _dfinite_expr_to_ann(expr, Ac)
    clear_denominators!(gens, Ac)
    return map_algebras(gens, Ac, A)
end

function dfinite_expr_to_ann(expr_in::Union{String,Expr}; ratvars::Union{Vector{Symbol},Vector{String}} = Symbol[])
    expr = expr_in isa String ? Meta.parse(expr_in, raise=true) : expr_in
    if expr isa Expr && expr.head == :incomplete
        error("Input expression is incomplete or malformed")
    end
    ratdiffvars = _infer_ratdiffvars(expr)
    isempty(ratdiffvars) && error("Could not infer any differential variables from expression")
    ratvar_names = ratvars isa Vector{Symbol} ? String.(ratvars) : ratvars
    filter!(v -> !(v in ratvar_names), ratdiffvars)
    A = dfinite_ore_alg(ratdiffvars; ratvars = ratvar_names)
    Ac = _dfinite_temp_closure_alg(A)
    gens = _dfinite_expr_to_ann(expr, Ac)
    clear_denominators!(gens, Ac)
    return map_algebras(gens, Ac, A), A
end

function _dfinite_expr_to_ann(expr::Expr, A::OreAlg)
    expr.head == :call || error("Unsupported expression node $(expr.head). Only :call nodes are supported")
    op = expr.args[1]
    
    if _dfinite_expr_is_rational(expr)
        return _dfinite_rational_ann(expr, A)
    elseif _dfinite_expr_is_alg(expr, A.inp.ratvars)
        return _dfinite_leaf_algebraic_ann(expr, A)
    elseif op == :+ || op == :-
        return _dfinite_parse_sum_ann(expr, A)
    elseif op == :*
        return _dfinite_parse_product_ann(expr, A)
    elseif op == :^
        return _dfinite_parse_power_ann(expr, A)
    elseif op == :/ || op == ://
        return _dfinite_parse_division_ann(expr, A)
    end

    return _dfinite_leaf_ann(expr, A)
end

function _dfinite_rational_ann(expr, A::OreAlg)
    g_rat = _dfinite_expr_to_ratfun(expr, A)
    nactive = div(A.nrdv, 2)
    gens = Vector{OrePoly{eltype1_ctx(ctx(A)),eltype_mo(A)}}(undef, nactive)

    @inbounds for j in 1:nactive
        xname = A.inp.ratdiffvars[1][j]
        dname = A.inp.ratdiffvars[2][j]
        idx = A.strvar_to_indexp[dname]
        dg = _ann_ratfun_derivative(g_rat, xname, A)
        gens[j] = sub!(makepoly(one(ctx(A)), makemon(idx, A)), makepoly(dg / g_rat, makemon(-1, A)), A)
    end

    return gens
end

function _dfinite_expr_to_ann(sym::Symbol, A::OreAlg)
    nactive = div(A.nrdv, 2)
    return [parse_OrePoly(A.inp.ratdiffvars[1][i] *"*"* A.inp.ratdiffvars[2][i]* " - 1", A) for i in 1:nactive]
end


function _dfinite_expr_to_ann(x::Number, A::OreAlg)
    nactive = div(A.nrdv, 2)
    nactive == 0 && (nactive = A.nrdv)
    nactive > 0 || error("Cannot build annihilator of constant without differential variables")
    return [parse_OrePoly(A.inp.ratdiffvars[2][i], A) for i in 1:nactive]
end

function _dfinite_expr_to_ann(x, ::OreAlg)
    error("Unsupported node type $(typeof(x))")
end

function dfinite_ore_alg(ratdiffvars::Vector{String};
                         ratvars::Vector{String} = String[],
                         ratdiffops::Vector{String} = ["d" * v for v in ratdiffvars],
                         order::String = "",
                         char::Int = 0,
                         fraction_free::Bool = false,
                         varord::String = "dright")
    length(ratdiffvars) == length(ratdiffops) || error("ratdiffvars and ratdiffops must have the same length")
    ord = isempty(order) ? "grevlex " * join(vcat(ratdiffops), " ") : order
    return OreAlg(order = ord,
                char = char,
                ratvars = copy(ratvars),
                ratdiffvars = (copy(ratdiffvars), copy(ratdiffops)),
                fraction_free = fraction_free,
                varord = varord)
end

function _dfinite_parse_sum_ann(expr::Expr, A::OreAlg)
    op = expr.args[1]
    if op == :- && length(expr.args) == 2
        return _dfinite_expr_to_ann(expr.args[2], A)
    end
    length(expr.args) >= 3 || error("$(op) has $(length(expr.args)-1) operands")
    cur = _dfinite_expr_to_ann(expr.args[2], A)
    for i in 3:length(expr.args)
        rhs = _dfinite_expr_to_ann(expr.args[i], A)
        cur = ann_sum(cur, rhs, A)
    end
    return cur
end

function _dfinite_parse_product_ann(expr::Expr, A::OreAlg)
    length(expr.args) >= 3 || error("* node must have at least 2 operands")
    cur = _dfinite_expr_to_ann(expr.args[2], A)
    for i in 3:length(expr.args)
        rhs = _dfinite_expr_to_ann(expr.args[i], A)
        cur = ann_product(cur, rhs, A)
    end
    return cur
end

function _dfinite_parse_power_ann(expr::Expr, A::OreAlg)
    length(expr.args) == 3 || error("^ expects two arguments")
    base = expr.args[2]
    pow_arg = expr.args[3]
    pow = _ann_parse_exponent(pow_arg, A.inp.ratvars)
    if pow isa Int
        return _dfinite_power_ann(_dfinite_expr_to_ann(base, A), pow, A)
    elseif pow isa Rational{Int}
        return _dfinite_leaf_algebraic_ann(expr, A)
    elseif pow isa Symbol
        String(pow) in A.inp.ratvars || error("Symbolic exponents must be a ratvar symbol")
        return ann_poly_power(base, pow, A)
    elseif _dfinite_expr_depends_on_ratvars(pow_arg, A)
        error("Polynomial powers currently support only bare symbolic exponents, not mixed expressions like $(pow_arg)")
    end
    error("Only integer, rational exponents are currently supported in parser power nodes")
end

function _dfinite_parse_division_ann(expr::Expr, A::OreAlg)
    length(expr.args) == 3 || error("/ expects two arguments")
    num = _dfinite_expr_to_ann(expr.args[2], A)
    den = expr.args[3]
    return _dfinite_divide_by_hyperexp_ann(num, den, A)
end

# todo: do something about negative powers
function _dfinite_power_ann(base::Vector{OrePoly{T,M}}, n::Int, A::OreAlg) where {T,M}
    n >= 0 || error("Negative powers not handelt yet")
    if n == 0
        return [parse_OrePoly(A.inp.ratdiffvars[2][1], A)]
    elseif n == 1
        return base
    end

    result = [parse_OrePoly(A.inp.ratdiffvars[2][1], A)]
    power = base
    exp = n

    while exp > 0
        if isodd(exp)
            result = ann_product(result, power, A)
        end
        exp >>= 1
        if exp > 0
            power = ann_product(power, power, A)
        end
    end

    return result
end

function _dfinite_leaf_ann(expr::Expr, A::OreAlg)
    f = expr.args[1]
    if _dfinite_is_hypergeometric_pfq_symbol(f)
        return _dfinite_leaf_hypergeometric_pfq_ann(expr, A)
    end
    if f isa Symbol && haskey(datab_LDE, Symbol(String(f)))
        return _dfinite_leaf_db_ann(expr, A)
    end
    return _dfinite_leaf_algebraic_ann(expr, A)
end

function _dfinite_leaf_db_ann(expr::Expr, A::OreAlg)
    f = expr.args[1]
    f isa Symbol || error("Leaf function name must be a Symbol")
    key = Symbol(String(f))
    haskey(datab_LDE, key) || error("Function $(key) is not in datab_LDE")

    db_expr = datab_LDE[key]
    nargs = length(expr.args) - 1
    param_syms = get(datab_LDE_params, key, Symbol[])
    nargs == length(param_syms) + 1 || error("Unsupported arity $(nargs) for function $(key)")

    db_params = Dict{Symbol,Any}()
    @inbounds for i in eachindex(param_syms)
        parg = expr.args[i + 1]
        if parg isa Integer
            db_params[param_syms[i]] = Int(parg)
        elseif _dfinite_expr_is_ratvar_rational(parg, A)
            db_params[param_syms[i]] = parg
        else
            error("Only integer literal parameters or rational expressions in ratvars are supported for $(key)")
        end
    end
    xarg = expr.args[end]

    if xarg isa Symbol
        xname = String(xarg)
        idx = findfirst(==(xname), A.inp.ratdiffvars[1])
        idx === nothing && error("Variable $(xname) is not in ratdiffvars $(A.inp.ratdiffvars[1])")
        dname = A.inp.ratdiffvars[2][idx]

        subs = Dict{Symbol,Any}(:_x => Symbol(xname), :_D => Symbol(dname))
        for (k, v) in db_params
            subs[k] = v
        end
        lexpr = _subs_expr(db_expr, subs)
        gens = OrePoly{eltype1_ctx(ctx(A)),eltype_mo(A)}[]
        push!(gens, parse_OrePoly(string(lexpr), A))

        n_active = div(A.nrdv,2)
        for i in 1:n_active
            dother = A.inp.ratdiffvars[2][i]
            dother == dname && continue
            push!(gens, parse_OrePoly(dother, A))
        end
        return gens
    end

    xarg isa Expr || error("Function $(key) expects a symbolic or expression argument")

    if _dfinite_expr_is_rational(xarg)
        comp = ann_comp_right_rat(key, xarg, A; params = db_params)
        return comp
    else #it should be algebraic : check or return an error
        comp = _dfinite_compose_db_with_algebraic_ann(key, xarg, A, db_params)
        return comp
    end
end

function _dfinite_leaf_algebraic_ann(expr::Expr, A::OreAlg)
    return _dfinite_leaf_algebraic_ann_from_minpoly(expr, A, Dict{Symbol,Int}())
end

function _dfinite_is_hypergeometric_pfq_symbol(f)
    f isa Symbol || return false
    s = String(f)
    return s == "hypergeometric_pfq" || s == "HypergeometricPFQ"
end

function _dfinite_parse_param_sequence(arg, A::OreAlg, key::Symbol)
    if arg isa Expr && (arg.head === :tuple || arg.head === :vect)
        vals = Any[arg.args[i] for i in eachindex(arg.args)]
    else
        vals = Any[arg]
    end

    out = Vector{Any}(undef, length(vals))
    @inbounds for i in eachindex(vals)
        v = vals[i]
        if v isa Integer
            out[i] = Int(v)
        elseif _dfinite_expr_is_ratvar_rational(v, A)
            out[i] = v
        else
            error("Only integer literal parameters or rational expressions in ratvars are supported for $(key)")
        end
    end
    return out
end

function _dfinite_pfq_term_product(xd, vals::Vector{Any}, shift_one::Bool, include_xd::Bool)
    prod = include_xd ? xd : 1
    @inbounds for v in vals
        term = shift_one ? :($xd + $v - 1) : :($xd + $v)
        prod = :($prod * $term)
    end
    return prod
end

function _dfinite_pfq_operator_expr(avalues::Vector{Any}, bvalues::Vector{Any})
    xd = :(_x * _D)
    left = _dfinite_pfq_term_product(xd, bvalues, true, true)
    right = _dfinite_pfq_term_product(xd, avalues, false, false)
    return :($left - _x * $right)
end

function _dfinite_leaf_operator_ann(op_expr::Expr, xarg, A::OreAlg)
    if xarg isa Symbol
        xname = String(xarg)
        idx = findfirst(==(xname), A.inp.ratdiffvars[1])
        idx === nothing && error("Variable $(xname) is not in ratdiffvars $(A.inp.ratdiffvars[1])")
        dname = A.inp.ratdiffvars[2][idx]

        lexpr = _subs_expr(op_expr, Dict{Symbol,Any}(:_x => Symbol(xname), :_D => Symbol(dname)))
        gens = OrePoly{eltype1_ctx(ctx(A)),eltype_mo(A)}[parse_OrePoly(string(lexpr), A)]

        n_active = div(A.nrdv,2)
        for i in 1:n_active
            dother = A.inp.ratdiffvars[2][i]
            dother == dname && continue
            push!(gens, parse_OrePoly(dother, A))
        end
        return gens
    end

    xarg isa Expr || error("Function expects a symbolic or expression argument")
    if _dfinite_expr_is_rational(xarg)
        nactive = div(A.nrdv, 2)
        g_rat = _dfinite_expr_to_ratfun(xarg, A)
        res = OrePoly{eltype1_ctx(ctx(A)),eltype_mo(A)}[]
        for j in 1:nactive
            tmp = ann_comp_right_rat_ev_LDE(A, xarg, j, op_expr)
            xname = A.inp.ratdiffvars[1][j]
            idx = A.strvar_to_indexp[A.inp.ratdiffvars[2][j]]
            dg = ctx(A) isa MRatFunCtx ? derivative(g_rat, ctx(A).vars[A.drvars_to_int[xname]]) : derivative(g_rat)
            push!(res, ann_comp_right_rat_ev_D(A, dg, idx, tmp))
        end
        return res
    end

    L = _dfinite_operator_first_active(op_expr, A)
    return ann_comp_right_alg(A, L, xarg)
end

function _dfinite_leaf_hypergeometric_pfq_ann(expr::Expr, A::OreAlg)
    length(expr.args) == 4 || error("hypergeometric_pfq expects three arguments")
    key = :hypergeometric_pfq
    avalues = _dfinite_parse_param_sequence(expr.args[2], A, key)
    bvalues = _dfinite_parse_param_sequence(expr.args[3], A, key)
    op_expr = _dfinite_pfq_operator_expr(avalues, bvalues)
    return _dfinite_leaf_operator_ann(op_expr, expr.args[4], A)
end

function _dfinite_compose_db_with_algebraic_ann(key::Symbol, expr::Expr, A::OreAlg, params::Dict{Symbol,Any})
    L = _dfinite_db_operator_first_active(key, A, params)
    return ann_comp_right_alg(A, L, expr)
end

function _dfinite_divide_by_hyperexp_ann(num::Vector{OrePoly{T,M}}, den :: Expr, A::OreAlg) where {T,M}
    return _dfinite_divide_by_hyperexp_ann_main(num, den, A)
end

function _subs_expr(x::Symbol, subs::Dict{Symbol,Any})
    return get(subs, x, x)
end

_subs_expr(x::Number, ::Dict{Symbol,Any}) = x

function _subs_expr(x::Expr, subs::Dict{Symbol,Any})
    return Expr(x.head, (_subs_expr(arg, subs) for arg in x.args)...)
end

_infer_ratdiffvars(expr :: Symbol) = [String(expr)]
_infer_ratdiffvars(::Number) = String[]
_infer_ratdiffvars(::Nothing) = String[]

function _infer_ratdiffvars(expr_in::Expr)
    vars = Set{String}()
    _infer_ratdiffvars!(vars, expr_in)
    return sort!(collect(vars))
end

function _infer_ratdiffvars!(vars::Set{String}, expr::Expr)
    if expr.head != :call
        return nothing
    end
    op = expr.args[1]
    if op == :+ || op == :- || op == :* || op == :^ || op == :/ || op == ://
        for i in 2:length(expr.args)
            arg = expr.args[i]
            _infer_vars_from_arg!(vars, arg)
        end
        return nothing
    end
    nargs = length(expr.args) - 1
    if nargs >= 1
        xarg = expr.args[end]
        _infer_vars_from_arg!(vars, xarg)
    end
    return nothing
end

function _infer_vars_from_arg!(vars::Set{String}, x)
    if x isa Symbol
        push!(vars, String(x))
        return nothing
    elseif x isa Number
        return nothing
    elseif x isa Expr
        if x.head == :call
            for i in 2:length(x.args)
                _infer_vars_from_arg!(vars, x.args[i])
            end
        else
            for a in x.args
                _infer_vars_from_arg!(vars, a)
            end
        end
        return nothing
    end
    return nothing
end

function _dfinite_expr_is_rational(x)
    if x isa Number || x isa Symbol
        return true
    elseif !(x isa Expr)
        return false
    end
    x.head == :call || return false
    op = x.args[1]
    if op == :+ || op == :* || op == :- || op == :/ || op == ://
        for i in 2:length(x.args)
            _dfinite_expr_is_rational(x.args[i]) || return false
        end
        return true
    elseif op == :^
        length(x.args) == 3 || return false
        _dfinite_expr_is_rational(x.args[2]) || return false
        expo = x.args[3]
        return (expo isa Integer) || (expo isa Rational{<:Integer})
    end
    return false
end

function _dfinite_expr_is_alg(x, ratvars::Vector{String} = String[])
    if x isa Number || x isa Symbol
        return true
    elseif !(x isa Expr)
        return false
    end
    x.head == :call || return false
    args = x.args
    op = args[1]
    if op == :+ || op == :- || op == :* || op == :x
        @inbounds for i in 2:length(args)
            _dfinite_expr_is_alg(args[i], ratvars) || return false
        end
        return true
    elseif op == :/ || op == ://
        length(args) == 3 || return false
        return _dfinite_expr_is_alg(args[2], ratvars) && _dfinite_expr_is_alg(args[3], ratvars)
    elseif op == :^
        length(args) == 3 || return false
        _dfinite_expr_is_alg(args[2], ratvars) || return false
        expo = _ann_parse_exponent(args[3], ratvars)
        return (expo isa Int) || (expo isa Rational{Int})
    end
    return false
end

function _dfinite_expr_is_ratvar_rational(x, A::OreAlg)
    if x isa Number
        return true
    elseif x isa Symbol
        return String(x) in A.inp.ratvars
    elseif !(x isa Expr)
        return false
    end
    x.head == :call || return false
    op = x.args[1]
    if op == :+ || op == :* || op == :- || op == :/ || op == ://
        for i in 2:length(x.args)
            _dfinite_expr_is_ratvar_rational(x.args[i], A) || return false
        end
        return true
    elseif op == :^
        length(x.args) == 3 || return false
        _dfinite_expr_is_ratvar_rational(x.args[2], A) || return false
        expo = x.args[3]
        return (expo isa Integer) || (expo isa Rational{<:Integer})
    end
    return false
end

function _dfinite_expr_depends_on_ratvars(x, A::OreAlg)
    if x isa Symbol
        return String(x) in A.inp.ratvars
    elseif x isa Expr
        for arg in x.args
            _dfinite_expr_depends_on_ratvars(arg, A) && return true
        end
    end
    return false
end

#todo: check what it does
function _dfinite_db_operator_first_active(key::Symbol, A::OreAlg, params::Dict{Symbol,Any})
    db_expr = datab_LDE[key]
    return _dfinite_operator_first_active(db_expr, A, params)
end

function _dfinite_operator_first_active(db_expr::Expr, A::OreAlg, params::Dict{Symbol,Any} = Dict{Symbol,Any}())
    xname = A.inp.ratdiffvars[1][1]
    dname = A.inp.ratdiffvars[2][1]
    subs = Dict{Symbol,Any}(:_x => Symbol(xname), :_D => Symbol(dname))
    for (k, v) in params
        subs[k] = v
    end
    lexpr = _subs_expr(db_expr, subs)
    return parse_OrePoly(string(lexpr), A)
end

function _dfinite_temp_closure_alg(A::OreAlg)
    # add one variable for closure by sum
    tvar = "_t"
    _dfinite_is_free_name(A, tvar) || error("Cannot create temporary closure algebra: variable name $(tvar) is already used")
    inp = deepcopy(A.inp)
    inp.fraction_free = false
    n = length(inp.ratdiffvars[1])
    # double the number of variable for closure by product
    for i in 1:n
        xv2 = inp.ratdiffvars[1][i] * "__2"
        dv2 = inp.ratdiffvars[2][i] * "__2"
        _dfinite_is_free_name(A, xv2) || error("Cannot create temporary closure algebra: variable name $(xv2) is already used")
        _dfinite_is_free_name(A, dv2) || error("Cannot create temporary closure algebra: variable name $(dv2) is already used")
        push!(inp.ratdiffvars[1], xv2)
        push!(inp.ratdiffvars[2], dv2)
    end
    push!(inp.polvars, tvar)
    push!(inp.nomul, tvar)
    if n == 0
        inp.order = "lex " * tvar
    else
        inp.order = "lex " * tvar * " > grevlex " * join(inp.ratdiffvars[2], " ")
    end
    return OreAlg(inp)
end

function _dfinite_is_free_name(A::OreAlg, s::String)
    inp = A.inp
    return !(s in inp.ratvars ||
             s in inp.ratdiffvars[1] ||
             s in inp.ratdiffvars[2] ||
             s in inp.poldiffvars[1] ||
             s in inp.poldiffvars[2] ||
             s in inp.polvars ||
             s in inp.locvars[1] ||
             s in inp.nomul)
end
