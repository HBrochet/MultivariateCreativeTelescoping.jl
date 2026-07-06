"""
    compute_module_groebner_basis(
        gens::Vector{<:OrePoly},
        A0::OreAlg,
        ratvars::Vector{String},
        polvars::Vector{String},
        trunc_order::Integer;
        ordering = nothing,
        paramgb_kwargs...
    )

Build the commutative module presentation attached to the annihilators `gens`
returned by `dfinite_expr_to_ann`, truncate it at differential order
`trunc_order`, then compute a Groebner basis with `ParamPunPam.paramgb`.

The module trick is encoded as an ideal:

- a module variable `e_{d^alpha}` is introduced for each differential monomial
  `d^alpha` of total order at most `trunc_order`;
- `ed0` denotes the basis vector corresponding to `d^0`;
- all quadratic products of module variables are set to zero;
- each prolonged operator `d^beta * g` is converted to a linear polynomial in
  the module variables.

`ratvars` and `polvars` describe the target coefficient ring
`Q(ratvars)[polvars]`. They must contain all coefficient variables used in `A0`.

The keyword `ordering` can either be:

- `nothing`, in which case a default `DegRevLex(polvars..., module_vars...)`
  ordering is used;
- a parsable block ordering string such as
  `"lex x > grevlex ed0 edx edx2"`;
- a `Groebner.jl`/`ParamPunPam.jl` ordering object.

The returned named tuple contains the Groebner basis together with the target
ring, generators, ordering and module-variable metadata.
"""
function compute_module_groebner_basis(
    gens::Vector{<:OrePoly},
    A0::OreAlg,
    ratvars::Vector{String},
    polvars::Vector{String},
    trunc_order::Integer;
    ordering = nothing,
    paramgb_kwargs...
)
    trunc_order < 0 && error("trunc_order must be non-negative")
    _validate_module_groebner_inputs(A0, ratvars, polvars)

    diffops = copy(A0.inp.ratdiffvars[2])
    ndiff = length(diffops)
    basis_exponents = _module_basis_exponents(ndiff, Int(trunc_order))
    module_var_names = [_module_var_name(diffops, exponents) for exponents in basis_exponents]
    _validate_module_var_collisions(ratvars, polvars, module_var_names)

    ring_data = _build_module_target_ring(ratvars, polvars, module_var_names)
    return_var_lookup = ring_data.return_var_lookup
    source_coeff_substitution = _build_source_coeff_substitution(A0, ring_data.symbol_lookup)
    module_var_map = Dict(basis_exponents[i] => ring_data.module_generators[i] for i in eachindex(basis_exponents))

    prolonged_generators = _prolong_generators(gens, A0, Int(trunc_order))
    linear_generators = [
        _convert_operator_to_module_polynomial(
            op,
            A0,
            ring_data.ring,
            source_coeff_substitution,
            module_var_map,
            Int(trunc_order),
        )
        for op in prolonged_generators
    ]

    module_vars = ring_data.module_generators
    square_zero = [module_vars[i] * module_vars[j] for i in eachindex(module_vars) for j in i:length(module_vars)]
    ideal_generators = vcat(square_zero, linear_generators)

    gb_ordering = _module_groebner_ordering(
        ordering,
        ring_data.polynomial_generators,
        ring_data.module_generators,
        ring_data.polynomial_var_lookup,
        ring_data.module_var_lookup,
    )
    basis = ParamPunPam.paramgb(ideal_generators; ordering = gb_ordering, paramgb_kwargs...)

    return (
        basis = basis,
        generators = ideal_generators,
        square_zero_generators = square_zero,
        linear_generators = linear_generators,
        prolonged_generators = prolonged_generators,
        ring = ring_data.ring,
        coefficient_field = ring_data.coefficient_field,
        parameter_ring = ring_data.parameter_ring,
        ordering = gb_ordering,
        module_var_names = module_var_names,
        module_monomials = basis_exponents,
        diffop_names = diffops,
        variables = return_var_lookup,
    )
end

compute_module_groebner(
    gens::Vector{<:OrePoly},
    A0::OreAlg,
    ratvars::Vector{String},
    polvars::Vector{String},
    trunc_order::Integer;
    kwargs...
) = compute_module_groebner_basis(gens, A0, ratvars, polvars, trunc_order; kwargs...)

"""
    compute_weyl_closure_module_groebner_basis(
        gens::Vector{<:OrePoly},
        A0::OreAlg,
        trunc_order::Integer;
        ratvars = A0.inp.ratvars,
        polvars = A0.inp.ratdiffvars[1],
        localizing_factor = :singular_locus,
        localization_var = "T",
        ordering = nothing,
        paramgb_kwargs...
    )

Build the localized commutative module presentation used by
`heuristic_weyl_closure_holonomic_ideal` and compute its parametric Groebner
basis with `ParamPunPam.paramgb`.

Compared with `compute_module_groebner_basis`, this adds the relations
`(f*T - 1)*e` for every module variable `e`, where `f` is either the square-free
singular-locus factor (`localizing_factor = :singular_locus`), an explicit
factor, or omitted with `localizing_factor = :none`. By default, the localized
module basis uses an elimination block `Lex(T, edt)` when `edt` is present.
"""
function compute_weyl_closure_module_groebner_basis(
    gens::Vector{<:OrePoly},
    A0::OreAlg,
    trunc_order::Integer;
    ratvars::Vector{String} = copy(A0.inp.ratvars),
    polvars::Vector{String} = copy(A0.inp.ratdiffvars[1]),
    localizing_factor = :singular_locus,
    localization_var::String = "T",
    ordering = nothing,
    paramgb_kwargs...
)
    trunc_order < 0 && error("trunc_order must be non-negative")
    _validate_module_groebner_inputs(A0, ratvars, polvars)

    diffops = copy(A0.inp.ratdiffvars[2])
    ndiff = length(diffops)
    basis_exponents = _module_basis_exponents(ndiff, Int(trunc_order))
    module_var_names = [_module_var_name(diffops, exponents) for exponents in basis_exponents]
    _validate_module_var_collisions(ratvars, polvars, module_var_names)

    source_localizing_factor = _module_source_localizing_factor(localizing_factor, gens, A0)
    use_localization = !_module_is_trivial_localizing_factor(source_localizing_factor)
    aux_var_names = String[]
    if use_localization
        used = Set(vcat(ratvars, polvars, module_var_names))
        push!(aux_var_names, _fresh_module_aux_var_name(localization_var, used))
    end

    ring_data = _build_module_target_ring(ratvars, polvars, module_var_names; aux_var_names)
    source_coeff_substitution = _build_source_coeff_substitution(A0, ring_data.symbol_lookup)
    module_var_map = Dict(basis_exponents[i] => ring_data.module_generators[i] for i in eachindex(basis_exponents))

    prolonged_generators = _prolong_generators(gens, A0, Int(trunc_order))
    linear_generators = [
        _convert_operator_to_module_polynomial(
            op,
            A0,
            ring_data.ring,
            source_coeff_substitution,
            module_var_map,
            Int(trunc_order),
        )
        for op in prolonged_generators
    ]

    module_vars = ring_data.module_generators
    square_zero = [module_vars[i] * module_vars[j] for i in eachindex(module_vars) for j in i:length(module_vars)]

    localization_generators = elem_type(ring_data.ring)[]
    target_localizing_factor = one(ring_data.ring)
    if use_localization
        target_localizing_factor = _module_localizing_factor_to_target_ring(
            source_localizing_factor,
            A0,
            ring_data.ring,
            source_coeff_substitution,
        )
        T = ring_data.aux_generators[1]
        loc_relation = target_localizing_factor * T - one(ring_data.ring)
        localization_generators = [loc_relation * e for e in module_vars]
    end

    ideal_generators = vcat(square_zero, localization_generators, linear_generators)
    gb_ordering = _module_weyl_closure_ordering(ordering, ring_data)
    basis = ParamPunPam.paramgb(ideal_generators; ordering = gb_ordering, paramgb_kwargs...)

    return (
        basis = basis,
        generators = ideal_generators,
        square_zero_generators = square_zero,
        localization_generators = localization_generators,
        linear_generators = linear_generators,
        prolonged_generators = prolonged_generators,
        ring = ring_data.ring,
        coefficient_field = ring_data.coefficient_field,
        parameter_ring = ring_data.parameter_ring,
        ordering = gb_ordering,
        module_var_names = module_var_names,
        module_monomials = basis_exponents,
        diffop_names = diffops,
        localization_var_name = use_localization ? aux_var_names[1] : nothing,
        localizing_factor = target_localizing_factor,
        source_localizing_factor = source_localizing_factor,
        variables = ring_data.return_var_lookup,
        ring_data = ring_data,
    )
end

"""
    heuristic_weyl_closure_holonomic_ideal(gens, A0; kwargs...)

Heuristically compute a holonomic ideal contained in the Weyl closure of the
annihilator `gens` returned by `dfinite_expr_to_ann`.

The function tries truncation orders from `min_order` to `max_order`. At each
order it computes a localized module Groebner basis with
`ParamPunPam.paramgb`, extracts the `T`-free linear module relations as Weyl
operators, computes their Weyl Groebner basis with `f4`, and stops as soon as
the result satisfies `isholonomic`.

The returned named tuple contains `ideal`, `generators`, `algebra`,
`trunc_order`, `holonomic`, and the underlying commutative `module_data`.
"""
function heuristic_weyl_closure_holonomic_ideal(
    gens::Vector{<:OrePoly},
    A0::OreAlg;
    min_order::Union{Nothing, Integer} = nothing,
    max_order::Union{Nothing, Integer} = nothing,
    ratvars::Vector{String} = copy(A0.inp.ratvars),
    polvars::Vector{String} = copy(A0.inp.ratdiffvars[1]),
    target_algebra::Union{Nothing, OreAlg} = nothing,
    target_order::Union{Nothing, String} = nothing,
    localizing_factor = :singular_locus,
    localization_var::String = "T",
    ordering = nothing,
    compute_weyl_groebner::Bool = true,
    return_last::Bool = false,
    paramgb_kwargs...
)
    _validate_module_groebner_inputs(A0, ratvars, polvars)
    start_order = isnothing(min_order) ? _max_operator_order(gens, A0) : Int(min_order)
    stop_order = isnothing(max_order) ? start_order : Int(max_order)
    stop_order >= start_order || error("max_order must be at least min_order")

    target_A = isnothing(target_algebra) ?
        _module_default_weyl_target_algebra(A0, ratvars, polvars; target_order) :
        target_algebra
    _validate_module_weyl_target_algebra(target_A, ratvars, polvars, A0.inp.ratdiffvars[2])

    last_result = nothing
    for trunc_order in start_order:stop_order
        module_data = compute_weyl_closure_module_groebner_basis(
            gens,
            A0,
            trunc_order;
            ratvars,
            polvars,
            localizing_factor,
            localization_var,
            ordering,
            paramgb_kwargs...
        )
        operators = _module_groebner_basis_to_weyl_operators(module_data, target_A)
        ideal = compute_weyl_groebner && !isempty(operators) ? f4(operators, target_A) : operators
        hol = !isempty(ideal) && isholonomic(ideal, target_A)
        result = (
            ideal = ideal,
            generators = operators,
            algebra = target_A,
            trunc_order = trunc_order,
            holonomic = hol,
            module_data = module_data,
        )
        hol && return result
        last_result = result
    end

    return_last && return last_result
    error("failed to find a holonomic ideal up to truncation order $(stop_order)")
end

function _validate_module_groebner_inputs(A0::OreAlg, ratvars::Vector{String}, polvars::Vector{String})
    length(unique(ratvars)) == length(ratvars) || error("ratvars must not contain duplicates")
    length(unique(polvars)) == length(polvars) || error("polvars must not contain duplicates")
    isempty(intersect(Set(ratvars), Set(polvars))) || error("ratvars and polvars must be disjoint")

    if (A0.npdv != 0) || (A0.npv != 0) || (A0.nlv != 0)
        error("compute_module_groebner_basis currently expects the algebra returned by automatic dfinite_expr_to_ann, i.e. only ratdiffvars and optional ratvars")
    end

    coeff_vars = unique(vcat(A0.inp.ratdiffvars[1], A0.inp.ratvars))
    supplied = Set(vcat(ratvars, polvars))
    missing = filter(v -> !(v in supplied), coeff_vars)
    isempty(missing) || error("missing coefficient variables in ratvars/polvars: $(join(missing, ", "))")
    return nothing
end

function _validate_module_var_collisions(ratvars::Vector{String}, polvars::Vector{String}, module_var_names::Vector{String})
    collisions = intersect(Set(vcat(ratvars, polvars)), Set(module_var_names))
    isempty(collisions) || error("module variable names collide with coefficient variables: $(join(collect(collisions), ", "))")
    return nothing
end

function _fresh_module_aux_var_name(preferred::String, used::Set{String})
    name = preferred
    ctr = 1
    while name in used
        name = preferred * "_" * string(ctr)
        ctr += 1
    end
    return name
end

function _module_basis_exponents(ndiff::Int, trunc_order::Int)
    ndiff > 0 || error("A0 must contain at least one differential variable")
    res = Vector{NTuple{ndiff, Int}}()
    current = zeros(Int, ndiff)

    function rec!(idx::Int, remaining::Int)
        if idx == ndiff
            current[idx] = remaining
            push!(res, Tuple(current))
            return nothing
        end
        for e in remaining:-1:0
            current[idx] = e
            rec!(idx + 1, remaining - e)
        end
        return nothing
    end

    for total_deg in 0:trunc_order
        rec!(1, total_deg)
    end
    return res
end

function _module_var_name(diffops::Vector{String}, exponents::NTuple)
    all(iszero, exponents) && return "ed0"
    parts = String[]
    for (dop, exponent) in zip(diffops, exponents)
        exponent == 0 && continue
        suffix = exponent == 1 ? "" : string(exponent)
        push!(parts, dop * suffix)
    end
    return "e" * join(parts, "_")
end

function _build_module_target_ring(
    ratvars::Vector{String},
    polvars::Vector{String},
    module_var_names::Vector{String};
    aux_var_names::Vector{String} = String[],
)
    parameter_var_names = isempty(ratvars) ? ["_ppp_dummy_param"] : ratvars
    parameter_ring, param_tuple = Nemo.polynomial_ring(Nemo.QQ, parameter_var_names)
    coefficient_field = Nemo.fraction_field(parameter_ring)
    rat_param_gens = collect(param_tuple)

    ring_var_names = vcat(aux_var_names, polvars, module_var_names)
    ring, ring_tuple = Nemo.polynomial_ring(coefficient_field, ring_var_names)
    ring_generators = collect(ring_tuple)
    naux = length(aux_var_names)
    npol = length(polvars)
    aux_gens = ring_generators[1:naux]
    pol_gens = ring_generators[(naux + 1):(naux + npol)]
    module_gens = ring_generators[(naux + npol + 1):end]
    RingElem = typeof(module_gens[1])

    aux_var_lookup = Dict(aux_var_names[i] => aux_gens[i] for i in eachindex(aux_var_names))
    polynomial_var_lookup = Dict(polvars[i] => pol_gens[i] for i in eachindex(polvars))
    module_var_lookup = Dict(module_var_names[i] => module_gens[i] for i in eachindex(module_var_names))
    symbol_lookup = Dict{String, RingElem}()
    return_var_lookup = Dict{String, RingElem}()

    for (name, param_gen) in zip(ratvars, rat_param_gens)
        const_poly = ring(coefficient_field(param_gen))
        symbol_lookup[name] = const_poly
        return_var_lookup[name] = const_poly
    end

    for (name, var) in pairs(polynomial_var_lookup)
        symbol_lookup[name] = var
        return_var_lookup[name] = var
    end
    for (name, var) in pairs(aux_var_lookup)
        return_var_lookup[name] = var
    end
    for (name, var) in pairs(module_var_lookup)
        return_var_lookup[name] = var
    end

    return (
        parameter_ring = parameter_ring,
        coefficient_field = coefficient_field,
        ring = ring,
        aux_generators = aux_gens,
        polynomial_generators = pol_gens,
        module_generators = module_gens,
        aux_var_lookup = aux_var_lookup,
        polynomial_var_lookup = polynomial_var_lookup,
        module_var_lookup = module_var_lookup,
        aux_var_names = aux_var_names,
        polynomial_var_names = polvars,
        module_var_names = module_var_names,
        symbol_lookup = symbol_lookup,
        return_var_lookup = return_var_lookup,
    )
end

function _module_groebner_ordering(
    ordering,
    polynomial_generators::Vector,
    module_generators::Vector,
    polynomial_var_lookup::Dict{String},
    module_var_lookup::Dict{String},
)
    if isnothing(ordering)
        vars = Any[v for v in polynomial_generators]
        append!(vars, module_generators)
        return ParamPunPam.DegRevLex(vars...)
    elseif ordering isa String
        return _module_groebner_ordering_from_string(ordering, merge(polynomial_var_lookup, module_var_lookup))
    else
        return ordering
    end
end

function _module_groebner_ordering_from_string(ordering::String, var_lookup::Dict{String})
    combined = nothing
    for block in split(ordering, ">")
        tokens = split(strip(block))
        isempty(tokens) && continue
        length(tokens) >= 2 || error("invalid ordering block: $(block)")
        ord_name = lowercase(tokens[1])
        vars = map(tokens[2:end]) do token
            get(var_lookup, token) do
                error("unknown variable $(token) in module Groebner ordering")
            end
        end
        block_order = if ord_name == "lex"
            ParamPunPam.Lex(vars...)
        elseif ord_name == "grevlex" || ord_name == "degrevlex"
            ParamPunPam.DegRevLex(vars...)
        elseif ord_name == "deglex"
            ParamPunPam.DegLex(vars...)
        else
            error("unsupported ordering block $(ord_name)")
        end
        combined = isnothing(combined) ? block_order : combined * block_order
    end
    isnothing(combined) && error("failed to parse ordering $(ordering)")
    return combined
end

function _module_weyl_closure_ordering(ordering, ring_data)
    if isnothing(ordering)
        elim_vars = Any[v for v in ring_data.aux_generators]
        if haskey(ring_data.module_var_lookup, "edt")
            push!(elim_vars, ring_data.module_var_lookup["edt"])
        end

        tail_vars = Any[v for v in ring_data.polynomial_generators]
        append!(
            tail_vars,
            [v for (name, v) in zip(ring_data.module_var_names, ring_data.module_generators) if name != "edt"],
        )
        tail_order = ParamPunPam.DegRevLex(tail_vars...)
        isempty(elim_vars) && return tail_order
        return ParamPunPam.Lex(elim_vars...) * tail_order
    elseif ordering isa String
        lookup = merge(ring_data.aux_var_lookup, ring_data.polynomial_var_lookup, ring_data.module_var_lookup)
        return _module_groebner_ordering_from_string(ordering, lookup)
    else
        return ordering
    end
end

function _prolong_generators(gens::Vector{<:OrePoly}, A0::OreAlg, trunc_order::Int)
    ndiff = length(A0.inp.ratdiffvars[2])
    one_coeff = one(ctx(A0))
    prolonged = OrePoly{eltype_co(A0), eltype_mo(A0)}[]

    for g in gens
        ord_g = _operator_order(g, ndiff)
        ord_g <= trunc_order || error("generator of order $(ord_g) exceeds truncation order $(trunc_order)")
        for beta in _module_basis_exponents(ndiff, trunc_order - ord_g)
            if all(iszero, beta)
                push!(prolonged, copy(g))
            else
                push!(prolonged, mul(makepoly(one_coeff, _source_diff_monomial(beta, A0)), g, A0))
            end
        end
    end
    return prolonged
end

function _operator_order(g::OrePoly, ndiff::Int)
    max_order = 0
    for i in 1:length(g)
        cur = 0
        for j in 1:ndiff
            cur += Int(mon(g, i)[j])
        end
        max_order = max(max_order, cur)
    end
    return max_order
end

function _source_diff_monomial(exponents::NTuple, A0::OreAlg)
    E = exptype(eltype_mo(A0))
    full = zeros(E, nvars(A0))
    for i in eachindex(exponents)
        full[i] = E(exponents[i])
    end
    return makemon(full)
end

function _convert_operator_to_module_polynomial(
    op::OrePoly,
    A0::OreAlg,
    target_ring,
    source_coeff_substitution::Vector,
    module_var_map::Dict,
    trunc_order::Int
)
    ndiff = length(A0.inp.ratdiffvars[2])
    acc = zero(target_ring)
    for i in 1:length(op)
        key = _module_key_from_source_monomial(mon(op, i), ndiff, A0)
        total_deg = sum(key)
        total_deg <= trunc_order || error("operator term of order $(total_deg) exceeds truncation order $(trunc_order)")
        coeff_target = _coefficient_to_target_ring(coeff(op, i), target_ring, source_coeff_substitution)
        acc += coeff_target * module_var_map[key]
    end
    return acc
end

function _module_key_from_source_monomial(m::OreMonVE, ndiff::Int, A0::OreAlg)
    for i in (ndiff + 1):nvars(A0)
        iszero(m[i]) || error("source algebra monomials must only involve differential operators")
    end
    return ntuple(i -> Int(m[i]), ndiff)
end

function _source_coefficient_names(A0::OreAlg)
    return vcat(A0.inp.ratdiffvars[1], A0.inp.ratvars)
end

function _build_source_coeff_substitution(A0::OreAlg, symbol_lookup::Dict{String, T}) where T
    names = _source_coefficient_names(A0)
    return T[
        get(symbol_lookup, name) do
            error("missing image for source coefficient variable $(name)")
        end for name in names
    ]
end

function _coefficient_to_target_ring(c::RatFun, target_ring, source_coeff_substitution::Vector)
    num = _commutative_poly_substitution(Nemo.numerator(c, false), source_coeff_substitution, target_ring)
    den = _commutative_poly_substitution(Nemo.denominator(c, false), source_coeff_substitution, target_ring)
    if isone(den)
        return num
    end
    is_constant(den) || error("coefficient denominator is not in the coefficient field after substitution")
    den_coeff = first(coefficients(den))
    return num * target_ring(inv(den_coeff))
end

function _coefficient_to_target_ring(c, target_ring, ::Vector)
    return target_ring(base_ring(target_ring)(c))
end

function _commutative_poly_substitution(p, source_coeff_substitution::Vector{T}, target_ring) where T
    acc = zero(target_ring)
    @inbounds for (coeff_term, exponents) in zip(coefficients(p), _exponent_vectors_as_vectors(p))
        term = target_ring(base_ring(target_ring)(coeff_term))
        for i in eachindex(exponents)
            expo = exponents[i]
            expo == 0 && continue
            term *= source_coeff_substitution[i]^expo
        end
        acc += term
    end
    return acc
end

function _exponent_vectors_as_vectors(p)
    try
        return exponent_vectors(Vector{Int}, p)
    catch err
        err isa MethodError || rethrow()
        return (_normalize_exponent_vector(exponents) for exponents in exponent_vectors(p))
    end
end

_normalize_exponent_vector(exponent::Integer) = [Int(exponent)]
_normalize_exponent_vector(exponents) = [Int(exponent) for exponent in exponents]

function _module_source_localizing_factor(::Nothing, ::Vector{<:OrePoly}, ::OreAlg)
    return nothing
end

function _module_source_localizing_factor(factor::Symbol, gens::Vector{<:OrePoly}, A0::OreAlg)
    if factor == :none
        return nothing
    elseif factor == :singular_locus
        sl = singular_locus(gens, A0, singular_locus_init(A0))
        isone(sl) && return sl
        return prod(fact[1] for fact in Nemo.factor(sl))
    end
    error("unsupported localizing_factor symbol $(factor)")
end

function _module_source_localizing_factor(factor::String, ::Vector{<:OrePoly}, A0::OreAlg)
    return _module_constant_operator_coefficient(parse_OrePoly(factor, A0), A0)
end

function _module_source_localizing_factor(factor::OrePoly, ::Vector{<:OrePoly}, A0::OreAlg)
    return _module_constant_operator_coefficient(factor, A0)
end

function _module_source_localizing_factor(factor, ::Vector{<:OrePoly}, ::OreAlg)
    return factor
end

function _module_constant_operator_coefficient(p::OrePoly, A::OreAlg)
    iszero(p, A) && error("localizing factor must be non-zero")
    one_mon = makemon(-1, A)
    for i in 1:length(p)
        mon(p, i) == one_mon || error("localizing factor must be a coefficient polynomial")
    end
    return coeff(p, 1)
end

_module_is_trivial_localizing_factor(::Nothing) = true
_module_is_trivial_localizing_factor(factor) = isone(factor)

function _module_localizing_factor_to_target_ring(factor::RatFun, A0::OreAlg, target_ring, source_coeff_substitution::Vector)
    return _coefficient_to_target_ring(factor, target_ring, source_coeff_substitution)
end

function _module_localizing_factor_to_target_ring(factor, A0::OreAlg, target_ring, source_coeff_substitution::Vector)
    return _commutative_poly_substitution(factor, source_coeff_substitution, target_ring)
end

function _max_operator_order(gens::Vector{<:OrePoly}, A0::OreAlg)
    ndiff = length(A0.inp.ratdiffvars[2])
    isempty(gens) && error("at least one generator is required")
    return maximum(_operator_order(g, ndiff) for g in gens)
end

function _module_default_weyl_target_algebra(
    A0::OreAlg,
    ratvars::Vector{String},
    polvars::Vector{String};
    target_order::Union{Nothing, String} = nothing,
)
    diffops = copy(A0.inp.ratdiffvars[2])
    order = isnothing(target_order) ? _module_default_weyl_order(polvars, diffops) : target_order
    return OreAlg(
        order = order,
        ratvars = copy(ratvars),
        poldiffvars = (copy(polvars), diffops),
        varord = A0.varord,
    )
end

function _module_default_weyl_order(polvars::Vector{String}, diffops::Vector{String})
    blocks = String[]
    !isempty(polvars) && push!(blocks, "grevlex " * join(polvars, " "))
    !isempty(diffops) && push!(blocks, "grevlex " * join(diffops, " "))
    isempty(blocks) && error("cannot build a Weyl target algebra without variables")
    return join(blocks, " > ")
end

function _validate_module_weyl_target_algebra(
    A::OreAlg,
    ratvars::Vector{String},
    polvars::Vector{String},
    diffops::Vector{String},
)
    for name in ratvars
        haskey(A.ratvars, name) || error("target algebra is missing parameter $(name)")
    end
    for name in vcat(polvars, diffops)
        haskey(A.strvar_to_indexp, name) || error("target algebra is missing Weyl variable $(name)")
    end
    return nothing
end

function _module_groebner_basis_to_weyl_operators(module_data, target_A::OreAlg)
    operators = OrePoly{eltype_co(target_A), eltype_mo(target_A)}[]
    for p in module_data.basis
        op = _module_relation_to_weyl_operator(p, module_data, target_A)
        isnothing(op) && continue
        iszero(op, target_A) && continue
        push!(operators, op)
    end
    return operators
end

function _module_relation_to_weyl_operator(p, module_data, target_A::OreAlg)
    ring_data = module_data.ring_data
    module_monomials = module_data.module_monomials
    naux = length(ring_data.aux_var_names)
    npol = length(ring_data.polynomial_var_names)
    nmod = length(ring_data.module_var_names)
    res = zero(target_A)

    for (coeff_term, exponents) in zip(coefficients(p), exponent_vectors(Vector{Int}, p))
        for i in 1:naux
            exponents[i] == 0 || return nothing
        end

        module_exponents = view(exponents, (naux + npol + 1):(naux + npol + nmod))
        sum(module_exponents) == 1 || return nothing
        module_index = findfirst(==(1), module_exponents)
        module_index === nothing && return nothing

        coeff_target = _paramgb_coeff_to_ore_coeff(coeff_term, target_A)
        mon_target = _module_term_monomial_to_weyl_monomial(
            exponents,
            ring_data,
            module_monomials[module_index],
            module_data.diffop_names,
            target_A,
        )
        res = add!(res, makepoly(coeff_target, mon_target), target_A)
    end
    return normalize!(res, target_A)
end

function _module_term_monomial_to_weyl_monomial(exponents, ring_data, diff_exponents, diffop_names, target_A::OreAlg)
    E = exptype(eltype_mo(target_A))
    full = zeros(E, nvars(target_A))
    naux = length(ring_data.aux_var_names)

    for (i, name) in enumerate(ring_data.polynomial_var_names)
        expo = exponents[naux + i]
        expo == 0 && continue
        full[target_A.strvar_to_indexp[name]] = E(expo)
    end

    for (diffop, expo) in zip(diffop_names, diff_exponents)
        expo == 0 && continue
        full[target_A.strvar_to_indexp[diffop]] += E(expo)
    end
    return makemon(full)
end

function _paramgb_coeff_to_ore_coeff(c, A::OreAlg)
    ct = ctx(A)
    if ct isa RatFunCtx
        return _paramgb_coeff_to_ratfun(c, ct)
    elseif ct isa QQCtx
        return _paramgb_coeff_to_qq(c)
    end
    error("conversion of ParamPunPam coefficients to $(typeof(ct)) is not implemented")
end

function _paramgb_coeff_to_qq(c::QQFieldElem)
    return c
end

function _paramgb_coeff_to_qq(c)
    num = Nemo.numerator(c, false)
    den = Nemo.denominator(c, false)
    return _constant_paramgb_poly_to_qq(num) / _constant_paramgb_poly_to_qq(den)
end

function _constant_paramgb_poly_to_qq(p::QQFieldElem)
    return p
end

function _constant_paramgb_poly_to_qq(p)
    is_constant(p) || error("coefficient depends on parameters, but target algebra has no ratvars")
    return constant_coefficient(p)
end

function _paramgb_coeff_to_ratfun(c::QQFieldElem, ct::RatFunCtx)
    return _qq_to_ratfun(c, ct)
end

function _paramgb_coeff_to_ratfun(c, ct::RatFunCtx)
    num = Nemo.numerator(c, false)
    den = Nemo.denominator(c, false)
    return _paramgb_parameter_poly_to_ratfun(num, ct) / _paramgb_parameter_poly_to_ratfun(den, ct)
end

function _paramgb_parameter_poly_to_ratfun(p::QQFieldElem, ct::RatFunCtx)
    return _qq_to_ratfun(p, ct)
end

function _paramgb_parameter_poly_to_ratfun(p, ct::RatFunCtx)
    acc = ct.F(0)
    for (coeff_term, exponents) in zip(coefficients(p), _exponent_vectors_as_vectors(p))
        length(exponents) <= length(ct.vars) || error("too many parameters in ParamPunPam coefficient")
        term = _qq_to_ratfun(coeff_term, ct)
        for i in eachindex(exponents)
            expo = exponents[i]
            expo == 0 && continue
            term *= ct.F(ct.vars[i])^expo
        end
        acc += term
    end
    return acc
end

function _qq_to_ratfun(q::QQFieldElem, ct::RatFunCtx)
    return ct.F(Nemo.numerator(q)) / ct.F(Nemo.denominator(q))
end
