"""
        picard_fuchs(ratfun; parameters = ["t"], ratvars = String[], rho = 1, use_trivial_syzygies = false, debug_param = pf_debug_param())

Compute a Picard-Fuchs operator for a rational integrand by the Griffiths-Dwork-Lairez reduction.
The input must be a parsable rational function. The unknown in parameters avec variables wrt 
one wants to compute differential variables, variables in ratvars are variables of the base field,
and the remaining variables are interpreted as the integration variables.
For homogeneous reductions, `rho` counts extra homogeneous strata computed after confinement.
Set `use_trivial_syzygies=true` to prune homogeneous pre-reduction candidates using
the leading monomial ideal of the trivial Koszul syzygies. It is disabled by default.
"""
struct PFDebugParam{D} end
pf_debug_param(;debug::Val{D} = Val(false)) where {D} = PFDebugParam{D}()
debug(::PFDebugParam{D}) where {D} = D

function picard_fuchs(
    ratfun :: String;
    parameters :: Vector{String} = ["t"],
    ratvars :: Vector{String} = String[],
    rho::Integer = 1,
    use_trivial_syzygies::Bool = false,
    debug_param::PFDebugParam = pf_debug_param())

    num, den, A, finalA = _pf_parse_and_homogenize(ratfun, parameters, ratvars)
    return _picard_fuchs(num, den, A, finalA; rho = rho, use_trivial_syzygies = use_trivial_syzygies, debug_param = debug_param)
end

# main internal function
function _picard_fuchs(
    num :: OrePoly,
    den :: OrePoly,
    A :: OreAlg,
    finalA :: OreAlg;
    rho::Integer = 1,
    use_trivial_syzygies::Bool = false,
    debug_param::PFDebugParam = pf_debug_param(),
)
    gb = _pf_gb_jac(den,A) 
    _pf_debug_polynomial_family_stats(debug_param, "Picard-Fuchs Jacobian Groebner basis computed", gb, A)
    red = _pf_der_red_map_many(gb, num, den, A; rho = rho, use_trivial_syzygies = use_trivial_syzygies, debug_param = debug_param)
    lde = _pf_find_LDE(red, num, den, A, finalA; debug_param = debug_param) 
    return lde, finalA
end






function _pf_parse_and_homogenize(ratfun :: String, parameters :: Vector{String}, ratvars :: Vector{String})
    num, den, affA = _pf_parse_rational(ratfun, parameters, ratvars)

    inp = deepcopy(affA.inp)
    inp.poldiffvars = (vcat(["x0"], inp.poldiffvars[1]), vcat(["dx0"], inp.poldiffvars[2]))
    inp.order = _pf_order(inp.ratdiffvars[2], inp.poldiffvars[1], inp.poldiffvars[2])
    A = OreAlg(inp)

    num = map_algebras(num, affA, A)
    den = map_algebras(den, affA, A)
    hdata = _pf_homogenize_fraction(num, den, A, length(affA.inp.poldiffvars[1]), 1; force_x0_factor = true)
    hnum = hdata.numerator
    hden = hdata.f
    finalA = OreAlg(order = _pf_order(A.inp.ratdiffvars[2], String[], String[]),
                    ratvars = copy(ratvars),
                    ratdiffvars = deepcopy(A.inp.ratdiffvars))
    _pf_forbid_derivative_multipliers!(A)
    return hnum, hden, A, finalA
end

function _pf_parse_rational(ratfun :: String, parameters :: Vector{String}, ratvars :: Vector{String})
    expr = Meta.parse(ratfun, raise = true)
    vars = _infer_ratdiffvars(expr)
    parameter_set = Set(parameters)
    ratvar_set = Set(ratvars)
    pvars = [v for v in vars if !(v in parameter_set) && !(v in ratvar_set)]
    pdops = ["d" * v for v in pvars]
    rdops = ["d" * v for v in parameters]

    A = OreAlg(order = _pf_order(rdops, pvars, pdops),
               ratvars = copy(ratvars),
               ratdiffvars = (copy(parameters), rdops),
               poldiffvars = (pvars, pdops))
    parseA = OreAlg(order = _pf_order(vcat(rdops, pdops), String[], String[]),
                    ratvars = copy(ratvars),
                    ratdiffvars = (vcat(parameters, pvars), vcat(rdops, pdops)))
    rat = _dfinite_expr_to_ratfun(expr, parseA)
    num = _pf_coeff_polynomial_to_orepoly(Nemo.numerator(rat, false), parseA, A, pvars, parameters)
    den = _pf_coeff_polynomial_to_orepoly(Nemo.denominator(rat, false), parseA, A, pvars, parameters)
    return num, den, A
end

function _pf_coeff_polynomial_to_orepoly(p, parseA::OreAlg, A::OreAlg, variables::Vector{String}, parameters::Vector{String})
    source_names = vcat(parseA.inp.ratdiffvars[1], parseA.inp.ratvars)
    variable_set = Set(variables)
    parameter_set = Set(parameters)
    res = zero(A)

    for (coeff_term, exponents) in zip(coefficients(p), exponent_vectors(Vector{Int}, p))
        coeff_target = convertn(coeff_term, ctx(A))
        mon_target = makemon(-1, A)
        for (name, expo) in zip(source_names, exponents)
            expo == 0 && continue
            if name in variable_set
                mon_target *= makemon(A.strvar_to_indexp[name], A)^expo
            elseif name in parameter_set || haskey(A.ratvars, name)
                coeff_target *= A.ratvars[name]^expo
            else
                error("cannot map variable $(name) to Picard-Fuchs algebra")
            end
        end
        res = add!(res, makepoly(coeff_target, mon_target), A)
    end
    return res
end

function _pf_order(rdops::Vector{String}, pvars::Vector{String}, pdops::Vector{String})
    blocks = String[]
    isempty(rdops) || push!(blocks, "lex " * join(rdops, " "))
    isempty(pvars) || push!(blocks, "grevlex " * join(pvars, " "))
    isempty(pdops) || push!(blocks, "grevlex " * join(pdops, " "))
    return join(blocks, " > ")
end

function _pf_gb_jac(den :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    gens = OrePoly{T,M}[]
    # construct the diff(den,xi) - dxi
    for i in A.nrdv+1:A.nrdv+A.npdv 
        push!(gens,sub(diff(den,i,A), makepoly(one(ctx(A)), makemon(i,A)),A))
    end 
    gb = f4(gens,A) 
    # add the the dti - diff(den,ti)
    for i in 1:A.nrdv
        push!(gb, sub(makepoly(one(ctx(A)), makemon(i,A)), diff(den,i,A),A))
    end
    return gb 
end 

function _pf_poly_max_degree(p::OrePoly)
    isempty(p) && return -1
    return maximum(Int(degree(m)) for m in mons(p))
end

function _pf_poly_max_x_degree(p::OrePoly, A::OreAlg)
    isempty(p) && return -1
    return maximum(_pf_x_degree(m, A) for m in mons(p))
end

function _pf_polynomial_family_stats(polys::Vector{<:OrePoly}, A::OreAlg)
    max_degree = -1
    max_x_degree = -1
    total_monomials = 0
    for p in polys
        total_monomials += length(p)
        max_degree = max(max_degree, _pf_poly_max_degree(p))
        max_x_degree = max(max_x_degree, _pf_poly_max_x_degree(p, A))
    end
    return (count = length(polys),
            max_degree = max_degree,
            max_x_degree = max_x_degree,
            total_monomials = total_monomials)
end

function _pf_debug_polynomial_family_stats(param::PFDebugParam, label::String, polys::Vector{<:OrePoly}, A::OreAlg)
    debug(param) || return nothing
    stats = _pf_polynomial_family_stats(polys, A)
    @debug label elements=stats.count max_degree=stats.max_degree max_x_degree=stats.max_x_degree total_monomials=stats.total_monomials
    return nothing
end

# todo: utiliser l'homogénéité pour réduire le cout de l'algèbre linéaire.
function _pf_der_red_map_many(
    gb :: Vector{OrePoly{T,M}},
    num :: OrePoly,
    den :: OrePoly,
    A :: OreAlg;
    rho::Integer = 1,
    use_trivial_syzygies::Bool = false,
    debug_param::PFDebugParam = pf_debug_param(),
) where {T,M}
    geob = GeoBucket(zero(A))
    tmp_poly = ReuseOrePoly(1, A)
    trivial_lms = use_trivial_syzygies ? _pf_trivial_syzygy_leading_monomials(den, A) : M[]
    spol, g1, red_dts, echelon = _pf_der_red_map_precomp(
        deepcopy(num),
        deepcopy(gb),
        A,
        geob,
        tmp_poly;
        rho = rho,
        homogeneous_degree = _pf_total_x_degree(den, A),
        trivial_lms = trivial_lms,
        debug_param = debug_param,
    )
    _pf_debug_polynomial_family_stats(debug_param, "Picard-Fuchs confinement finished", echelon, A)
    der_maps, basis = find_der_red_map_many(spol, g1, red_dts, echelon, A, geob, tmp_poly)
    debug(debug_param) && @debug "Picard-Fuchs confinement basis computed" confinement_dim=length(basis) approximate_lde_order=length(basis)
    return (dops = copy(A.inp.ratdiffvars[2]),
            red_dts = red_dts,
            basis = basis,
            der_maps = der_maps,
            der_mats = der_maps_to_matrices(der_maps, basis, A),
            spol = spol)
end

function _pf_find_LDE(red :: NamedTuple, num :: OrePoly, den :: OrePoly, A :: OreAlg, finalA :: OreAlg; debug_param::PFDebugParam = pf_debug_param())
    debug(debug_param) && @debug "Picard-Fuchs LDE computation starts" confinement_dim=length(red.basis) approximate_lde_order=length(red.basis)
    return map_algebras(find_LDE_direct_many(red, A), A, finalA)
end

function _pf_der_red_map_precomp(
    spol :: OrePoly{T,M},
    gb :: Vector{OrePoly{T,M}},
    A :: OreAlg,
    geob :: GeoBucket,
    tmp_poly :: ReuseOrePoly;
    rho::Integer = 1,
    homogeneous_degree::Integer,
    trivial_lms::Union{Nothing,Vector{M}} = nothing,
    debug_param::PFDebugParam = pf_debug_param(),
) where {T,M}
    red_dts = find_red_dts(gb, A)
    Ls = [red_dt_tail(red_dt, A) for red_dt in red_dts]
    hdeg = Int(homogeneous_degree)
    g1, g2 = separate(gb, A)
    tsyzlm = trivial_lms === nothing ? M[] : trivial_lms

    if isempty(g2)
        return spol, g1, red_dts, OrePoly{T,M}[]
    end

    spol = GD_reduction1!(spol, g1, A, geob, tmp_poly)
    s0 = Int(degree(mon(spol,1)))
    s = s0
    rho >= 0 || error("rho must be non-negative")

    echelon, frontiers, seen, built_degree = _pf_homogeneous_prereduction_init(
        g2,
        g1,
        s0,
        hdeg,
        tsyzlm,
        A,
        geob,
        tmp_poly,
        debug_param,
    )
    spol = reduce_with_echelon!(echelon, spol, A, geob, tmp_poly)

    while true
        done = SortedSet{M}(order(A))
        todo = SortedSet{M}(order(A))
        append!(todo, mons(spol))
        restart = false
        while !isempty(todo)
            m = poplast!(todo)
            if _pf_x_degree(m, A) > s
                s += 1
                target_degree = _pf_next_homogeneous_prereduction_degree(s, s0, hdeg)
                built_degree = _pf_homogeneous_prereduction_advance_to_degree!(
                    echelon,
                    frontiers,
                    g1,
                    g2,
                    seen,
                    built_degree,
                    target_degree,
                    s0,
                    hdeg,
                    tsyzlm,
                    A,
                    geob,
                    tmp_poly,
                    debug_param,
                )
                spol = reduce_with_echelon!(echelon, spol, A, geob, tmp_poly)
                restart = true
                break
            end
            push!(done, m)
            for L in Ls
                isempty(L) && continue
                geob = addmul_geobucket!(geob, one(ctx(A)), m, L, A)
                Lm = GD_reduction1!(g1, A, geob, tmp_poly)
                Lm = reduce_with_echelon!(echelon, Lm, A, geob, tmp_poly)
                for mo in mons(Lm)
                    mo in done || push!(todo, mo)
                end
            end
        end
        restart && continue
        break 
    end
    target_degree = built_degree + Int(rho) * hdeg
    built_degree = _pf_homogeneous_prereduction_advance_to_degree!(
                    echelon,
                    frontiers,
                    g1,
                    g2,
                    seen,
                    built_degree,
                    target_degree,
                    s0,
                    hdeg,
                    tsyzlm,
                    A,
                    geob,
                    tmp_poly,
                    debug_param,
                )
    spol = reduce_with_echelon!(echelon, spol, A, geob, tmp_poly)
    return spol, g1, red_dts, echelon
end

function _pf_homogeneous_prereduction_init(
    g2::Vector{OrePoly{T,M}},
    g1::Vector{OrePoly{T,M}},
    s0::Int,
    homogeneous_degree::Int,
    trivial_lms::Vector{M},
    A::OreAlg,
    geob::GeoBucket,
    tmp_poly::ReuseOrePoly,
    debug_param::PFDebugParam,
) where {T,M}
    echelon = OrePoly{T,M}[]
    frontiers = Dict{Int,Vector{OrePoly{T,M}}}()
    seen = Set{M}()
    built_degree = typemax(Int)

    for g in g2
        current_degree = _pf_homogeneous_prereduction_monomial_degree(lm(g), A, homogeneous_degree)
        built_degree = min(built_degree, current_degree)
        current = get!(frontiers, current_degree) do
            OrePoly{T,M}[]
        end
        _pf_add_homogeneous_prereduction_candidate!(
            echelon,
            current,
            seen,
            g,
            g1,
            g2,
            trivial_lms,
            A,
            geob,
            tmp_poly;
            check_trivial = false,
            add_to_echelon = _pf_is_homogeneous_prereduction_degree(current_degree, s0, homogeneous_degree),
            candidate_degree = current_degree,
            debug_param = debug_param,
        )
    end
    built_degree == typemax(Int) && (built_degree = s0)
    return echelon, frontiers, seen, built_degree
end

function _pf_homogeneous_prereduction_advance_to_degree!(
    echelon::Vector{OrePoly{T,M}},
    frontiers::Dict{Int,Vector{OrePoly{T,M}}},
    g1::Vector{OrePoly{T,M}},
    g2::Vector{OrePoly{T,M}},
    seen::Set{M},
    built_degree::Int,
    target_degree::Int,
    s0::Int,
    homogeneous_degree::Int,
    trivial_lms::Vector{M},
    A::OreAlg,
    geob::GeoBucket,
    tmp_poly::ReuseOrePoly,
    debug_param::PFDebugParam,
) where {T,M}
    target_degree <= built_degree && return built_degree
    empty_frontier = OrePoly{T,M}[]
    for next_degree in built_degree+1:target_degree
        previous = get(frontiers, next_degree - 1, empty_frontier)
        isempty(previous) && continue
        current = get!(frontiers, next_degree) do
            OrePoly{T,M}[]
        end
        add_to_echelon = _pf_is_homogeneous_prereduction_degree(next_degree, s0, homogeneous_degree)
        for g in previous
            for i in xvars_range(A)
                p = mul(makepoly(one(ctx(A)), makemon(i, A)), g, A)
                _pf_add_homogeneous_prereduction_candidate!(
                    echelon,
                    current,
                    seen,
                    p,
                    g1,
                    g2,
                    trivial_lms,
                    A,
                    geob,
                    tmp_poly;
                    check_trivial = true,
                    add_to_echelon = add_to_echelon,
                    candidate_degree = next_degree,
                    debug_param = debug_param,
                )
            end
        end
    end
    return target_degree
end

function _pf_add_homogeneous_prereduction_candidate!(
    echelon::Vector{OrePoly{T,M}},
    frontier::Vector{OrePoly{T,M}},
    seen::Set{M},
    p::OrePoly{T,M},
    g1::Vector{OrePoly{T,M}},
    g2::Vector{OrePoly{T,M}},
    trivial_lms::Vector{M},
    A::OreAlg,
    geob::GeoBucket,
    tmp_poly::ReuseOrePoly;
    check_trivial::Bool,
    add_to_echelon::Bool,
    candidate_degree::Int = -1,
    debug_param::PFDebugParam = pf_debug_param(),
) where {T,M}
    isempty(p) && return
    m = lm(p)
    m in seen && return
    push!(seen, m)

    # Trivial syzygies are detected by membership in their leading monomial ideal.
    if check_trivial && _pf_is_trivial_syzygy_leading_monomial(m, trivial_lms, A)
        return
    end

    if add_to_echelon
        red = GD_reduction1!(p, g1, A, geob, tmp_poly)
        red = reduce_with_echelon!(echelon, red, A, geob, tmp_poly)
        if length(red) > 0
            add_echelon!(echelon, red, A)
            debug(debug_param) && @debug "Picard-Fuchs echelon element added" candidate_degree=candidate_degree max_degree=_pf_poly_max_degree(red) max_x_degree=_pf_poly_max_x_degree(red, A) monomials=length(red) echelon_size=length(echelon)
        end
    end
    push!(frontier, p)
    return
end

@inline function _pf_is_homogeneous_prereduction_degree(deg::Int, s0::Int, homogeneous_degree::Int)
    return deg >= s0 && (deg - s0) % homogeneous_degree == 0
end

@inline function _pf_next_homogeneous_prereduction_degree(deg::Int, s0::Int, homogeneous_degree::Int)
    deg <= s0 && return s0
    return s0 + cld(deg - s0, homogeneous_degree) * homogeneous_degree
end

@inline function _pf_homogeneous_prereduction_monomial_degree(m::OreMonVE, A::OreAlg, homogeneous_degree::Int)
    deg = _pf_x_degree(m, A)
    dx_weight = homogeneous_degree - 1
    if dx_weight != 0
        for i in dxvars_range(A)
            deg += dx_weight * Int(m[i])
        end
    end
    return deg
end

function _pf_trivial_syzygy_leading_monomials(den::OrePoly{T,M}, A::OreAlg) where {T,M}
    dxs = collect(dxvars_range(A))
    dfs = [diff(den, i, A) for i in dxs]
    dxpolys = [makepoly(one(ctx(A)), makemon(i, A)) for i in dxs]
    lms = M[]
    for i in eachindex(dxs), j in eachindex(dxs)
        i == j && continue
        syz = sub(mul(dfs[i], dxpolys[j], A), mul(dfs[j], dxpolys[i], A), A)
        isempty(syz) || push!(lms, lm(syz))
    end
    return _pf_minimal_monomial_generators(lms, A)
end

function _pf_minimal_monomial_generators(candidates::Vector{M}, A::OreAlg) where {M}
    minimal = M[]
    for m in candidates
        any(g -> divide(g, m), minimal) && continue
        filter!(g -> !divide(m, g), minimal)
        push!(minimal, m)
    end
    return minimal
end

function _pf_is_trivial_syzygy_leading_monomial(m::M, trivial_lms::Vector{M}, A::OreAlg) where {M}
    for lmtriv in trivial_lms
        divide(lmtriv, m) && return true
    end
    return false
end
