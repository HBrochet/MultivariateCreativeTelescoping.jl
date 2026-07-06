struct PFCauchyConfig
    rho::Int
    homogeneous_degree::Int
    use_trivial_syzygies::Bool
    debug_param::PFDebugParam
end

function evaluate_parameter_many(
    cfg::PFCauchyConfig,
    v::Vector{UInt},
    ::OreAlg;
    denisone::Val{B} = Val(false),
) where {B}
    return fill(cfg, length(v))
end

mutable struct PFZeroTrace{M}
    zero_lms::Set{M}
    echelon_lms::Vector{M}
end

PFZeroTrace(::Type{M}) where {M} = PFZeroTrace{M}(Set{M}(), M[])

function _support_signature(tr::PFZeroTrace{M}, A::OreAlg) where {M}
    zeros = collect(tr.zero_lms)
    echs = copy(tr.echelon_lms)
    sort!(zeros, order = order(A))
    sort!(echs, order = order(A))
    return (zeros, echs)
end

struct PFRedData{T,M}
    spol::OrePoly{T,M}
    red_dts::Vector{OrePoly{T,M}}
    basis::Vector{M}
    der_maps::Vector{Dict{M,OrePoly{T,M}}}
end

PFRedData(red::NamedTuple) = PFRedData(red.spol, red.red_dts, red.basis, red.der_maps)

function Base.:(==)(a::PFRedData, b::PFRedData)
    return a.spol == b.spol &&
           a.red_dts == b.red_dts &&
           a.basis == b.basis &&
           a.der_maps == b.der_maps
end

function _support_signature(red::PFRedData, A::OreAlg)
    return (
        _support_signature(red.spol, A),
        _support_signature(red.red_dts, A),
        red.basis,
        [_support_signature(d, A) for d in red.der_maps],
    )
end

function cauchy_interpolation(
    vec::Vector{<:PFRedData},
    randpoints::Vector{Int},
    A::OreAlg;
    prd::Union{fpPolyRingElem,Nothing} = nothing,
)
    isempty(vec) && error("cannot interpolate empty Picard-Fuchs reduction data")
    spol = cauchy_interpolation([red.spol for red in vec], randpoints, A; prd = prd)
    red_dts = cauchy_interpolation([red.red_dts for red in vec], randpoints, A; prd = prd)
    nder = length(vec[1].der_maps)
    der_maps = [
        cauchy_interpolation([red.der_maps[i] for red in vec], randpoints, A; prd = prd)
        for i in 1:nder
    ]
    return PFRedData(spol, red_dts, copy(vec[1].basis), der_maps)
end

function _pf_red_named_tuple(red::PFRedData, A::OreAlg)
    return (
        dops = copy(A.inp.ratdiffvars[2]),
        red_dts = red.red_dts,
        basis = red.basis,
        der_maps = red.der_maps,
        spol = red.spol,
    )
end

function _pf_final_algebra_for_characteristic(A::OreAlg)
    return OreAlg(
        order = _pf_order(A.inp.ratdiffvars[2], String[], String[]),
        char = char(A),
        ratvars = copy(A.inp.ratvars),
        ratdiffvars = deepcopy(A.inp.ratdiffvars),
    )
end

"""
    picard_fuchs_crt(ratfun; parameters = ["t"], rho = 1, ...)

Compute the one-parameter `picard_fuchs` operator with CRT over primes and
Cauchy interpolation in the parameter. The optional trace remembers homogeneous
pre-reduction candidates whose row reduced to zero, then reuses that information
at later evaluation points.
"""
function picard_fuchs_crt(
    ratfun::String;
    parameters::Vector{String} = ["t"],
    ratvars::Vector{String} = String[],
    rho::Integer = 1,
    use_trivial_syzygies::Bool = false,
    tracer::Bool = true,
    crt_comp::Symbol = :medium,
    ci_comp::Symbol = :medium,
    crt_parallel::Bool = false,
    crt_batch_size::Integer = Threads.nthreads(),
    debug_param::PFDebugParam = pf_debug_param(),
)
    length(parameters) == 1 || error("picard_fuchs_crt currently supports exactly one parameter")
    isempty(ratvars) || error("picard_fuchs_crt currently supports no extra rational variables")

    num, den, A, finalA = _pf_parse_and_homogenize(ratfun, parameters, ratvars)
    op = compute_with_CRT(
        _pf_picard_fuchs_modp_cauchy,
        A,
        num,
        den,
        A,
        Int(rho),
        ci_comp,
        use_trivial_syzygies,
        debug_param;
        param = crt_param(tracer = Val(tracer), comp = Val(crt_comp)),
        parallel = crt_parallel,
        batch_size = crt_batch_size,
    )
    return op, finalA
end

function _pf_picard_fuchs_modp_cauchy(
    num::OrePoly,
    den::OrePoly,
    A::OreAlg,
    rho::Integer,
    ci_comp::Symbol,
    use_trivial_syzygies::Bool = false,
    debug_param::PFDebugParam = pf_debug_param();
    tracer::Val{B} = Val(false),
) where {B}
    gb = _pf_gb_jac(den, A)
    _pf_debug_polynomial_family_stats(debug_param, "Picard-Fuchs CRT Jacobian Groebner basis computed", gb, A)
    cfg = PFCauchyConfig(Int(rho), Int(_pf_total_x_degree(den, A)), use_trivial_syzygies, debug_param)
    if B
        red, trace = _pf_compute_with_cauchy_zero_trace(
            _pf_picard_fuchs_red_eval_at_parameter,
            A,
            num,
            den,
            gb,
            cfg;
            comp = ci_comp,
            learn_trace = Val(true),
        )
        op = _pf_find_LDE(_pf_red_named_tuple(red, A), num, den, A, _pf_final_algebra_for_characteristic(A); debug_param = debug_param)
        return op, trace
    end

    red = _pf_compute_with_cauchy_zero_trace(
        _pf_picard_fuchs_red_eval_at_parameter,
        A,
        num,
        den,
        gb,
        cfg;
        comp = ci_comp,
        learn_trace = Val(false),
    )
    return _pf_find_LDE(_pf_red_named_tuple(red, A), num, den, A, _pf_final_algebra_for_characteristic(A); debug_param = debug_param)
end

function _pf_picard_fuchs_modp_cauchy(
    trace::PFZeroTrace,
    num::OrePoly,
    den::OrePoly,
    A::OreAlg,
    rho::Integer,
    ci_comp::Symbol,
    use_trivial_syzygies::Bool = false,
    debug_param::PFDebugParam = pf_debug_param(),
)
    gb = _pf_gb_jac(den, A)
    _pf_debug_polynomial_family_stats(debug_param, "Picard-Fuchs CRT Jacobian Groebner basis computed", gb, A)
    cfg = PFCauchyConfig(Int(rho), Int(_pf_total_x_degree(den, A)), use_trivial_syzygies, debug_param)
    red = _pf_compute_with_cauchy_zero_trace(
        _pf_picard_fuchs_red_eval_at_parameter,
        A,
        num,
        den,
        gb,
        cfg;
        comp = ci_comp,
        trace = trace,
    )
    return _pf_find_LDE(_pf_red_named_tuple(red, A), num, den, A, _pf_final_algebra_for_characteristic(A); debug_param = debug_param)
end

function _pf_compute_with_cauchy_zero_trace(
    f::F,
    A::OreAlg,
    args...;
    comp::Symbol = :medium,
    learn_trace::Val{B} = Val(false),
    trace = nothing,
) where {F<:Function,B}
    randpoints = Int[]
    ev_res = PFRedData[]
    traces = PFZeroTrace[]
    npoints = 0
    bound = 1
    bnd = _ci_bnd_init(comp, bound)
    succeeded = false
    local prev_res

    glen = guess_length(args...)
    nargs = length(args)
    nA = evaluate_parameter_algebra(1, A)
    vpoints, vargs = evaluate_parameter_many(glen, randpoints, nA, args...)
    vctr = 1
    target_initial = (trace === nothing && B) ? 3 : 1

    while npoints < target_initial
        if vctr > glen
            vpoints, vargs = evaluate_parameter_many(glen, randpoints, nA, args...)
            vctr = 1
        end
        ev_args = _ci_eval_args_at(vctr, vargs, nargs)
        if trace === nothing && B
            tmp, tr = f(nA, ev_args...; tracer = Val(true))
            push!(traces, tr)
        elseif trace === nothing
            tmp = f(nA, ev_args...)
        else
            tmp = f(nA, trace, ev_args...)
        end
        push!(ev_res, tmp)
        push!(randpoints, vpoints[vctr])
        vctr += 1
        npoints += 1
        globalstats.counters[:number_evaluation] += 1
    end

    learned_trace = trace
    if trace === nothing && B
        _pf_majority_trace!(traces, ev_res, randpoints, A)
        length(traces) >= 2 || error("too many bad evaluation points for Picard-Fuchs trace")
        learned_trace = traces[1]
        npoints = length(randpoints)
        while bnd <= npoints
            bound += 1
            bnd = _ci_bnd_next(comp, bound, bnd)
        end
    end

    t = ctx(A).vars[1]
    while true
        if vctr > glen
            vpoints, vargs = evaluate_parameter_many(glen, randpoints, nA, args...)
            vctr = 1
        end
        ev_args = _ci_eval_args_at(vctr, vargs, nargs)
        tmp = learned_trace === nothing ? f(nA, ev_args...) : f(nA, learned_trace, ev_args...)
        push!(ev_res, tmp)
        push!(randpoints, vpoints[vctr])
        vctr += 1
        npoints += 1
        globalstats.counters[:number_evaluation] += 1

        if succeeded
            majority_test!(ev_res, randpoints, A)
            prd = _recompute_prd(randpoints, t, A)
            res = cauchy_interpolation(ev_res, randpoints, A; prd = prd)
            if res == prev_res
                return trace === nothing && B ? (res, learned_trace) : res
            end
            prev_res = res
            succeeded = false
        elseif npoints == bnd
            try
                majority_test!(ev_res, randpoints, A)
                prd = _recompute_prd(randpoints, t, A)
                prev_res = cauchy_interpolation(ev_res, randpoints, A; prd = prd)
                succeeded = true
            catch
            end
            bound += 1
            bnd = _ci_bnd_next(comp, bound, bnd)
        end
    end
end

function _pf_majority_trace!(
    traces::Vector{<:PFZeroTrace},
    ev_res::Vector{<:PFRedData},
    randpoints::Vector{Int},
    A::OreAlg,
)
    sigs = [_support_signature(tr, A) for tr in traces]
    counts = Dict{typeof(sigs[1]),Int}()
    for sig in sigs
        counts[sig] = get(counts, sig, 0) + 1
    end
    best = sigs[1]
    best_count = 0
    for (sig, count) in counts
        if count > best_count
            best = sig
            best_count = count
        end
    end
    best_count * 2 > length(traces) || error("too many bad evaluation points for Picard-Fuchs trace")

    new_traces = eltype(traces)[]
    new_res = eltype(ev_res)[]
    new_points = Int[]
    for i in eachindex(sigs)
        if sigs[i] == best
            push!(new_traces, traces[i])
            push!(new_res, ev_res[i])
            push!(new_points, randpoints[i])
        end
    end
    empty!(traces); append!(traces, new_traces)
    empty!(ev_res); append!(ev_res, new_res)
    empty!(randpoints); append!(randpoints, new_points)
    return
end

function _pf_picard_fuchs_red_eval_at_parameter(
    A::OreAlg,
    num::OrePoly,
    den::OrePoly,
    gb::Vector{<:OrePoly},
    cfg::PFCauchyConfig;
    tracer::Val{B} = Val(false),
) where {B}
    if B
        red, trace = _pf_der_red_map_many_zero_trace(
            gb,
            num,
            den,
            A;
            rho = cfg.rho,
            homogeneous_degree = cfg.homogeneous_degree,
            use_trivial_syzygies = cfg.use_trivial_syzygies,
            debug_param = cfg.debug_param,
            tracer = Val(true),
        )
        return PFRedData(red), trace
    end
    red = _pf_der_red_map_many_zero_trace(
        gb,
        num,
        den,
        A;
        rho = cfg.rho,
        homogeneous_degree = cfg.homogeneous_degree,
        use_trivial_syzygies = cfg.use_trivial_syzygies,
        debug_param = cfg.debug_param,
    )
    return PFRedData(red)
end

function _pf_picard_fuchs_red_eval_at_parameter(
    A::OreAlg,
    trace::PFZeroTrace,
    num::OrePoly,
    den::OrePoly,
    gb::Vector{<:OrePoly},
    cfg::PFCauchyConfig,
)
    red = _pf_der_red_map_many_zero_trace(
        gb,
        num,
        den,
        A;
        rho = cfg.rho,
        homogeneous_degree = cfg.homogeneous_degree,
        use_trivial_syzygies = cfg.use_trivial_syzygies,
        debug_param = cfg.debug_param,
        trace = trace,
    )
    return PFRedData(red)
end

function _pf_der_red_map_many_zero_trace(
    gb::Vector{OrePoly{T,M}},
    num::OrePoly,
    den::OrePoly,
    A::OreAlg;
    rho::Integer = 1,
    homogeneous_degree::Integer,
    use_trivial_syzygies::Bool = false,
    debug_param::PFDebugParam = pf_debug_param(),
    tracer::Val{B} = Val(false),
    trace::Union{Nothing,PFZeroTrace{M}} = nothing,
) where {T,M,B}
    geob = GeoBucket(zero(A))
    tmp_poly = ReuseOrePoly(1, A)
    learned = B ? PFZeroTrace(M) : trace
    trivial_lms = use_trivial_syzygies ? _pf_trivial_syzygy_leading_monomials(den, A) : M[]
    spol, g1, red_dts, echelon = _pf_der_red_map_precomp_zero_trace(
        deepcopy(num),
        deepcopy(gb),
        A,
        geob,
        tmp_poly;
        rho = rho,
        homogeneous_degree = homogeneous_degree,
        trivial_lms = trivial_lms,
        debug_param = debug_param,
        trace = learned,
        record_trace = B,
    )
    _pf_debug_polynomial_family_stats(debug_param, "Picard-Fuchs CRT confinement finished", echelon, A)
    if B && learned !== nothing
        empty!(learned.echelon_lms)
        append!(learned.echelon_lms, [lm(p) for p in echelon])
    end
    der_maps, basis = find_der_red_map_many(spol, g1, red_dts, echelon, A, geob, tmp_poly)
    debug(debug_param) && @debug "Picard-Fuchs CRT confinement basis computed" confinement_dim=length(basis) approximate_lde_order=length(basis)
    red = (dops = copy(A.inp.ratdiffvars[2]),
           red_dts = red_dts,
           basis = basis,
           der_maps = der_maps,
           spol = spol)
    return B ? (red, learned) : red
end

function _pf_der_red_map_precomp_zero_trace(
    spol::OrePoly{T,M},
    gb::Vector{OrePoly{T,M}},
    A::OreAlg,
    geob::GeoBucket,
    tmp_poly::ReuseOrePoly;
    rho::Integer = 1,
    homogeneous_degree::Integer,
    trivial_lms::Union{Nothing,Vector{M}} = nothing,
    debug_param::PFDebugParam = pf_debug_param(),
    trace::Union{Nothing,PFZeroTrace{M}} = nothing,
    record_trace::Bool = false,
) where {T,M}
    red_dts = find_red_dts(gb, A)
    Ls = [red_dt_tail(red_dt, A) for red_dt in red_dts]
    hdeg = Int(homogeneous_degree)
    g1, g2 = separate(gb, A)
    tsyzlm = trivial_lms === nothing ? M[] : trivial_lms
    isempty(g2) && return spol, g1, red_dts, OrePoly{T,M}[]

    spol = GD_reduction1!(spol, g1, A, geob, tmp_poly)
    s0 = Int(degree(mon(spol, 1)))
    s = s0
    rho >= 0 || error("rho must be non-negative")

    echelon, frontiers, seen, built_degree = _pf_homogeneous_prereduction_init_zero_trace(
        g2, g1, s0, hdeg, tsyzlm, A, geob, tmp_poly, debug_param, trace, record_trace)
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
                built_degree = _pf_homogeneous_prereduction_advance_to_degree_zero_trace!(
                    echelon, frontiers, g1, g2, seen, built_degree, target_degree,
                    s0, hdeg, tsyzlm, A, geob, tmp_poly, debug_param, trace, record_trace)
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
    _pf_homogeneous_prereduction_advance_to_degree_zero_trace!(
        echelon, frontiers, g1, g2, seen, built_degree, target_degree,
        s0, hdeg, tsyzlm, A, geob, tmp_poly, debug_param, trace, record_trace)
    spol = reduce_with_echelon!(echelon, spol, A, geob, tmp_poly)
    return spol, g1, red_dts, echelon
end

function _pf_homogeneous_prereduction_init_zero_trace(
    g2::Vector{OrePoly{T,M}},
    g1::Vector{OrePoly{T,M}},
    s0::Int,
    homogeneous_degree::Int,
    trivial_lms::Vector{M},
    A::OreAlg,
    geob::GeoBucket,
    tmp_poly::ReuseOrePoly,
    debug_param::PFDebugParam,
    trace::Union{Nothing,PFZeroTrace{M}},
    record_trace::Bool,
) where {T,M}
    echelon = OrePoly{T,M}[]
    frontiers = Dict{Int,Vector{OrePoly{T,M}}}()
    seen = Set{M}()
    built_degree = typemax(Int)
    for g in g2
        current_degree = _pf_homogeneous_prereduction_monomial_degree(lm(g), A, homogeneous_degree)
        built_degree = min(built_degree, current_degree)
        current = get!(() -> OrePoly{T,M}[], frontiers, current_degree)
        _pf_add_homogeneous_prereduction_candidate_zero_trace!(
            echelon, current, seen, g, g1, g2, trivial_lms, A, geob, tmp_poly, trace;
            check_trivial = false,
            add_to_echelon = _pf_is_homogeneous_prereduction_degree(current_degree, s0, homogeneous_degree),
            candidate_degree = current_degree,
            debug_param = debug_param,
            record_trace = record_trace,
        )
    end
    built_degree == typemax(Int) && (built_degree = s0)
    return echelon, frontiers, seen, built_degree
end

function _pf_homogeneous_prereduction_advance_to_degree_zero_trace!(
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
    trace::Union{Nothing,PFZeroTrace{M}},
    record_trace::Bool,
) where {T,M}
    target_degree <= built_degree && return built_degree
    empty_frontier = OrePoly{T,M}[]
    for next_degree in built_degree+1:target_degree
        previous = get(frontiers, next_degree - 1, empty_frontier)
        isempty(previous) && continue
        current = get!(() -> OrePoly{T,M}[], frontiers, next_degree)
        add_to_echelon = _pf_is_homogeneous_prereduction_degree(next_degree, s0, homogeneous_degree)
        for g in previous
            for i in xvars_range(A)
                p = mul(makepoly(one(ctx(A)), makemon(i, A)), g, A)
                _pf_add_homogeneous_prereduction_candidate_zero_trace!(
                    echelon, current, seen, p, g1, g2, trivial_lms, A, geob, tmp_poly, trace;
                    check_trivial = true,
                    add_to_echelon = add_to_echelon,
                    candidate_degree = next_degree,
                    debug_param = debug_param,
                    record_trace = record_trace,
                )
            end
        end
    end
    return target_degree
end

function _pf_add_homogeneous_prereduction_candidate_zero_trace!(
    echelon::Vector{OrePoly{T,M}},
    frontier::Vector{OrePoly{T,M}},
    seen::Set{M},
    p::OrePoly{T,M},
    g1::Vector{OrePoly{T,M}},
    g2::Vector{OrePoly{T,M}},
    trivial_lms::Vector{M},
    A::OreAlg,
    geob::GeoBucket,
    tmp_poly::ReuseOrePoly,
    trace::Union{Nothing,PFZeroTrace{M}};
    check_trivial::Bool,
    add_to_echelon::Bool,
    candidate_degree::Int = -1,
    debug_param::PFDebugParam = pf_debug_param(),
    record_trace::Bool = false,
) where {T,M}
    isempty(p) && return
    m = lm(p)
    m in seen && return
    push!(seen, m)
    if check_trivial && _pf_is_trivial_syzygy_leading_monomial(m, trivial_lms, A)
        return
    end
    if add_to_echelon
        if trace !== nothing && m in trace.zero_lms
            push!(frontier, p)
            return
        end
        red = GD_reduction1!(p, g1, A, geob, tmp_poly)
        red = reduce_with_echelon!(echelon, red, A, geob, tmp_poly)
        if length(red) == 0
            record_trace && trace !== nothing && push!(trace.zero_lms, m)
        else
            add_echelon!(echelon, red, A)
            debug(debug_param) && @debug "Picard-Fuchs CRT echelon element added" candidate_degree=candidate_degree max_degree=_pf_poly_max_degree(red) max_x_degree=_pf_poly_max_x_degree(red, A) monomials=length(red) echelon_size=length(echelon)
        end
    end
    push!(frontier, p)
    return
end
