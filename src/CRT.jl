struct CRTParam{A,B} end 

@inline function _modp_args(ctxA::RatFunCtx, prime::Int, A::OreAlg, args...)
    nA = change_alg_char_ratfun(prime, A)
    R = Native.GF(prime)
    return nA, change_coefficient_field(R, nA, args...)
end

@inline function _modp_args(ctxA::QQCtx, prime::Int, A::OreAlg, args...)
    nA = change_alg_char_QQ(prime, A)
    return nA, change_coefficient_field(nA, args...)
end

function crt_param(;tracer :: Val{A} = Val(false),
                    comp :: Val{B} = Val(:medium)
                    ) where {A,B}
    return CRTParam{A,B}()
end

tracer(p :: CRTParam{A,B}) where {A,B} = A 
comp(p :: CRTParam{A,B}) where {A,B} = B


@inline _support_signature(p::OrePoly, ::OreAlg) = mons(p)
@inline _support_signature(v::Vector{OrePoly{T,M}}, ::OreAlg) where {T,M}= [mons(p) for p in v]
@inline _support_signature(tr::F4Trace, ::OreAlg) = tr.spairs
@inline function _support_signature(d::Dict{M,OrePoly{T,M}}, A::OreAlg) where {T,M}
    ks = collect(keys(d))
    sort!(ks, order = order(A))
    return [(k, mons(d[k])) for k in ks]
end
@inline _support_signature(t::Tuple{Dict{M,OrePoly{T,M}},OrePoly{T,M}}, A::OreAlg) where {T,M} = (_support_signature(t[1], A), _support_signature(t[2], A))


function majority_test!(res_ev::Vector{T}, primes_::Vector{Int}, A::OreAlg) where {T}
    n = length(res_ev)
    if n < 3
        error("not enough primes for the majority test")
    end
    first_sig = _support_signature(res_ev[1], A)
    counts = Dict{typeof(first_sig),Int}()
    counts[first_sig] = 1
    for i in 2:n
        sig = _support_signature(res_ev[i], A)
        counts[sig] = get(counts, sig, 0) + 1
    end
    maxsig = length(res_ev)
    maxcount = 0
    for (sig, c) in counts
        if c > maxcount
            maxsig = sig
            maxcount = c
        end
    end
    if maxcount * 2 <= n
        error("too many bad primes for the majority test")
    end
    keep = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        keep[i] = _support_signature(res_ev[i], A) == maxsig
    end
    if all(keep)
        return length(primes_)
    end
    new_res = eltype(res_ev)[]
    new_primes = Int[]
    @inbounds for i in 1:n
        if keep[i]
            push!(new_res, res_ev[i])
            push!(new_primes, primes_[i])
        end
    end
    empty!(res_ev)
    append!(res_ev, new_res)
    empty!(primes_)
    append!(primes_, new_primes)
    return length(primes_)
end


# the tracer assume f can be called with f(...;tracer=Val(true)) to learn and return a trace  
# and f(trace,...) to apply the trace at subsequent computations

function compute_with_CRT(f::F, A::OreAlg, args...; param::CRTParam = crt_param()) where {F<:Function}
    nprime = 1
    globalstats.counters[:number_primes] += 1
    primes_ = Int[primes[1]]
    prime = primes[1]
    bound = comp(param) == :slow ? 3 : 2

    bnd = 3
    succeeded = false # after reconstructing the result we try one more prime 

    ctxA = ctx(A)
    if tracer(param)
        nA, mod_args = _modp_args(ctxA, prime, A, args...)
        tmp, trace = f(mod_args...;tracer=Val(true))
        res_modp = [tmp]
        tab_tr = [trace]
        for _ in 2:3
            nprime += 1
            globalstats.counters[:number_primes] += 1
            prime = primes[nprime]
            push!(primes_, prime)
            @debug "computing the function for $(nprime)th prime"
            nA, mod_args = _modp_args(ctxA, prime, A, args...)
            tmp, trace = f(mod_args...;tracer=Val(true))
            push!(res_modp, tmp)
            push!(tab_tr, trace)
        end
        majority_test!(tab_tr,[1,2,3], A)
        if length(tab_tr) < 2 
            error("too many bad primes")
        end
        trace = tab_tr[1]
    else 
        nA, mod_args = _modp_args(ctxA, prime, A, args...)
        res_modp = [f(mod_args...)]
    end

    local prev_res
    while true 
        if succeeded 
            @debug "trying to reconstruct result via CRT"
            majority_test!(res_modp, primes_, A)

            res = crt(res_modp, primes_, A)
            clear_denominators!(res,A)
            if res == prev_res 
                return res 
            else 
                @debug "reconstructions are not consistent, trying more primes"
                prev_res = res
            end
        elseif nprime == bnd
            @debug "trying to reconstruct result via CRT"
            majority_test!(res_modp, primes_, A)
            try
                prev_res = crt(res_modp, primes_, A)
                clear_denominators!(prev_res,A)
                succeeded = true 
                @debug "success, trying one more prime"
            catch
                @debug "failed trying more primes"
                if comp(param) == :fast 
                    bnd = 2^bound + 1
                elseif comp(param) == :medium 
                    bnd = bound^2 + 1
                else # comp(param) = :slow 
                    bnd = bound + 1
                end 
                bound += 1

            end
        end
        # if nprime == 100
        #     for i in 1:5 
        #         prettyprint(res_modp[i],A)
        #     end
        #     prev_res = crt(res_modp, primes_, A)
        #     clear_denominators!(prev_res,A)
        #     succeeded = true 
        #     @debug "success, trying one more prime"
        #     error("fin")
        # end

        nprime += 1 
        globalstats.counters[:number_primes] += 1
        prime = primes[nprime]
        push!(primes_,prime)
        @debug "computing the function for $(nprime)th prime"
        if tracer(param)
            nA, mod_args = _modp_args(ctxA, prime, A, args...)
            push!(res_modp, f(trace, mod_args...))
        else
            nA, mod_args = _modp_args(ctxA, prime, A, args...)
            push!(res_modp, f(mod_args...))
        end
    end
end


function crt(vec :: Vector{OrePoly{K,M}}, primes_ :: Vector{Int}, A :: OreAlg) where {K <: RatFunModp, M}
    nvec = length(vec)
    ncoeff = length(vec[1])
    ctxA = ctx(A)
    primes = [ZZ(prime) for prime in primes_]
    prodp = prod(primes)
    cos = Vector{RatFunQQ}(undef, ncoeff)
    cs = Vector{typeof(ZZ(0))}(undef, nvec)
    for i in 1:ncoeff
        num = zero(ctxA.F)
        term1 = Nemo.numerator(coeff(vec[1], i))
        for j in 1:length(term1)
            @inbounds for k in 1:nvec
                cs[k] = ZZ(Nemo.coeff(Nemo.numerator(coeff(vec[k], i)), j).data)
            end
            c = Nemo.reconstruct(Nemo.crt(cs, primes), prodp)
            num += ctxA.F(ctxA.R([Nemo.numerator(c)], [exponent_vector(term1, j)])) / ctxA.F(Nemo.denominator(c))
        end

        den = zero(ctxA.F)
        term1 = Nemo.denominator(coeff(vec[1], i))
        for j in 1:length(term1)
            @inbounds for k in 1:nvec
                cs[k] = ZZ(Nemo.coeff(Nemo.denominator(coeff(vec[k], i)), j).data)
            end
            c = Nemo.reconstruct(Nemo.crt(cs, primes), prodp)
            den += ctxA.F(ctxA.R([Nemo.numerator(c)], [exponent_vector(term1, j)])) / ctxA.F(Nemo.denominator(c))
        end

        cos[i] = num / den
    end
    return OrePoly(cos, deepcopy(mons(vec[1])))
end

function crt(vec :: Vector{OrePoly{K,M}}, primes_ :: Vector{Int}, A :: OreAlg) where {K <: UnivRatFunModp, M}
    nvec = length(vec)
    ncoeff = length(vec[1])
    ctxA = ctx(A)
    primes = [ZZ(prime) for prime in primes_]
    prodp = prod(primes)
    cos = Vector{UnivRatFunQQ}(undef, ncoeff)
    cs = Vector{typeof(ZZ(0))}(undef, nvec)
    for i in 1:ncoeff
        num = zero(ctxA.F)
        term1 = numerator(coeff(vec[1], i), false)
        for j in 0:length(term1)-1
            @inbounds for k in 1:nvec
                cs[k] = ZZ(Nemo.coeff(numerator(coeff(vec[k], i), false), j).data)
            end
            c = Nemo.reconstruct(Nemo.crt(cs, primes), prodp)
            num += ctxA.F(numerator(c, false) * ctxA.vars[1]^j) / ctxA.F(Nemo.denominator(c, false))
        end

        den = zero(ctxA.F)
        term1 = Nemo.denominator(coeff(vec[1], i), false)
        for j in 0:length(term1)-1
            @inbounds for k in 1:nvec
                cs[k] = ZZ(Nemo.coeff(Nemo.denominator(coeff(vec[k], i), false), j).data)
            end
            c = Nemo.reconstruct(Nemo.crt(cs, primes), prodp)
            den += ctxA.F(numerator(c, false) * ctxA.vars[1]^j) / ctxA.F(Nemo.denominator(c, false))
        end

        cos[i] = num / den
    end
    return OrePoly(cos, deepcopy(mons(vec[1])))
end

function crt(vec :: Vector{OrePoly{K,M}}, primes_ :: Vector{Int},:: OreAlg{T,C,M,O}) where {K,M,T,O,C<: QQCtx}
    nvec = length(vec)
    ncoeff = length(vec[1])
    primes = [ZZ(prime) for prime in primes_]
    cos = Vector{QQFieldElem}(undef, ncoeff)
    prodp = prod(primes)
    cs = Vector{typeof(ZZ(0))}(undef, nvec)
    for i in 1:ncoeff
        @inbounds for k in 1:nvec
            cs[k] = ZZ(coeff(vec[k], i))
        end
        c = Nemo.crt(cs, primes)
        cos[i] = Nemo.reconstruct(c, prodp)
    end
    return OrePoly(cos, deepcopy(mons(vec[1])))
end

function crt(vec :: Vector{Vector{OrePoly{K,M}}}, primes :: Vector{Int}, A :: OreAlg) where {K, M}
    tmp = [crt([vec[j][i] for j in 1:length(vec)],primes,A) for i in 1:length(vec[1])]
    return tmp
end
