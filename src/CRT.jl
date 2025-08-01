struct CRTParam{A,B} end 

function crt_param(;tracer :: Val{A} = Val(false),
                    comp :: Val{B} = Val(:medium)
                    ) where {A,B}
    return CRTParam{A,B}()
end

tracer(p :: CRTParam{A,B}) where {A,B} = A 
comp(p :: CRTParam{A,B}) where {A,B} = B


function compute_with_CRT(f :: Function, A :: OreAlg, args...;param ::CRTParam = crt_param())
    nprime = 1
    globalstats.counters[:number_primes] += 1
    primes_ = [primes[1]]
    prime = primes[1]
    bound = 1
    bnd = 1
    succeeded = false # after reconstructing the result we try one more prime 

    if ctx(A) isa RatFunCtx
        nA = change_alg_char_ratfun(prime,A)
        R = Native.GF(prime)
        if tracer(param)
            tmp, trace = f(change_coefficient_field(R,nA, args...)...;tracer=Val(true))
            res_modp = [tmp]
        else 
            res_modp = [f(change_coefficient_field(R,nA, args...)...)]
        end
    elseif ctx(A) isa QQCtx
        nA = change_alg_char_QQ(prime,A)
        if tracer(param)
            tmp, trace = f(change_coefficient_field(nA, args...)...;tracer=Val(true))
            res_modp = [tmp]
        else 
            res_modp = [f(change_coefficient_field(nA, args...)...)]
        end
        res_modp = [f(change_coefficient_field(nA,args...)...)]
    end
    let prev_res
    while true 
        if succeeded 
            @debug "trying to reconstruct result via CRT"
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
        # if nprime == 5
        #     # for i in 1:5 
        #     #     prettyprint(res_modp[i],A)
        #     # end
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
        @debug "computing ggd for $(nprime)th prime"
        if ctx(A) isa RatFunCtx
            nA = change_alg_char_ratfun(prime,A)
            R = Native.GF(prime)
            if tracer(param)
                push!(res_modp, f(trace,change_coefficient_field(R,nA, args...)...))
            else 
                push!(res_modp, f(change_coefficient_field(R,nA, args...)...))
            end
        elseif ctx(A) isa QQCtx
            nA = change_alg_char_QQ(prime,A)
            if tracer(param)
                push!(res_modp, f(trace,change_coefficient_field(nA,args...)...))
            else 
                push!(res_modp, f(change_coefficient_field(nA,args...)...))
            end
        end
    end
    end
end


function crt(vec :: Vector{OrePoly{K,M}}, primes_ :: Vector{Int}, A :: OreAlg) where {K <: RatFunModp, M}
    primes = [ZZ(prime) for prime in primes_]
    cos = RatFunQQ[]
    prodp = prod(primes)
    for i in 1:length(vec[1])
        num = zero(ctx(A).F)
        term1 = Nemo.numerator(coeff(vec[1],i))
        for j in 1:length(term1)
            cs = [ZZ(Nemo.coeff(Nemo.numerator(coeff(vec[k],i)), j).data) for k in 1:length(vec)]
            c = Nemo.crt(cs, primes)
            c = Nemo.reconstruct(c,prodp)
            num += ctx(A).F(ctx(A).R([Nemo.numerator(c)],[exponent_vector(term1,j)])) / ctx(A).F(Nemo.denominator(c))
        end

        den = zero(ctx(A).F)
        term1 = Nemo.denominator(coeff(vec[1],i))
        for j in 1:length(term1)
            cs = [ZZ(Nemo.coeff(Nemo.denominator(coeff(vec[k],i)), j).data) for k in 1:length(vec)]
            c = Nemo.crt(cs, primes)
            c = Nemo.reconstruct(c,prodp)
            den += ctx(A).F(ctx(A).R([Nemo.numerator(c)],[exponent_vector(term1,j)])) / ctx(A).F(Nemo.denominator(c))
        end

        push!(cos, num/den)
    end
    return OrePoly(cos,deepcopy(mons(vec[1])))
end

function crt(vec :: Vector{OrePoly{K,M}}, primes_ :: Vector{Int}, A :: OreAlg) where {K <: UnivRatFunModp, M}
    primes = [ZZ(prime) for prime in primes_]
    cos = UnivRatFunQQ[]
    prodp = prod(primes)
    for i in 1:length(vec[1])
        num = zero(ctx(A).F)
        term1 = numerator(coeff(vec[1],i),false)
        for j in 0:length(term1)-1
            cs = [ZZ(Nemo.coeff(numerator(coeff(vec[k],i),false), j).data) for k in 1:length(vec)]
            c = Nemo.crt(cs, primes)
            c = Nemo.reconstruct(c,prodp)
            num += ctx(A).F(numerator(c,false)*ctx(A).vars[1]^j) / ctx(A).F(Nemo.denominator(c,false))
        end

        den = zero(ctx(A).F)
        term1 = Nemo.denominator(coeff(vec[1],i),false)
        for j in 0:length(term1)-1
            cs = [ZZ(Nemo.coeff(Nemo.denominator(coeff(vec[k],i),false), j).data) for k in 1:length(vec)]
            c = Nemo.crt(cs, primes)
            c = Nemo.reconstruct(c,prodp)
            den += ctx(A).F(numerator(c,false)*ctx(A).vars[1]^j) / ctx(A).F(Nemo.denominator(c,false))
        end

        push!(cos, num/den)
    end
    return OrePoly(cos,deepcopy(mons(vec[1])))
end

function crt(vec :: Vector{OrePoly{K,M}}, primes_ :: Vector{Int},:: OreAlg{T,C,M,O}) where {K,M,T,O,C<: QQCtx}
    primes = [ZZ(prime) for prime in primes_]
    cos = Vector{QQFieldElem}(undef,length(vec[1]))
    prodp = prod(primes)
    for i in 1:length(vec[1])
        cs = [ZZ(coeff(vec[k],i)) for k in 1:length(vec)]
        c = Nemo.crt(cs, primes)
        cos[i] = Nemo.reconstruct(c,prodp)
    end
    return OrePoly(cos,deepcopy(mons(vec[1])))
end

function crt(vec :: Vector{Vector{OrePoly{K,M}}}, primes :: Vector{Int}, A :: OreAlg) where {K, M}
    tmp = [crt([vec[j][i] for j in 1:length(vec)],primes,A) for i in 1:length(vec[1])]
    return tmp
end
