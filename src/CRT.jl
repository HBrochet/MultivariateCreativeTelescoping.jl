function compute_with_CRT(f :: Function, A :: OreAlg, args...)
    nprime = 1
    globalstats.counters[:number_primes] += 1
    primes_ = [primes[1]]
    prime = primes[1]
    bound = 1
    succeeded = false # after reconstructing the result we try one more prime 

    if ctx(A) isa RatFunCtx
        nA = change_alg_char_ratfun(prime,A)
        R = Native.GF(prime)
        res_modp = [f(change_coefficient_field(R,nA, args...)...)]
    elseif ctx(A) isa QQCtx
        nA = change_alg_char_QQ(prime,A)
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
        elseif nprime == bound^2
            @debug "trying to reconstruct result via CRT"
            try
                prev_res = crt(res_modp, primes_, A)
                clear_denominators!(prev_res,A)
                succeeded = true 
                @debug "success, trying one more prime"
            catch
                @debug "failed trying more primes"
                bound += 1 
            end
        end
        # if nprime == 20
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
            push!(res_modp, f(change_coefficient_field(R,nA, args...)...))
        elseif ctx(A) isa QQCtx
            nA = change_alg_char_QQ(prime,A)
            push!(res_modp, f(change_coefficient_field(nA,args...)...))
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
        term1 = numerator(coeff(vec[1],i))
        for j in 1:length(term1)
            cs = [ZZ(Nemo.coeff(numerator(coeff(vec[k],i)), j).data) for k in 1:length(vec)]
            c = Nemo.crt(cs, primes)
            c = Nemo.reconstruct(c,prodp)
            num += ctx(A).F(ctx(A).R([numerator(c)],[exponent_vector(term1,j)])) / ctx(A).F(Nemo.denominator(c))
        end

        den = zero(ctx(A).F)
        term1 = Nemo.denominator(coeff(vec[1],i))
        for j in 1:length(term1)
            cs = [ZZ(Nemo.coeff(Nemo.denominator(coeff(vec[k],i)), j).data) for k in 1:length(vec)]
            c = Nemo.crt(cs, primes)
            c = Nemo.reconstruct(c,prodp)
            den += ctx(A).F(ctx(A).R([numerator(c)],[exponent_vector(term1,j)])) / ctx(A).F(Nemo.denominator(c))
        end

        push!(cos, num/den)
    end
    return OrePoly(cos,deepcopy(mons(vec[1])))
end

function crt(vec :: Vector{OrePoly{K,M}}, primes_ :: Vector{Int},:: OreAlg{T,C,M}) where {K,M,T,C<: QQCtx}
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
