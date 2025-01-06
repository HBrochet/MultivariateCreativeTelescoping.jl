# normalization


function normalize!(P :: OrePoly, A :: OreAlg)
    sort!(P, rev=true, order = order(A))

    # collapse like terms, remove zeros
    w = 0                   # read and write indices
    r = 1
    lenp = length(P)
    while r <= lenp
        m = mon(P, r)
        c = coeff(P, r)

        while true
            r += 1
            if r > lenp || mon(P, r) != m
                break
            end
            c = add(c,coeff(P, r),ctx(A))
        end

        if !iszero(c,A)
            w += 1
            P[w] = (c,m)
        end
    end

    resize!(coeffs(P), w)
    resize!(mons(P), w)
    return P
end




### multiplication for OrePolynomials with context Context

function add2!(P :: OrePoly, Q :: OrePoly, A :: OreAlg;normal :: Bool =true)
    append!(coeffs(P),coeffs(Q))
    append!(mons(P),mons(Q))
    if normal
        normalize!(P,A)
    end
    return P
end

function add!(P :: OrePoly, Q :: OrePoly, A :: OreAlg; normal :: Bool = true)
    lenp = length(P)
    lenq = length(Q)
    len = lenp + lenq
    res = undefOrePoly(len,A)
    ip = 1
    iq = 1
    s = 1
    while ip <= lenp && iq <= lenq
        if lt(order(A),mon(Q,iq),mon(P,ip))
            res[s] = P[ip]
            ip += 1
            s += 1
        elseif mon(P,ip) == mon(Q,iq)
            c = add(coeff(P,ip), coeff(Q,iq), ctx(A))
            if !iszero(c,A)
                res[s] = (c, mon(P,ip))
                s +=1
            end
            ip +=1
            iq +=1
        else 
            res[s] = Q[iq]
            iq += 1
            s += 1
        end
    end
    while iq <= lenq
        res[s] = Q[iq]
        iq +=1
        s +=1
    end
    while ip <= lenp 
        res[s] = P[ip]
        ip +=1
        s +=1
    end

    resize!(coeffs(res), s-1)
    resize!(mons(res), s-1)
    # @assert all(!(c==zero(ctx(A))) for c in coeffs(P))
    # @assert all(!(c==zero(ctx(A))) for c in coeffs(Q))
    # @assert all(!(c==zero(ctx(A))) for c in coeffs(res))
    return res
end


# function add_in_place!(P ::OrePoly, Q :: OrePoly, A::OreAlg)
#     ctx = ctx(A)
#     lenp = length(P)
#     lenq = length(Q)

#     ipr = 1 # reading index in P
#     ipw = 1 # writting index in P
#     iq = 1 # reading index in Q 
#     # R is a temporary polynomial that we store at the beginning of Q 
#     irr = 1 # reading index of R in Q
#     irw = 1 # writing index of R in Q

    # while ipr <= lenp && iqr <= lenq 
    #     if lt(order(A), mon(Q,iq),mon(P,ip))
    #         P[ipw] = P[ipr]
    #         ipw += 1
    #         ipr += 1
    #     elseif mon(P,ip) == mon(Q,iq) 
    #         c = add(P.coeffs[ipr],Q.coeffs[iq],ctx)
    #         if !iszero(c)
    #             P[ipw] = (c,mp) 
    #             ipw += 1
    #         end
    #         ipr += 1
    #         iq += 1
    #     elseif ipw < ipr  # we know that mq > mp
    #         P[ipw]= Q[iq]
    #         iq += 1
    #         ipw += 1 
    #     else
    #         tmp = Q[iq]
    #         Q[irw] = P[ipr]
    #         P[ipw] = tmp

    #         irw += 1
    #         iq += 1
    #         ipw += 1
    #         ipr += 1

    #         # by construction we allways have mr > mp
    #         while irw - irr > 0 && ipr <= lenp && iqr <= lenq
    #             mr = mon(Q,irr)
    #             mq = mon(Q,iq)
    #             if lt(order(A), mr, mq)
    #                 tmp = (coeff(Q,iq), mq)
    #                 Q[irw] = P[ipr]
    #                 P[ipw] = tmp
        
    #                 iq += 1
    #                 ipw += 1
    #                 ipr += 1
    #                 irw += 1

    #             elseif mr == mq 
    #                 c = coeff(Q,iq) + coeff(Q,irr)
    #                 if !iszero(c)
    #                     tmp = P[ipr]
    #                     P[ipw] = (c,mr) 
    #                     Q[ipw] = tmp
    #                     ipw += 1
    #                     ipr += 1
    #                     irw +=1
    #                 end
    #                 iq += 1 
    #                 irr += 1
    #             else # mq < mr 
    #                 c = coeff(Q,irr)
                    

    #     end
    # end


"""
    add(p :: OrePoly, q :: OrePoly, A :: OreAlg)

Return the Ore polynomial p + q.
"""
function add(PP :: OrePoly, Q :: OrePoly, A :: OreAlg)
    P = copy(PP)
    P = add!(P,Q,A)
    return P
end

function add2(PP :: OrePoly, Q :: OrePoly, A :: OreAlg)
    P = copy(PP)
    P = add2!(P,Q,A)
    return P
end

function add(v_ :: Vector{OrePoly{T,M}},A :: OreAlg) where {T,M}
    v = v_
    len = length(v)

    if len == 0 
        return zero(A)
    end
    while len > 1
        n = div(len,2)
        tmp = OrePoly{T,M}[]
        for i in 1:n
            push!(tmp,add!(v[2*i-1],v[2*i],A))
        end
        if isodd(len)
            push!(tmp,v[end])
        end
        v = tmp
        len = length(v)
    end
    return v[1]
end


function mul(c :: T, P :: OrePoly, A :: OreAlg{T, C, M,O}) where {T, C, M,O}
    res = zero(A)
    append!(mons(res),mons(P))
    resize!(coeffs(res),length(coeffs(P)))
    cs = coeffs(res)
    for i in eachindex(coeffs(P))
        cs[i] = mul(c, coeff(P,i), ctx(A))
    end
    return res
end

function mul!(c :: T, P :: OrePoly, A :: OreAlg{T, C, M,O}) where {T, C, M,O}
    cs = coeffs(P)
    for i in eachindex(coeffs(P))
        cs[i] = mul(c, coeff(P,i), ctx(A))
    end
    return P
end

function makemonic!(P :: OrePoly, A :: OreAlg)
    mul!(inv(coeff(P,1),ctx(A)),P,A)
    return P
end

function sub!(P :: OrePoly, Q :: OrePoly, A :: OreAlg; normal :: Bool = true)
    lenp = length(P)
    lenq = length(Q)
    len = lenp + lenq
    res = undefOrePoly(len,A)
    ip = 1
    iq = 1
    s = 1

    while ip <= lenp && iq <= lenq
        if lt(order(A),mon(Q,iq),mon(P,ip))
            res[s] = P[ip]
            ip += 1
            s += 1
        elseif mon(P,ip) == mon(Q,iq)
            c = sub(coeff(P,ip), coeff(Q,iq), ctx(A))
            if !iszero(c,A)
                res[s] = (c, mon(P,ip))
                s +=1
            end
            ip +=1
            iq +=1
        else 
            res[s] = (opp(coeff(Q,iq),ctx(A)), mon(Q,iq))
            iq += 1
            s += 1
        end
    end
    while iq <= lenq
        res[s] = (opp(coeff(Q,iq),ctx(A)), mon(Q,iq))
        iq +=1
        s +=1
    end
    while ip <= lenp 
        res[s] = P[ip]
        ip +=1
        s +=1
    end

    resize!(coeffs(res), s-1)
    resize!(mons(res), s-1)
    return res
end

function sub2!(P :: OrePoly, Q :: OrePoly, A :: OreAlg)
    append!(coeffs(P), [opp(coeff(Q,i),ctx(A)) for i in 1:length(Q)])
    append!(mons(P),mons(Q))
    return normalize!(P,A)
end

"""
    sub(p :: OrePoly, q :: OrePoly, A :: OreAlg)

Return the Ore polynomial p - q.
"""
function sub(PP :: OrePoly, Q :: OrePoly, A :: OreAlg)
    P = copy(PP)
    return sub!(P,Q,A)
end


"""
    mul(p :: OrePoly, q :: OrePoly, A :: OreAlg)

Return the Ore polynomial pq.
"""
function mul(P :: OrePoly, Q :: OrePoly, A :: OreAlg; normal = true)
    vec = typeof(P)[]
    for i in 1:length(P)
        for j in 1:length(Q)
            push!(vec,mul(P[i],Q[j],A))
        end
    end
    return add(vec,A)
end

function mul(T ::Tuple{K,M}, Q :: OrePoly, A :: OreAlg; normal = true) where {K,M}
    vec = typeof(Q)[]
    for j in 1:length(Q)
        push!(vec,mul(T,Q[j],A))
    end
    return add(vec,A)
end

function mul(P :: OrePoly, T ::Tuple{K,M}, A :: OreAlg; normal = true) where {K,M}
    vec = typeof(P)[]
    for i in 1:length(P)
        push!(vec,mul(P[i],T,A))
    end
    return add(vec,A)
end

function mul(m :: OreMonVE, p :: OrePoly, A :: OreAlg; normal = true)
    return mul((one(ctx(A)),m),p,A,normal = normal)
end


function mul( T1 ::Tuple{K,M}, T2 :: Tuple{K2,M}, A :: OreAlg; normal=true) where {M,K,K2}
    res = mul_(T1,T2,A)
    if normal
        normalize!(res,A)
    end
    return res
end

function mul_( T1 ::Tuple{K,M}, T2 :: Tuple{K2,M}, A :: OreAlg) where {M,K,K2}
    #noncomm = Int[]
    s = 1 # size of the product before the call to normalize!

    for l in 1:A.nlv # for each loc vars
        if T1[2][l + A.nrdv+2*A.npdv] > 0
            for k in 1:A.nrdv + A.npdv
                if T2[2][k] > 0
                    indT = l + A.nrdv+2*A.npdv+A.npv
                    T = makemon(indT,A)
                    dt = makemon(k, A)

                    Tpower = T^T1[2][indT]
                    fact1 = (T1[1], T1[2]/ Tpower)
                    fact2 = makepoly(T2[1],dt*Tpower)

                    tmp = mul(A.diff_pols_loc[l][k], makepoly(T2[1]*convertn(Int(T1[2][indT]),ctx(A)),T^(T1[2][indT]+1)),A,normal=true)
                    fact2 = add!(fact2, tmp, A,normal=true)


                    fact3 = (one(ctx(A)), T2[2]/dt)
                    tmp = mul(fact1,mul(fact2,fact3,A,normal=true),A,normal=true)
                    return tmp
                end
            end
        end
    end

    if ctx(A) isa RatFunCtx
        for l in 1:A.nrdv
            if T1[2][l] == 0; continue end
            if ctx(A) isa MRatFunCtx
                c = derivative(T2[1],ctx(A).vars[l])
            else 
                c = derivative(T2[1])
            end
            if iszero(c,A); continue end
            res = undefOrePoly(T1[2][l]+1,A)
            dt = makemon(l,A)
            res[1] = (T2[1], T2[2] * dt^T1[2][l])
            ctr = 2
            for i in T1[2][l]-1:-1:0
                res[ctr] = (c*binomial(Int(T1[2][l]),i), T2[2] * dt^i)
                if ctx(A) isa MRatFunCtx
                    c = derivative(c,ctx(A).vars[l])
                else 
                    c = derivative(c)
                end
                ctr += 1
            end
            return mul((T1[1],T1[2] / dt^T1[2][l]), res, A,normal=true)
        end 
    end

    if A.varord == "dleft"
        for l in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv
            tmp = min(T1[2][l],T2[2][l-A.npdv])
            s = s * (tmp + 1)
            if s < 0
                println(T1[2][l], " ", T2[2][l-A.npdv])
            end
        end

        res = undefOrePoly(s,A)
        res[1] = (mul(T1[1],T2[1],ctx(A)), T1[2] * T2[2])

        w = 2 # write index in res
        r = 1

    
        for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv
            if T1[2][i] == 0 || T2[2][i-A.npdv] == 0; continue end

            fact = convertn(1,ctx(A)) # to store decreasing factorials
            for j in 1:min(T1[2][i],T2[2][i-A.npdv])
                fact = mul(fact,convertn(T1[2][i]-j+1,ctx(A)),ctx(A))
                fact = opp(fact,ctx(A))
                for l in 1:r
                    c = mul(mul(fact,
                                convertn(binomial(Int(T2[2][i-A.npdv]), j), ctx(A)),
                                ctx(A)),
                            res[l][1], ctx(A))

                    m = res[l][2] / (makemon(i,A) * makemon(i-A.npdv,A))^j
                    res[w] = (c, m)

                    w = w + 1
                end
            end
            r = w-1
        end
    else 
        for l in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv
            tmp = min(T1[2][l-A.npdv],T2[2][l])
            s = s * (tmp + 1)
            if s < 0
                println(T1[2][l], " ", T2[2][l-A.npdv])
            end
        end

        res = undefOrePoly(s,A)
        res[1] = (mul(T1[1],T2[1],ctx(A)), T1[2] * T2[2])

        w = 2 # write index in res
        r = 1

        for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv
            if T1[2][i-A.npdv] == 0 || T2[2][i] == 0; continue end

            fact = convertn(1,ctx(A)) # to store decreasing factorials
            for j in 1:min(T1[2][i-A.npdv],T2[2][i])
                fact = mul(fact,convertn(T1[2][i-A.npdv]-j+1,ctx(A)),ctx(A))
                for l in 1:r
                    c = mul(mul(fact,
                                convertn(binomial(Int(T2[2][i]), j), ctx(A)),
                                ctx(A)),
                            res[l][1], ctx(A))

                    m = res[l][2] / (makemon(i-A.npdv,A) * makemon(i,A))^j
                    res[w] = (c, m)

                    w = w + 1
                end
            end
            r = w-1
        end
    end
    return res
end

function shift(P :: OrePoly{K,M}, m :: M, A :: alg) where {K,M,alg <: OreAlg}
    return mul(m/mon(P,1),P,A)
end
