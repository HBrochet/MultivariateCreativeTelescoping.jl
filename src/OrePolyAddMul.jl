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

        if !iszero(c)
            w += 1
            P[w] = (c,m)
        end
    end

    resize!(coeffs(P), w)
    resize!(mons(P), w)
    return P
end




### multiplication for OrePolynomials with context Context

function add!(P :: OrePoly, Q :: OrePoly, A :: OreAlg;normal :: Bool =true)
    append!(coeffs(P),coeffs(Q))
    append!(mons(P),mons(Q))
    if normal
        normalize!(P,A)
    end
    return P
end

function add2!(P :: OrePoly, Q :: OrePoly, A :: OreAlg)
    lenp = length(P)
    lenq = length(Q)
    len = lenp + lenq
    res = undefOrePoly(len,A)
    ip = 1
    iq = 1
    s = 1
    while ip <= lenp && iq <= lenq
        if lt(order(A),mon(P,ip),mon(Q,iq))
            c = coeff(Q,iq)
            m = mon(Q,iq)
            iq +=1
        else
            c = coeff(P,ip)
            m = mon(P,ip)
            ip +=1
        end
        while iq <= lenq && m == mon(Q,iq)
            c = add(c, coeff(Q, iq),ctx(A))
            iq +=1
        end
        while ip <= lenp && m == mon(P,ip)
            c = add(c, coeff(P, ip),ctx(A))
            ip +=1
        end
        # println("partial res ")
        # println(c,m)
        if !iszero(c,ctx(A))
            res[s] = (c,m)
            s += 1
        end
    end
    while ip <= lenp
        res[s] = (coeff(P,ip),mon(P,ip))
        s += 1
        ip += 1
    end
    while iq <= lenq
        res[s] = (coeff(Q,iq),mon(Q,iq))
        s += 1
        iq += 1
    end
    resize!(coeffs(res), s-1)
    resize!(mons(res), s-1)
    return res
end

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


function mul(c :: T, P :: OrePoly, A :: OreAlg{T, C, M}) where {T, C, M}
    res = zero(A)
    append!(mons(res),mons(P))
    resize!(coeffs(res),length(coeffs(P)))
    cs = coeffs(res)
    for i in eachindex(coeffs(P))
        cs[i] = mul(c, coeff(P,i), ctx(A))
    end
    return res
end

function mul!(c :: T, P :: OrePoly, A :: OreAlg{T, C, M}) where {T, C, M}
    cs = coeffs(P)
    for i in eachindex(coeffs(P))
        cs[i] = mul(c, coeff(P,i), ctx(A))
    end
end

function makemonic!(P :: OrePoly, A :: OreAlg)
    mul!(inv(coeff(P,1),ctx(A)),P,A)
end

function sub!(P :: OrePoly, Q :: OrePoly, A :: OreAlg)
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
    res = zero(A)
    for i in 1:length(P)
        for j in 1:length(Q)
            tmp = mul(P[i],Q[j],A,normal=false)
            res = add!(res,tmp,A,normal = false)
        end
    end
    if normal
        normalize!(res,A)
    end
    return res
end

function mul(T ::Tuple{K,M}, Q :: OrePoly, A :: OreAlg; normal = true) where {K,M}
    res = zero(A)
    for i in 1:length(Q)
        res = add!(res,mul(T,Q[i],A,normal = false),A,normal=false)
    end
    if normal
        normalize!(res,A)
    end
    return res
end

function mul(P :: OrePoly, T ::Tuple{K,M}, A :: OreAlg; normal = true) where {K,M}
    res = zero(A)
    for i in 1:length(P)
        res = add!(res,mul(P[i],T,A,normal = false),A,normal=false)
    end
    if normal
        normalize!(res,A)
    end
    return res
end

function mul(m :: OreMonVE, p :: OrePoly, A :: OreAlg; normal = true)
    return mul((one(ctx(A)),m),p,A,normal = normal)
end


function mul( T1 ::Tuple{K,M}, T2 :: Tuple{K,M}, A :: OreAlg; normal=true) where {M,K}
    res = mul_(T1,T2,A)
    if normal
        normalize!(res,A)
    end
    return res
end

function mul_( T1 ::Tuple{K,M}, T2 :: Tuple{K,M}, A :: OreAlg) where {M,K}
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

                    tmp = mul(A.diff_pols_loc[l][k], makepoly(T2[1]*convertn(T1[2][indT],ctx(A)),T^(T1[2][indT]+1)),A,normal=false)
                    fact2 = add!(fact2, tmp, A,normal=false)


                    fact3 = (one(ctx(A)), T2[2]/dt)
                    tmp = mul(fact1,mul(fact2,fact3,A,normal=false),A,normal=false)
                    return tmp
                end
            end
        end
    end

    if ctx(A) isa RatFunCtx
        for l in 1:A.nrdv
            if T1[2][l] == 0; continue end
            c = derivative(T2[1],ctx(A).vars[l])
            if iszero(c,ctx(A)); continue end
            res = undefOrePoly(T1[2][l]+1,A)
            dt = makemon(l,A)
            res[1] = (T2[1], T2[2] * dt^T1[2][l])
            ctr = 2
            for i in T1[2][l]-1:-1:0
                res[ctr] = (c*binomial(Int(T1[2][l]),i), T2[2] * dt^i)
                c = derivative(c,ctx(A).vars[l])
                ctr += 1
            end
            return mul((T1[1],T1[2] / dt^T1[2][l]), res, A,normal=false)
        end 
    end

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
    return res
end

function shift(P :: OrePoly{K,M}, m :: M, A :: alg) where {K,M,alg <: OreAlg}
    return mul(m/mon(P,1),P,A)
end
