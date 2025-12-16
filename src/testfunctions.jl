function ann_inv_pol(p :: OrePoly{K,M},A :: OreAlg) where {K,M}
    g = OrePoly{K,M}[]
    for i in 1:A.nrdv+A.npdv 
        d = makepoly(one(ctx(A)),makemon(i,A))
        push!(g, mul(d,p,A))
    end
    indT = A.nrdv + A.npdv*2 + A.npv + 1
    T = makepoly(one(ctx(A)),makemon(indT,A))
    rel = mul(p,T,A)
    rel = sub!(rel,one(A),A)
    push!(g,rel)
    return saturation(g,indT,A,wc_trace())
end


function OrePoly_to_flint_poly(pol :: OrePoly{T,M},A) where {T <: Integer, M} 
    str = [s for s in A.indexp_to_strvar[A.npdv+1:A.npdv*2]]
    R,vars = polynomial_ring(Nemo.QQ,str)
    F = fraction_field(R)
    p = Int(ctx(A).char) 
    res = F(0) 
    for (co,mo) in pol
        tmp = F(Nemo.reconstruct(Int(co),p))
        for k in 1:A.npdv
            tmp = tmp * F(vars[k])^mo[A.npdv + k]
        end
        res += tmp 
    end
    return res, F ,vars
end


function check_ann(g :: Vector{OrePoly{T,K}},p_ :: OrePoly, A :: OreAlg) where {T,K}
    p, F, vars = OrePoly_to_flint_poly(p_,A)
    prime = Int(ctx(A).char) 
    for (l,pol) in enumerate(g) 
        tmp = F(0)
        for (co,mo) in pol 
            tmp2 = F(Nemo.reconstruct(Int(co),prime)) / p
            for j in A.npdv+1:2*A.npdv 
                tmp2 = tmp2 * F(vars[j-A.npdv])^mo[j]
            end
            for j in 1:A.npdv 
                for k in 1:mo[j]
                    tmp2 = derivative(tmp2,j)
                end           
            end
            tmp += tmp2
        end
        if tmp != 0 
            println("OrePoly $l does not annihilate 1/p")
            return false
        end
    end
    return true 
end
