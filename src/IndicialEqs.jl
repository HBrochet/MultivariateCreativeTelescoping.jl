

function weight(m :: OreMonVE{N,E}, A :: OreAlg) where {N,E}
    let w, ind 
        ord = order(A)
        if ord isa Weight
            w = order(A).weights
            ind = SVector{N,E}(i for i in 1:N)
        else 
            while ord isa Block 
                if ord.o1 isa Weight
                    w = ord.o1.weights
                    ind = ord.ind1
                    break 
                end
                ord = ord.o2 
            end
        end
        return sum(w[i]*m[ind[i]] for i in 1:length(w))
    end
end

function initialform(p :: OrePoly, A::OreAlg) 
    max = weight(mon(p,1),A)
    for i in 2:length(p)
        if max > weight(mon(p,i),A)
            mons = [mon(p,j) for j in 1:i-1]
            coeffs = [coeff(p,j) for j in 1:i-1]
            return OrePoly(coeffs,mons)
        end
    end
    return p 
end

function initialforms(v :: Vector{OrePoly{K,M}}, A :: OreAlg) where {K,M} 
    return [initialform(p,A) for p in v]
end

function initialformweights(v :: Vector{OrePoly{K,M}}, A :: OreAlg) where {K,M} 
    return [weight(mon(p,1),A) for p in v]
end



# it assumes that d's are on the right 
function indicial_eq(t :: Tuple{K,OreMonVE{N,E}},A :: OreAlg) where {K,N,E}
    c,m = t 
    m2 = OreMonVE(SVector{N,E}(i <= A.nrdv + A.npdv || i > A.nrdv + 2*A.npdv ? 0 : m[i] - m[i-A.npdv] for i in 1:N))
    for i in 1:A.npdv
        expd = m[A.nrdv + i] 
        expx = m[A.nrdv + A.npdv + i]
        a = gen(ctx(A).R,i)
        for i in expx-expd+1:expx
            c = mul(c,convertn(a+i,ctx(A)),ctx(A))
        end
    end
    if isodd(sum(m[i] for i in A.nrdv+1:A.nrdv+A.npdv))
        return opp(c, ctx(A)),m2
    else
        return c,m2
    end
end


function indicial_eq(p :: OrePoly,A :: OreAlg)
    res = undefOrePoly(length(p),A)
    for i in 1:length(p) 
        res[i] =  indicial_eq(p[i],A)
    end
    normalize!(res,A)
    return res
end

function indicial_eqs(g :: Vector{OrePoly{K,M}}, A:: OreAlg) where {K,M}
    res = Vector{OrePoly{K,M}}(undef,length(g))
    for i in 1:length(g) 
        res[i] = indicial_eq(g[i],A)
    end
    return res
end


# It assumes that m is a monomial in x only and g is a vector of indicial equations returned by indicial_eqs
# function is_reducible(g :: Vector{OrePoly{K,M}},m::OreMonVE, A :: OreAlg) where {K,M}
#     for f in g 
#         for i in 1:length(f)
#             if divides(mon(f,i),m) 
#                 n = m/mon(f,i) 
                
#                     return true 
#                 end