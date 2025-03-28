abstract type GBParam end 

struct DefaultParam <: GBParam end 
use_geobucket(:: DefaultParam) = true
stophol(:: DefaultParam)= false
stat(:: DefaultParam) = false
debug(::DefaultParam) = false


# Reduce mon r of pol f with pol g 
function reduce!(f :: OrePoly, r :: Int, g :: OrePoly, A :: OreAlg;param :: GBParam = DefaultParam()) 
    mo = mon(f,r) / mon(g,1)
    co = mul(opp(coeff(f,r),ctx(A)), inv(coeff(g,1),ctx(A)),ctx(A))
    np = makepoly(co,mo)
    res = mul(np, g, A)
    tmp = add!(f, res, A) 

    if stat(param)
        globalstats.counters[:bbg_size_reducer] += length(res)
        globalstats.counters[:bbg_size_m]+= sum(mo)
        globalstats.counters[:bbg_deg_reducer] += maxdeg(res)
    end
    return tmp
end



function div!(f :: OrePoly, g :: Vector{OrePoly{K,M}} ,A :: OreAlg; full :: Val{B} = Val(true),param :: GBParam = DefaultParam()) where {K,M,B}
    r=1 
    while r <= length(f) 
        div = false
        for i in 1:length(g)
            if divide(mon(g[i],1), mon(f,r),A) && iscompatible(mon(g[i],1), mon(f,r),A)
                div = true
                if stat(param)
                    globalstats.counters[:bbg_nb_div] += 1
                    add_reducer_globalstats!(i, mon(f,r)/mon(g[i],1))
                end
                
                f = reduce!(f,r, g[i], A,param = param)

                break
            end
        end
        if !div 
            if !B 
                return f
            end
            r = r + 1
        end
    end
    return f
end




function spair(f :: OrePoly, g :: OrePoly, A :: OreAlg)
    l = lcm(mon(f,1), mon(g,1))
    gg = makepoly( opp(coeff(f,1),ctx(A)), l / mon(g,1))
    gg = mul(gg, g, A)
    ff = makepoly(coeff(g,1),l / mon(f,1))
    ff = mul(ff, f, A)

    ff = add!(ff, gg, A)
    return ff
end





function buchberger(f :: Vector{OrePoly{K,M}}, A :: OreAlg) where {K,M}
    g = deepcopy(f)
    leng = length(g)
    pairs = [(i,j) for i in 1:length(g) for j in i+1:length(g)]
    npairs = div(leng * (leng - 1), 2)
    while npairs > 0 
        pair = pop!(pairs)
        npairs = npairs - 1
        r = spair(g[pair[1]],g[pair[2]],A)
        r = div!(r,g,A)
        if length(r) > 0 
            leng = leng + 1 
            makemonic!(r,A)
            push!(g, r)
            append!(pairs,[(i,leng) for i in 1:leng-1])
            npairs = npairs + leng - 1
            # ntotpairs = ntotpairs + leng -1
        else 
            # nredzero = nredzero + 1
        end
    end
    return g
end

# It assumes that every pivots have been found
function reducebasis!(f :: Vector{OrePoly{K,M}}, A :: Alg;param :: GBParam = DefaultParam()) where {K,M, Alg <:OreAlg}
    i = 1

    lenf = length(f)
    
    while i <= lenf
        makemonic!(f[i],A)
        tmp = f[lenf]
        f[lenf] = f[i]
        f[i] = tmp
        f[lenf] = div!(f[lenf],f[1:lenf-1],A,param=param)
        if length(f[lenf]) == 0 
            lenf = lenf -1 
            pop!(f)
        else
            tmp = f[lenf]
            f[lenf] = f[i]
            f[i] = tmp
            i = i + 1
        end
    end

    sort!(f, lt = (x,y) -> lt(order(A),x[1][2], y[1][2]), rev = true)
end

function divbyGB(v :: Vector{OrePoly{T,M}}, gb :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T, M}
    vv = deepcopy(v)
    for i in 1:length(vv)
        div!(vv[i],gb,A)
    end
    return vv
end