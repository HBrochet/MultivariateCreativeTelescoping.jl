struct F5Param{B,C,D,E} <: GBParam end 

f5_param(;geobucket ::Val{B} = Val(true),
         stophol  :: Val{C} = Val(false),
         stat :: Val{D} = Val(false),
         debug :: Val{E} = Val(false)) where {B,C,D,E} = F5Param{B,C,D,E}() 

use_geobucket(:: F5Param{B,C,D,E}) where {B,C,D,E} = B
stophol(:: F5Param{B,C,D,E}) where {B,C,D,E}= C
stat(:: F5Param{B,C,D,E}) where {B,C,D,E} = D 
debug(::F5Param{B,C,D,E}) where {B,C,D,E} = E



function prebasis(gen :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T, M}
    N = nvars(M)
    sgb = SigPair{N, T, M}[]
    for i in 1:length(gen)
        lmon = lm(gen[i])
        sig = Signature(lmon,i)
        push!(sgb, SigPair(gen[i], sig))
    end
    return sgb
end



function updateQ!(Q :: SortedSet, sgb :: Vector{SigPair{N,T,M}}, 
                    g :: SigPair{N,T,M}, sorder :: SigOrder, A :: OreAlg) where {N, T, M}
    for i in 1:length(sgb)
        if !iscompatible(mon(sgb[i],1),mon(g,1),A)
            continue
        end
        thelcm = lcm(mon(sgb[i],1), mon(g,1))
        si = mul(thelcm/mon(sgb[i],1),sgb[i].sig)
        sg = mul(thelcm/mon(g,1),g.sig)
        if lt(sorder,si, sg)
            push!(Q, sg)
        elseif lt(sorder,sg, si)
            push!(Q, si)
        end
    end
    todelete = Set{Signature{N}}()
    for sig in Q
        for tau in Q    
            if sig != tau && divide(tau,sig,A)
                push!(todelete,sig)
            end
        end
    end
    
    globalstats.counters[:f5_eliminated_signatures] += length(todelete)

    for s in todelete
        delete!(Q,s)
    end
end


function selectreductant(sgb :: Vector{SigPair{N,T,M}},sigma :: Signature{N},A :: OreAlg) where {T,M,N}
    for i in length(sgb):-1:1
        if iscompatible(sgb[i].sig,sigma,A) && divide(sgb[i].sig,sigma,A)
            mon = sigma.mon/sgb[i].sig.mon
            return mul(mon, sgb[i], A)
        end
    end
    return zeroSigPair(A)
end

function selectreductant(geob :: GeoBucket, sgb :: Vector{SigPair{N,T,M}},sigma :: Signature{N},A :: OreAlg) where {T,M,N}
    for i in length(sgb):-1:1
        if iscompatible(sgb[i].sig,sigma,A) && divide(sgb[i].sig,sigma,A)
            mon = sigma.mon/sgb[i].sig.mon
            addmul_geobucket!(geob, one(ctx(A)), mon, sgb[i].op, A)
            res = SigPair(normalform(geob,A),sigma) 
            return res
        end
    end
    return zeroSigPair(A)
end

# version without GeoBucket
function sigdiv!(redf :: SigPair{N,K,M}, sgb :: Vector{SigPair{N,K,M}}, so :: SigOrder, A :: OreAlg,param :: F5Param) where {N,K,M}
    res = redf
    r = 1
    while r <= length(res.op)
        div = false
        for i in length(sgb):-1:1
            if iscompatible(mon(sgb[i],1), mon(res,r),A) && divide(mon(sgb[i],1), mon(res,r),A) && lt(so,mul(mon(res,r)/mon(sgb[i],1), sgb[i].sig), res.sig)
                div = true
                if stat(param)
                    add_reducer_globalstats!(i, mon(res,r)/mon(sgb[i],1))
                    globalstats.counters[:f5_number_divisions]+=1
                end
                res = reduce!(res,r, sgb[i], A,param)
                break
            end
        end
        if !div 
            r = r + 1
        end

    end
    return res  
end

function add_reducer_globalstats!(i :: Int, m ::OreMonVE)
    s = Symbol(i,string(m.exp))
    if haskey(globalstats.reducers,s) 
        globalstats.reducers[s] += 1 
    else
        globalstats.reducers[s] = 1 
    end 
end

# version without GeoBucket
# Reduce mon r of pol f with pol g 
function reduce!(f :: SigPair, r :: Int, g :: SigPair, A :: OreAlg,param :: F5Param) 
    themon = mon(f,r)/mon(g,1)
    thectx = ctx(A)
    thecoeff = opp(mul(coeff(f,r),inv(coeff(g,1),thectx),thectx),thectx)
    np = OrePoly([thecoeff],[themon])
    np=mul(np, g.op, A)
    if stat(param)
        globalstats.counters[:f5_size_reducer] += length(np)
        globalstats.counters[:f5_size_m]+= sum(themon)
        globalstats.counters[:f5_deg_reducer] += maxdeg(np)
    end
    return SigPair(add!(f.op, np, A),f.sig)
end

# version with GeoBucket
function sigdiv!(geob :: GeoBucket, tmp_poly :: ReuseOrePoly, redf :: SigPair{N,K,M}, sgb :: Vector{SigPair{N,K,M}}, so :: SigOrder, A :: OreAlg, param :: F5Param) where {N,K,M,}
    sig = redf.sig

    geob = init_GeoBucket!(geob, redf.op)
    while true 
        div = false
        lco,lmon = lt(geob,A)
        if iszero(lco)
            break
        end
        for i in length(sgb):-1:1
            lmon2 = lm(sgb[i])
            if iscompatible(lmon2, lmon,A) && divide(lmon2, lmon,A) && lt(so,mul(lmon/lmon2, sgb[i].sig), sig)
                div = true
                geob = reduce_geob!(geob, (lco,lmon), sgb[i], A, Val(false))
                if stat(param)
                    globalstats.counters[:f5_number_divisions]+=1
                end
                break
            end
        end
        if !div
            push!(tmp_poly,lco,lmon)
            rem_lt!(geob,lmon)
        end
    end
    # res = normalform(geob,A)
    res = copy_to_OrePoly!(tmp_poly,A)
    return SigPair(res,sig)
end

# version SigPairGB
function reduce_geob!(g :: GeoBucket, t:: Tuple{K,M}, f :: SigPair, A :: OreAlg, :: Val{B})  where {K, M, B}
    (c,m) = t # lc and lm of g   
    themon = m/lm(f) 
    thectx = ctx(A)
    thecoeff = opp(mul(c,inv(lc(f),thectx),thectx),thectx)

    addmul_geobucket!(g, thecoeff, themon, f.op, A, mod_der = Val(B))
    return g 
end




function delete_op_with_T!(v :: Vector{SigPair{N, C, M}},:: OreAlg) where {N,C,M}
    todelete = Int[]
    for (i,g) in enumerate(v) 
        if mon(g,1)[N] > 0 
            push!(todelete,i)
        end
    end
    reverse!(todelete)
    for s in todelete
        deleteat!(v,s)
    end
end

function sigbasis(gen :: Vector{OrePoly{K,M}}, A :: OreAlg, param :: F5Param) where {K,M} 
    N = nvars(A)
    sgb = prebasis(copy(gen),A)

    sigorder = TOP{typeof(order(A))}(order(A))
    Q = SortedSet{Signature{N}}(sigorder,Signature{N}[]) # signature queue

    for i in 1:length(sgb)
        updateQ!(Q, sgb, sgb[i], sigorder,A)
    end

    if use_geobucket(param) 
        geob = GeoBucket(zero(A))
        tmp_poly = ReuseOrePoly(1,A)
    end
    # main loop 
    while length(Q) > 0 
        if stophol(param) && isholonomic(sgb,A)
            debug(param) && @debug "The current partial gb generates a holonomic ideal"
            delete_op_with_T!(sgb,A)
        end
        debug(param) && @debug "size of the signature queue: $(length(Q))"

        sig = pop!(Q)
        if stat(param)
            globalstats.counters[:f5_candidate_signatures] += 1
        end
        # @debug "dealing with signature: $((sig.ind,sig.mon.exp))"
        if use_geobucket(param)
            f = selectreductant(geob,sgb,sig,A)

        else
            f = selectreductant(sgb,sig,A)
        end
        
        iszero(f) && continue  # is that really supposed to happen ? yes if stophol = true 

        lmf = lm(f)
        if use_geobucket(param)
            f  = sigdiv!(geob,tmp_poly, f, sgb,sigorder, A,param)
        else 
            f  = sigdiv!(f, sgb,sigorder, A,param)
        end
        if iszero(f)
            debug(param) && @debug "reduced to zero"  
        elseif lm(f) != lmf
            updateQ!(Q, sgb, f, sigorder,A) 
            push!(sgb, f)
            debug(param) && @debug "new operator in the gb with leading monomial $(mon(f,1).exp)"
        else
            debug(param) && @debug "no reducible operator found for the given signature"
        end
    end
    return sgb
end


function sgbtogb(sgb :: Vector{SigPair{N,T,M}}, A :: OreAlg,param :: F5Param) where {N,T,M}
    gb = OrePoly{T,M}[sgb[i].op for i in 1:length(sgb)]
    reducebasis!(gb, A,param = param)
    return gb
end




"""
    f5(gens :: Vector{OrePoly{T,M}}, A :: Alg)

Return a reduced Gröbner basis of the left ideal generated by gens for the monomial order defined in A.
"""
function f5(gen :: Vector{OrePoly{K,M}}, A :: OreAlg; param :: F5Param = f5_param()) where {K,M}
    reducebasis!(gen,A,param=param) 
    sgb = sigbasis(gen,A,param)
    return sgbtogb(sgb,A,param)
end

function f5(A :: OreAlg,gen :: Vector{OrePoly{K,M}}; param :: F5Param = f5_param()) where {K,M} 
    return f5(gen,A,param = param)
end


