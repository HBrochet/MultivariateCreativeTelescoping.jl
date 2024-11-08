function prebasis(gen :: Vector{OrePoly{T,M}}) where {T, M}
    N = nvars(M)
    sgb = SigPair{N, T, M}[]
    for i in 1:length(gen)
        lmon = lm(gen[i])
        sig = Signature([lmon[i] for j in 1:N],i)
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
    return zero(A)
end

# version without GeoBucket
function sigdiv!(redf :: SigPair{N,K,M}, sgb :: Vector{SigPair{N,K,M}}, so :: SigOrder, A :: OreAlg, :: Val{false}) where {N,K,M}
    res = redf
    r = 1
    divhappened = false
    ctr = 0
    while r <= length(res.op)
        div = false
        for i in length(sgb):-1:1
            if iscompatible(mon(sgb[i],1), mon(res,r),A) && divide(mon(sgb[i],1), mon(res,r),A) && lt(so,mul(mon(res,r)/mon(sgb[i],1), sgb[i].sig), res.sig)
                div = true
                res = reduce!(res,r, sgb[i], A)
                globalstats.counters[:f5_number_divisions]+=1
                break
            end
        end
        if !div 
            break # to be removed 
            r = r + 1
        else
            divhappened = true
        end

    end
    return res, divhappened  
end

# version without GeoBucket
# Reduce mon r of pol f with pol g 
function reduce!(f :: SigPair, r :: Int, g :: SigPair, A :: OreAlg) 
    themon = mon(f,r)/mon(g,1)
    thectx = ctx(A)
    thecoeff = opp(mul(coeff(f,r),inv(coeff(g,1),thectx),thectx),thectx)
    np = OrePoly([thecoeff],[themon])
    np=mul(np, g.op, A)
    globalstats.counters[:f5_field_operations] += length(np)
    globalstats.counters[:f5_size_m]+= sum(themon)
    if is_polynomial(themon,A)
        globalstats.counters[:f5_m_is_pol] += 1
    end
    if ctx(A) isa RatFunCtx
        globalstats.counters[:f5_size_coeff] += length(numerator(thecoeff,false)) + length(Nemo.denominator(thecoeff,false))
    end
    return SigPair(add!(f.op, np, A),f.sig)
end

# version with GeoBucket
function sigdiv!(redf :: SigPair{N,K,M}, sgb :: Vector{SigPair{N,K,M}}, so :: SigOrder, A :: OreAlg, :: Val{true}) where {N,K,M,}
    sig = redf.sig
    res = GeoBucket(redf.op)
    divhappened = false
    while true 
        div = false
        lco,lmon = lt(res,A)
        if iszero(lco)
            break 
        end
        for i in length(sgb):-1:1
            lmon2 = lm(sgb[i])
            if iscompatible(lmon2, lmon,A) && divide(lmon2, lmon,A) && lt(so,mul(lmon/lmon2, sgb[i].sig), sig)
                div = true
                res = reduce!(res, (lco,lmon), sgb[i], A)
                globalstats.counters[:f5_number_divisions]+=1
                break
            end
        end
        if !div
            break 
        else 
            divhappened = true 
            # @assert iszero(lt(res,A)[1]) || lmon != lm(res,A)
        end
    end
    return SigPair(normalform(res,A),sig), divhappened 
end

# version SigPairGB
function reduce!(g :: GeoBucket, t:: Tuple{K,M}, f :: SigPair, A :: OreAlg)  where {K, M}
    @assert nozero(g,A)
    (c,m) = t # lc and lm of g   
    themon = m/lm(f) 
    thectx = ctx(A)
    thecoeff = opp(mul(c,inv(lc(f),thectx),thectx),thectx)


    # cp = deepcopy(g) 
    addmul_geobucket!(g, thecoeff, themon, f.op, A)

    # gg = normalform(cp,A) 
    # gg = add!(gg, mul(OrePoly([thecoeff],[themon]),f.op,A),A)

    # cp2 = deepcopy(g)
    # gg2 = normalform(cp2,A)
    
    # @assert nozero(g,A)
    # @assert gg == gg2
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

function sigbasis(gen :: Vector{OrePoly{K,M}}, sigorderQ :: SigOrder, A :: OreAlg; stophol :: Bool =false, geobucket ::Val{B} = Val{false}) where {K,M, B} 
    N = nvars(A)
    sgb = prebasis(deepcopy(gen))

    Q = SortedSet{Signature{N}}(sigorderQ,Signature{N}[]) # signature queue
    sigorder = TOP{typeof(order(A))}(order(A))

    for i in 1:length(sgb)
        updateQ!(Q, sgb, sgb[i], sigorder,A)
    end

    # main loop 
    while length(Q) > 0 
        if stophol && isholonomic(sgb,A)
            delete_op_with_T!(sgb,A)
        end
        @debug "size of the signature queue: $(length(Q))"
        sig = pop!(Q)
        globalstats.counters[:f5_candidate_signatures] += 1
        @debug "dealing with signature: $((sig.ind,sig.mon.exp))"
        f = selectreductant(sgb,sig,A)
        if iszero(f) 
            globalstats.counters[:f5_eliminated_signatures_stophol] += 1
            continue
        end
        
        lmf = lm(f)
        f, divhappened = sigdiv!(f, sgb,sigorder, A, Val(B))

        if iszero(f)
            @debug "reduced to zero"  
        elseif lm(f) != lmf
            updateQ!(Q, sgb, f, sigorder,A) 
            push!(sgb, f)
            @debug "new operator in the gb with leading monomial $(mon(f,1))"
        else
            @debug "no reducible operator found for the given signature"
        end
    end
    return sgb
end


function sgbtogb(sgb :: Vector{SigPair{N,T,M}}, A :: OreAlg) where {N,T,M}
    gb = OrePoly{T,M}[sgb[i].op for i in 1:length(sgb)]
    reducebasis!(gb, A)
    return gb
end


"""
    f5(gens :: Vector{OrePoly{T,M}}, A :: Alg)

Return a reduced Gröbner basis of the left ideal generated by gens for the monomial order defined in A.
"""
function f5(gen :: Vector{OrePoly{K,M}}, A :: OreAlg; stophol :: Bool = false, geobucket ::Val{B} = Val(false),sigorder = nothing) where {K,M,B} 
    if isnothing(sigorder)
        sigorder = TOP{typeof(order(A))}(order(A))
    end
    sgb = sigbasis(gen,sigorder,A,stophol=stophol,geobucket = geobucket)
    return sgbtogb(sgb,A)
end

function f5(A :: OreAlg,gen :: Vector{OrePoly{K,M}}; stophol :: Bool = false) where {K,M} 
    return f5(gen,A)
end


