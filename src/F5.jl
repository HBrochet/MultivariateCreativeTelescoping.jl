
function prebasis(gen :: Vector{OrePoly{T,M}}) where {T, M}
    N = nvars(M)
    sgb = SigPair{N, T, M}[]
    for i in 1:length(gen)
        sig = Signature([0 for j in 1:N],i)
        push!(sgb, SigPair(gen[i], sig))
    end
    return sgb
end



function updateQ!(Q :: SortedSet, sgb :: Vector{SigPair{N,T,M}}, 
                    g :: SigPair{N,T,M}, sorder :: SigOrder, A :: OreAlg) where {N, T, M}
    for i in 1:length(sgb)
        if !iscompatible(sgb[i],g,A)
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
        if divide(sgb[i].sig,sigma,A)
            mon = sigma.mon/sgb[i].sig.mon
            return mul(mon, sgb[i], A)
        end
    end
    # prettyprint(sigma,A)
    # prettyprint(sgb,A)
end


function sigdiv!(redf :: SigPair, sgb :: Vector{SigPair{N,T,M}}, so :: SigOrder, A :: OreAlg) where {N,T,M}
    res = redf
    r = 1
    divhappened = false
    while r <= length(res.op)
        div = false
        for i in length(sgb):-1:1
            if divide(mon(sgb[i],1), mon(res,r),A) && lt(so,mul(mon(res,r)/mon(sgb[i],1), sgb[i].sig), res.sig)
                div = true
                res = reduce!(res,r, sgb[i], A)
                globalstats.counters[:f5_number_divisions]+=1
                break
            end
        end
        if !div 
            r = r + 1
        else
            divhappened = true
        end

    end
    return res, divhappened  
end

# Reduce mon r of pol f with pol g 
function reduce!(f :: SigPair, r :: Int, g :: SigPair, A :: OreAlg) 
    themon = mon(f,r)/mon(g,1)
    thectx = ctx(A)
    thecoeff = opp(mul(coeff(f,r),inv(coeff(g,1),thectx),thectx),thectx)
    np = OrePoly([thecoeff],[themon])
    np=mul(np, g.op, A)
    globalstats.counters[:f5_field_operations] += length(np)
    if ctx(A) isa RatFunCtx
        globalstats.counters[:f5_size_coeff] += length(numerator(thecoeff,false)) + length(Nemo.denominator(thecoeff,false))
    end
    return SigPair(add!(f.op, np, A),f.sig)
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

function sigbasis(gen :: Vector{OrePoly{K,M}}, A :: OreAlg; stophol :: Bool =false) where {K,M} 
    N = nvars(A)
    sgb = prebasis(deepcopy(gen))

    sigorder = TOP{typeof(order(A))}(order(A))
    Q = SortedSet{Signature{N}}(sigorder,Signature{N}[]) # signature queue


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
        @debug "dealing with signature: $(sig)"
        f = selectreductant(sgb,sig,A)

        if isnothing(f)
            globalstats.counters[:f5_eliminated_signatures_stophol] += 1
            continue
        end
        

        if length(f) == 0 
            println("zero")
        end
        f, divhappened = sigdiv!(f, sgb,sigorder, A)

        if divhappened
            if length(f.op) > 0
                updateQ!(Q, sgb, f, sigorder,A) 
                push!(sgb, f)
                @debug "new operator in the gb with leading monomial $(mon(f,1))"
            else
                @debug "reduced to zero"
            end
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
function f5(gen :: Vector{OrePoly{K,M}}, A :: OreAlg; stophol :: Bool = false) where {K,M} 
    sgb = sigbasis(gen,A,stophol=stophol)
    return sgbtogb(sgb,A)
end

