struct RepInIntMod{T,M}
    echelon :: Vector{OrePoly{T,M}}
    next_incr :: Vector{OrePoly{T,M}}
    g1 :: Vector{OrePoly{T,M}}
end

Base.show(io :: IO, ::RepInIntMod) = print(io,"RepInIntMod")

"""
    representative_in_integral_module_precomp(gb :: Vector{OrePoly{T,M}},sigma :: Int ,A :: OreAlg) where {T,M}

Precomputation step for the representative_in_integral_module and irreducible_monomials functions. 
"""
function representative_in_integral_module_precomp(gb :: Vector{OrePoly{T,M}},sigma :: Int ,A :: OreAlg) where {T,M}
    g1, g2 = separate(gb, A)
    return RepInIntMod(GD_prereduction_init(g2, g1,sigma, A)...,g1)
end

function representative_in_integral_module_precomp_incr!(precom :: RepInIntMod, l :: Int,A :: OreAlg)
    for i in l 
        tmp = GD_prereduction_increment!(precom.echelon, precomp.next_incr, precomp.g1, SortedSet{eltype_mo(A)}(order(A)), A)
        empty!(precomp.next_incr)
        append!(precomp.next_incr, tmp)
    end
end

"""
    representative_in_integral_module(precomp :: RepInIntMod,pol :: OrePoly, A :: OreAlg)
Compute a smaller representative of the operator pol w.r.t. the order of A in the integral of the module defined by precomp 
"""
function representative_in_integral_module(precomp :: RepInIntMod,pol :: OrePoly, A :: OreAlg)
    npol = deepcopy(pol)
    npol = GD_reduction1!(npol,precomp.g1,A)
    npol = reduce_with_echelon!(precomp.echelon,npol,A)
    return npol
end



function GD_prereduction_init(G2 :: Vector{OrePoly{T,M}}, G1 :: Vector{OrePoly{T,M}},sigma :: Int, A :: OreAlg) where {T,M}
    if length(G2) == 0 
        return OrePoly{T,M}[], OrePoly{T,M}[]
    end 
    s = maximum([maximum([sum(m[i] for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv) - sum(m[i] for i in A.nrdv+1:A.nrdv+A.npdv) for m in mons(p)]) for p in G2])
    if  sigma <= s 
        @warn "parameter sigma might be chosen too small" 
    end
    echelon = OrePoly{T,M}[] 
    G = OrePoly{T,M}[]
    donotadd = SortedSet{M}(order(A))

    for g in G2 
        tmpG = OrePoly{T,M}[]
        push!(tmpG,g)
        push!(donotadd,mon(g,1))
        newrel = deepcopy(g)
        newrel = GD_reduction1!(newrel,G1,A)
        newrel = reduce_with_echelon!(echelon,newrel,A)
        add_echelon!(echelon,newrel,A)
        k = maximum([sum(m[i] for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv) - sum(m[i] for i in A.nrdv+1:A.nrdv+A.npdv) for m in mons(g)])
        # println("g")
        # prettyprint(g,A)
        # println("$(k) $(s)")
        for i in k+1:sigma
            tmpG = GD_prereduction_increment!(echelon,tmpG,G1,donotadd,A)
        end
        append!(G,tmpG)
    end
    return echelon, G
end

function GD_prereduction_increment!(echelon :: Vector{OrePoly{T,M}}, G :: Vector{OrePoly{T,M}},G1 :: Vector{OrePoly{T,M}},donotadd :: SortedSet{M}, A :: OreAlg; newlm = nothing) where {T,M}
    newG = OrePoly{T,M}[]
    for g in G
        for i in A.nrdv + A.npdv +1:A.nrdv + 2*A.npdv  
            new_lmon = mon(g,1) * makemon(i,A) 

            # first criterion 
            if new_lmon in donotadd
                @label next_it
                continue 
            end

            # second criterion 
            for h in G1 
                if divide(mon(h,1),new_lmon,A)
                    @goto next_it
                end
            end

            newrel = mul(makepoly(one(ctx(A)), makemon(i,A)),g,A)
            push!(newG,deepcopy(newrel))
            push!(donotadd, new_lmon)

            newrel = GD_reduction1!(newrel,G1,A)
            newrel = reduce_with_echelon!(echelon,newrel,A)
            add_echelon!(echelon,newrel,A)
            if !isnothing(newlm) && length(newrel) > 0 
                push!(newlm, mon(newrel,1))
            end
        end 
    end
    return newG
end

function reduce_with_echelon!(echelon :: Vector{OrePoly{T,M}}, P :: OrePoly{T,M}, A :: OreAlg; augmented = false, echelonvect = Vector{T}[], debug = false, unitary = true) where {T,M}
    res = P
    if augmented 
        vect = T[zero(ctx(A)) for i in 1:length(echelon)]
        push!(vect,one(ctx(A)))
    end

    r = 1 
    e = 1 
    while r <= length(res)
        c,m = res[r]
        div = false
        for i in e:length(echelon)
            if lt(order(A), mon(echelon[i],1), m)
                e = i 
                break
            elseif m == mon(echelon[i],1)
                div = true
                if !unitary
                    mul!(coeff(echelon[i],1),res,A)
                end
                res = sub!(res,mul(c,echelon[i],A),A)
                if augmented 
                    if !unitary 
                        for l in 1:length(echelon)+1
                            vect[l] = mul(coeff(echelon[i],1),vect[l],ctx(A))
                        end
                    end
                    for l in 1:length(echelonvect[i])
                        vect[l] = sub(vect[l],mul(c,echelonvect[i][l], ctx(A)),ctx(A))
                    end
                end
                e = i + 1
                break
            end
        end
        if !div
            r += 1
        end
    end
    if augmented
        return res, vect 
    end
    return res
end

function add_echelon!(echelon :: Vector{OrePoly{T,M}}, P :: OrePoly{T,M}, A :: OreAlg; augmented = false, echelonvect = Vector{T}[], vect = T[],unitary = true ) where {T,M}
    if length(P) == 0
        return 
    end
    if augmented
        if unitary
            for i in 1:length(vect)
                vect[i] = mul(vect[i], inv(coeff(P,1),ctx(A)),ctx(A))
            end
        end

        for v in echelonvect 
            push!(v, zero(ctx(A)))
        end
    end
    if unitary
        makemonic!(P,A)
    end
    for i in length(echelon):-1:1 
        if lt(order(A),mon(P,1),mon(echelon[i],1))
            insert!(echelon,i+1,P)
            if augmented 
                insert!(echelonvect,i+1,vect)
            end
            return
        end 
    end 
    insert!(echelon,1,P)
    if augmented
        insert!(echelonvect,1,vect)
    end
end


function mod_derivatives!(pol :: OrePoly{T,M},A :: OreAlg) where {T,M}
    toremove = Int[]
    for (i,(c,m)) in enumerate(pol) 
        for j in A.nrdv + 1: A.nrdv + A.npdv 
            if m[j] != 0 
                push!(toremove,i)
                break
            end
        end
    end
    deleteat!(coeffs(pol),toremove)
    deleteat!(mons(pol),toremove)
end




function GD_reduction1!(pol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}},A :: OreAlg) where {T,M}
    mod_derivatives!(pol,A)

    r=1 
    while r <= length(pol) 
        div = false
        for i in 1:length(gb)
            if divide(mon(gb[i],1), mon(pol,r),A)
                div = true
                pol = reduce!(pol,r, gb[i], A)
                mod_derivatives!(pol,A)
                break
            end
        end
        if !div 
            r = r + 1
        end
    end
    return pol
end

function separate(gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T, M}
    g1 = OrePoly{T,M}[]
    g2 = OrePoly{T,M}[]
    for g in gb 
        deriv = false
        m = mon(g,1)
        for i in A.nrdv + 1:A.nrdv + A.npdv 
            if m[i] > 0 
                deriv = true 
                break 
            end
        end
        if deriv 
            push!(g2,g)
        else 
            push!(g1,g)
        end
    end
    return g1, g2 
end

"""
    irreducible_monomials(precomp ::RepInIntMod, l :: Int, A :: OreAlg)
    
Returns every monomial of degree up to l that do not have smaller representatives w.r.t the order of A 
in the integral of the module defined in precomp. 

Warning: If the parameter sigma chosen for precomp is not large enough, this function may return some reducible monomials.
"""
function irreducible_monomials(precomp ::RepInIntMod, l :: Int, A :: OreAlg)
    m = makemon(-1,A)
    T = eltype_co(A)
    M = eltype_mo(A)
    understair = SortedSet{M}(order(A),[m])
    lm_g1 = [mon(p,1) for p in precomp.g1]
    lm_ech = SortedSet{M}(order(A),[mon(p,1) for p in precomp.echelon])
    irred = OrePoly{T,M}[]

    if m in lm_g1 
        return irred 
    end
    if !(m in lm_ech)
        push!(irred, OrePoly([one(ctx(A))],[m]))
    end

    for i in 1:l
        understair, nirred = next_slice(understair,lm_g1,lm_ech,A)
        for m in nirred 
            push!(irred,OrePoly([one(ctx(A))],[m]))
        end
    end
    return irred
end

function next_slice(understair :: SortedSet{M}, lm_g1 :: Vector{M}, lm_ech :: SortedSet{M}  , A :: OreAlg) where {M}
    nunderstair = SortedSet{M}(order(A))
    irred = M[]
    for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv 
        x = makemon(i,A)
        for m in understair
            prod = x*m
            
            if prod in nunderstair || any(divide(p,prod,A) for p in lm_g1) 
                continue 
            end
            push!(nunderstair,prod)
            if !(prod in lm_ech)
                push!(irred,prod)
            end
        end
    end 
    return nunderstair,irred
end