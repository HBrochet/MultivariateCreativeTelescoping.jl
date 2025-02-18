

# duplicate using geobucket
function GD_reduction1!(p :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}},A :: OreAlg, geob :: GeoBucket,tmp_poly ::ReuseOrePoly) where {T,M}
    init_GeoBucket!(geob,p)
    return GD_reduction1!(gb,A,geob,tmp_poly)
end

function GD_reduction1!(gb :: Vector{OrePoly{T,M}},A :: OreAlg, geob :: GeoBucket,tmp_poly ::ReuseOrePoly) where {T,M}
    mod_derivatives!(geob,A)
    while true 
        div = false
        lco,lmon = lt(geob,A)
        if iszero(lco)
            break
        end
        for i in length(gb):-1:1
            lmon2 = lm(gb[i])
            if iscompatible(lmon2, lmon,A) && divide(lmon2, lmon,A)
                div = true
                geob = reduce_geob!(geob, (lco,lmon), gb[i], A, Val(true))
                break
            end
        end
        if !div
            push!(tmp_poly,lco,lmon)
            rem_lt!(geob,lmon)
        end
    end
    return copy_to_OrePoly!(tmp_poly,A)
end


function mod_derivatives!(geob :: GeoBucket, A :: OreAlg)
    for i in 1:length(geob)
        mod_derivatives!(geob,i,A)
    end
    return geob 
end
        
function mod_derivatives!(geob :: GeoBucket,i :: Int, A :: OreAlg)
    tab = active_tab(geob,i) 
    len = geob.indices[i] 
    r = 1
    w = 1 
    while r <= len  
        m = mon(tab,r)
        if any(m[i] > 0 for i in A.nrdv+1:A.nrdv+A.npdv)
            r += 1 
            break 
        end
        w += 1
        r += 1 
    end

    while r <= len 
        m = mon(tab,r)
        if all(m[i] == 0 for i in A.nrdv+1:A.nrdv+A.npdv)
            c = coeff(tab,r)
            tab[w] = (c,m)
            w += 1 
        end
        r += 1 
    end
    geob.indices[i] = w - 1
    return geob
end




function GD_prereduction_init(G2 :: Vector{OrePoly{T,M}}, G1 :: Vector{OrePoly{T,M}},sigma :: Int, A :: OreAlg, geob :: GeoBucket,tmp_poly ::ReuseOrePoly;tracer :: Val{B} =  Val(false)) where {T,M,B}
    if B 
        set_trace = SortedSet{M}(order(A))
    end
    
    if length(G2) == 0 
        if B
            return OrePoly{T,M}[], OrePoly{T,M}[], set_trace
        else 
            return OrePoly{T,M}[], OrePoly{T,M}[]
        end
    end 
    donotadd = SortedSet{M}(order(A))

    echelon = OrePoly{T,M}[] 
    G = OrePoly{T,M}[]

    for g in G2 
        tmpG = OrePoly{T,M}[]
        push!(tmpG,g)
        m = mon(g,1)

        push!(donotadd,m)
        newrel = GD_reduction1!(g,G1,A,geob,tmp_poly)
        newrel = reduce_with_echelon!(echelon,newrel,A,geob,tmp_poly)
        if B && length(newrel) == 0 
            push!(set_trace,m)
        end
        add_echelon!(echelon,newrel,A)
        k = sum(m[i] for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv)
        for i in k+1:sigma
            if B
                tmpG = GD_prereduction_increment!(echelon,tmpG,G1,donotadd,A,geob,tmp_poly,set_trace=set_trace)
            else 
                tmpG = GD_prereduction_increment!(echelon,tmpG,G1,donotadd,A,geob,tmp_poly)
            end
        end
        append!(G,tmpG)
    end
    if B
        return echelon, G, set_trace
    else 
        return echelon, G
    end
end


function GD_prereduction_init_apply(G2 :: Vector{OrePoly{T,M}}, G1 :: Vector{OrePoly{T,M}},sigma :: Int,_donotadd ::SortedSet{M}, A :: OreAlg, geob :: GeoBucket,tmp_poly ::ReuseOrePoly) where {T,M}
    if length(G2) == 0 
        return OrePoly{T,M}[], OrePoly{T,M}[]
    end 
    donotadd = deepcopy(_donotadd)
    donotadd = SortedSet{M}(order(A))

    echelon = OrePoly{T,M}[] 
    G = OrePoly{T,M}[]

    for g in G2 
        tmpG = OrePoly{T,M}[]
        push!(tmpG,g)
        m = mon(g,1)
        if !(m in donotadd)
            push!(donotadd,mon(g,1))
            newrel = GD_reduction1!(g,G1,A,geob,tmp_poly)
            newrel = reduce_with_echelon!(echelon,newrel,A,geob,tmp_poly)
            add_echelon!(echelon,newrel,A)
        end
        k = sum(m[i] for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv)
        for i in k+1:sigma
            tmpG = GD_prereduction_increment!(echelon,tmpG,G1,donotadd,A,geob,tmp_poly)
        end
        append!(G,tmpG)
    end
    return echelon, G
end

function GD_prereduction_increment!(echelon :: Vector{OrePoly{T,M}}, G :: Vector{OrePoly{T,M}},G1 :: Vector{OrePoly{T,M}},donotadd :: SortedSet{M}, A :: OreAlg, geob :: GeoBucket,tmp_poly ::ReuseOrePoly;set_trace = nothing) where {T,M}
    newG = OrePoly{T,M}[]
    for g in G
        for i in A.nrdv + A.npdv +1:A.nrdv + 2*A.npdv  
            m = makemon(i,A) 
            new_lmon = mon(g,1) * m

            # first criterion 
            if new_lmon in donotadd
                @label next_it
                continue 
            end

            # second criterion 
            for h in G1 
                if divide(mon(h,1),new_lmon)
                    @goto next_it
                end
            end

            newrel = addmul_geobucket!(geob,one(ctx(A)),m,g,A)
            newrel = normalform(geob,A)

            push!(newG,newrel)

            push!(donotadd, new_lmon)


            red = GD_reduction1!(newrel,G1,A,geob,tmp_poly)
            red = reduce_with_echelon!(echelon,red,A,geob,tmp_poly)
            if !isnothing(set_trace) && length(red) == 0 
                push!(set_trace,new_lmon)
            end
            add_echelon!(echelon,red,A)
        end 
    end
    return newG
end