function separate(gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T, M}
    g1 = OrePoly{T,M}[]
    g2 = OrePoly{T,M}[]
    for g in gb 
        deriv = false
        m = mon(g,1)
        for i in 2:1 + A.npdv 
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

function GD_hol_prereduction_init(G2 :: Vector{OrePoly{T,M}}, G1 :: Vector{OrePoly{T,M}},sigma :: Int, A :: OreAlg) where {T,M}
    if length(G2) == 0 
        return OrePoly{T,M}[], OrePoly{T,M}[]
    end

    s= maximum([sum(mon(g,1)) for g in G2]) + sigma # todo best strategy ? 
    gb = OrePoly{T,M}[] 
    G = OrePoly{T,M}[]
    donotadd = SortedSet{M}(order(A))
    for g in G2 
        tmpG = OrePoly{T,M}[]
        push!(tmpG,g)
        push!(donotadd,mon(g,1))
        newrel = copy(g)
        newrel = GD_reduction1!(newrel,G1,A)
        if !isempty(newrel)
            push!(gb,newrel)
        end

        for i in sum(mon(g,1))+1:s
            tmpG = GD_hol_prereduction_increment!(gb,tmpG,G1,donotadd,A)
        end
        append!(G,tmpG)
    end
    for p in gb
        if length(p) == 0 
            error("zero vector in basis in ggd")
        end
    end

    append!(A.nomul,[A.nrdv + A.npdv + i for i in 1:A.npdv])
    res = f5(gb,A)
    resize!(A.nomul,length(A.nomul)-A.npdv)
    return res, G
end

function GD_hol_prereduction_increment!(gb :: Vector{OrePoly{T,M}}, G :: Vector{OrePoly{T,M}},G1 :: Vector{OrePoly{T,M}},donotadd :: SortedSet{M}, A :: OreAlg) where {T,M}
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
            if !isempty(newrel)
                push!(gb,newrel)
            end
        end 
    end
    return newG 
end

function reduce_with_echelon!(echelon :: Vector{OrePoly{T,M}}, P :: OrePoly{T,M}, A :: OreAlg; augmented = false, echelonvect = Vector{T}[], debug = false) where {T,M}
    if augmented 
        vect = T[zero(ctx(A)) for i in 1:length(echelon)]
        push!(vect,one(ctx(A)))
    end

    r = 1 
    e = 1 
    if debug
        println("reducing")
        prettyprint(P,A)
        # println("echelon")
        # prettyprint(echelon,A)
    end
    while r <= length(P)
        c,m = P[r]
        div = false
        if debug
            print("reduction mon ")
            printmon(m,A)
            println()
            println(e," ",length(echelon))
        end
        for i in e:length(echelon)
            if debug 
                print("by ")
                printmon(mon(echelon[i],1), A)
                println(" ? ")
            end
            if lt(order(A), mon(echelon[i],1), m)
                e = i 
                break
            elseif m == mon(echelon[i],1)
                div = true
                @assert coeff(echelon[i],1) == one(ctx(A))
                sub!(P,mul(c,echelon[i],A),A)
                if debug 
                    println("yes, result :")
                    prettyprint(P,A)
                end
                if augmented 
                    for l in 1:length(echelonvect[i])
                        vect[l] = vect[l] - c * echelonvect[i][l]
                    end
                end
                e = i + 1
                break
            end
        end
        if !div
            r += 1
        end
        # if debug 
        #     error("fini 0")
        # end
    end
    if augmented
        return vect 
    end
end

function find_stable_mon_set_hol(starting_pol :: OrePoly{T,M}, G1 :: Vector{OrePoly{T,M}},echelon :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T, M}
    im = Dict{M,OrePoly{T,M}}()
    toadd = SortedSet{M}(order(A),mons(starting_pol))
    while !isempty(toadd)
        m = pop!(toadd)
        dtm = mul(makemon(1,A),makepoly(one(ctx(A)),m),A)
        dtm = GD_reduction1!(dtm, G1, A)
        reduce_with_echelon!(echelon, dtm, A, debug = false)
        im[m] = dtm 
        for m in mons(dtm)
            if !haskey(im,m)
                push!(toadd,m)
            end
        end 

    end
    return im 
end

function find_cbl(map_ :: Dict{M,OrePoly{T,M}}, starting_pol :: OrePoly{T,M}, A :: OreAlg, finalA :: OreAlg) where {T,M}
    nrel = starting_pol
    @debug "looking for the LDE"
    echelon_derivatives = OrePoly{T,M}[]
    echelonvect = Vector{T}[]
    ord = 0
    while true  
        @debug "dealing with $(ord)th derivative"
        nrel_red = copy(nrel)
        v = reduce_with_echelon!(echelon_derivatives,nrel_red,A,augmented = true,echelonvect = echelonvect)
        if length(nrel_red) == 0 
            mons = [makemon(1,finalA)^(i-1) for i in 1:length(v)]
            return OrePoly(v,mons)
        end
        add_echelon!(echelon_derivatives, nrel_red, A, augmented = true, echelonvect = echelonvect, vect = v)
        nrel = compute_next_rel(map_, nrel,A)
        ord += 1 
    end
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

function add_echelon!(echelon :: Vector{OrePoly{T,M}}, P :: OrePoly{T,M}, A :: OreAlg; augmented = false, echelonvect = Vector{T}[], vect = T[] ) where {T,M}
    if length(P) == 0
        return 
    end
    if augmented 
        for i in 1:length(vect)
            vect[i] = mul(vect[i], inv(coeff(P,1),ctx(A)),ctx(A))
        end

        for v in echelonvect 
            push!(v, zero(ctx(A)))
        end
    end
    makemonic!(P,A)
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

function compute_next_rel(map_ :: Dict{M,OrePoly{T,M}}, pol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    cs = copy(coeffs(pol))
    ms = copy(mons(pol))

    for i in 1:length(cs)
        cs[i] = derivative(cs[i],ctx(A).vars[1])
    end
    res = OrePoly(cs,ms)
    for (c,m) in pol 
        res = add!(res,mul(c,map_[m],A) ,A)
    end
    return res 
end

function MCT_hol(starting_pol :: OrePoly, g :: Vector{OrePoly{T,M}},sigma :: Int, A :: OreAlg, finalA :: OreAlg) where {T,M}
    sp = deepcopy(starting_pol)
    if !isholonomic(g,A)
        error("input does not define an holonomic ideal or is not a GB")
    end
    g1, g2 = separate(g, A)

    echelon, next_incr = GD_hol_prereduction_init(g2, g1, sigma, A)

    GD_reduction1!(sp,g1,A)
    reduce_with_echelon!(echelon,sp,A)
    map_ = find_stable_mon_set_hol(sp, g1, echelon, A)
    res = find_cbl(map_, sp, A, finalA)
    normalize!(res,finalA)
    den = ctx(finalA).F(denominator(res, finalA))
    mul!(den,res,finalA)
    return res
end