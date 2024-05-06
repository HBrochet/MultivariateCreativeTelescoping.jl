function maxdegx(p :: OrePoly, A:: OreAlg)
    return maximum(sum(m[i] for i in 2+A.npdv:1+2*A.npdv) for m in mons(p))
end
function maxdegx(m :: OreMonVE, A:: OreAlg)
    return sum(m[i] for i in 2+A.npdv:1+2*A.npdv)
end


function update_stairs!(stairs :: Vector{SortedSet{M}}, v :: SortedSet{M}, d :: Int, A:: OreAlg) where M 
    for m in v 
        deg = maxdegx(m,A)
        if deg > d
            delete!(stairs[deg+1],m)
        end
    end
end 

function termination_criterion_der_red_map_precomp(stairs :: Vector{SortedSet{M}},d :: Int ,l :: Int) where M
    k = 0
    for i in d+1:length(stairs)
        if length(stairs[i]) == 0 
            k +=1 
        else
            k = 0 
        end
        if k == l 
            return true
        end
    end
    return false
end 


function der_red_map_precomp(spol :: OrePoly, gb :: Vector{OrePoly{T,M}},A :: OreAlg) where {T,M}
    red_dt = find_red_dt(gb,A)
    g1,g2 = separate(gb,A)
    lm_g1 = [mon(g,1) for g in g1]
    spol = GD_reduction1!(spol,g1,A)
    d = maxdegx(spol,A)
    if length(g2) > 0
        tmp =  maximum([maximum([sum(m[i] for i in 2+A.npdv:1+2*A.npdv) - sum(m[i] for i in 2:1+A.npdv) for m in mons(p)]) for p in g2])
    else
        tmp = 0
    end
    d = max(d,tmp) + 2
    echelon, next_incr = GD_prereduction_init(g2, g1, d, A)
    lm_ech = SortedSet(order(A),[mon(g,1) for g in echelon])
    spol = reduce_with_echelon!(echelon,spol,A)
    nd = maxdegx(spol,A)
    stairs = SortedSet{M}[]
    # maximal increase in degree when multiplying by dt
    l = maximum([sum(mon(red_dt,j)[i] for i in 2+A.npdv:1+2*A.npdv) - sum(mon(red_dt,j)[i] for i in 2:1+A.npdv) for j in 2:length(red_dt)])
    if l <= 0 
        for i in 1:5
            next_incr = GD_prereduction_increment!(echelon,next_incr,g1,SortedSet{M}(order(A)),A)
        end
        spol = reduce_with_echelon!(echelon,spol,A)
        return spol,g1, red_dt, echelon
    end

    understair = SortedSet{M}(order(A),[makemon(1,A)^i for i in 0:mon(red_dt,1)[1]-1])
    push!(stairs, SortedSet{M}(order(A)))
    for i in 1:nd
        understair, _ = next_slice(understair, lm_g1, lm_ech, A)
        push!(stairs, SortedSet{M}(order(A)))
    end

    for i in nd+1:d
        understair, irred = next_slice(understair, lm_g1, lm_ech, A)
        push!(stairs, SortedSet{M}(order(A),irred))
    end

    itt = 0
    while true
        itt += 1
        newlm = SortedSet{M}(order(A))
        next_incr = GD_prereduction_increment!(echelon,next_incr,g1,SortedSet{M}(order(A)),A;newlm = newlm)
        union!(lm_ech, newlm)
        understair, irred = next_slice(understair, lm_g1, lm_ech, A)
        push!(stairs,SortedSet{M}(order(A),irred))
        update_stairs!(stairs,newlm,nd,A)
        if termination_criterion_der_red_map_precomp(stairs,nd,l)
            break
        end
    end
    ## bonus for minimality 
    for i in 1:5
        next_incr = GD_prereduction_increment!(echelon,next_incr,g1,SortedSet{M}(order(A)),A)
    end
    spol = reduce_with_echelon!(echelon,spol,A)
    return spol, g1, red_dt, echelon
end

function find_red_dt(g :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T,M}
    for (j,p) in enumerate(g) 
        if all(mon(p,1)[i] == 0 for i in 2:nvars(A))
            tmp = p 
            deleteat!(g,j)
            return tmp 
        end
    end
    # prettyprint(g,A)
    error("reducer dt not found")
end

function find_der_red_map(spol :: OrePoly, g1 :: Vector{OrePoly{T,M}}, red_dt :: OrePoly{T,M}, echelon :: Vector{OrePoly{T,M}},A :: OreAlg) where {T, M}
    im = Dict{M,OrePoly{T,M}}()
    toadd = SortedSet{M}(order(A),mons(spol))
    dt = OrePoly([one(ctx(A))], [makemon(1,A)])

    while !isempty(toadd)
        m = pop!(toadd)
        dtm = mul(dt,OrePoly([one(ctx(A))],[m]),A)
        dtm = div!(dtm,[red_dt],A)
        dtm = GD_reduction1!(dtm, g1, A)
        dtm = reduce_with_echelon!(echelon, dtm, A)
        im[m] = dtm 
        for m in mons(dtm)
            if !haskey(im,m)
                push!(toadd,m)
            end
        end 
    end
    return im 
end


function der_red_map(A :: OreAlg, spol :: OrePoly, gb :: Vector{OrePoly{T,M}}) where {T,M}
    spol, g1, red_dt, echelon = der_red_map_precomp(spol,gb,A)
    return find_der_red_map(spol,g1,red_dt,echelon,A), spol
end