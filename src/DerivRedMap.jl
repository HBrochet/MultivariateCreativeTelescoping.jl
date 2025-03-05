struct TracerDerRedMap{M} 
    s :: SortedSet{M}
    i :: Int 
end


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
    if length(g2) == 0
        return spol, g1, red_dt, OrePoly{T,M}[]
    end

    lm_g1 = [mon(g,1) for g in g1]
    spol = GD_reduction1!(spol,g1,A)
    d = maxdegx(spol,A)
    tmp =  maximum([maximum([sum(m[i] for i in 2+A.npdv:1+2*A.npdv) - sum(m[i] for i in 2:1+A.npdv) for m in mons(p)]) for p in g2])
    d = max(d,tmp) + 2
    t = time()
    echelon, next_incr = GD_prereduction_init(g2, g1, d, A)
    lm_ech = SortedSet(order(A),[mon(g,1) for g in echelon])
    spol = reduce_with_echelon!(echelon,spol,A)
    nd = maxdegx(spol,A)
    stairs = SortedSet{M}[]
    # maximal increase in degree when multiplying by dt
    l = maximum([sum(mon(red_dt,j)[i] for i in 2+A.npdv:1+2*A.npdv) - sum(mon(red_dt,j)[i] for i in 2:1+A.npdv) for j in 2:length(red_dt)])
    if l <= 0 
        #is it overkill ? 
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

    t = time()
    for i in nd+1:d
        understair, irred = next_slice(understair, lm_g1, lm_ech, A)
        push!(stairs, SortedSet{M}(order(A),irred))
    end

    itt = 0
    t = time()
    while !termination_criterion_der_red_map_precomp(stairs,nd,l)
        itt += 1
        newlm = SortedSet{M}(order(A))
        next_incr = GD_prereduction_increment!(echelon,next_incr,g1,SortedSet{M}(order(A)),A;newlm = newlm)
        union!(lm_ech, newlm)
        understair, irred = next_slice(understair, lm_g1, lm_ech, A)
        push!(stairs,SortedSet{M}(order(A),irred))
        update_stairs!(stairs,newlm,nd,A)
    end

    spol = reduce_with_echelon!(echelon,spol,A)
    return spol, g1, red_dt, echelon
end


# duplicate using geobucket


function der_red_map_precomp2(spol :: OrePoly, gb :: Vector{OrePoly{T,M}},A :: OreAlg,geob :: GeoBucket,tmp_poly :: ReuseOrePoly;tracer ::Val{B} = Val(false)) where {T,M,B}
    red_dt= find_red_dt(gb,A)    
    L = undefOrePoly(length(red_dt)-1,A)
    for i in 2:length(red_dt)
        L[i-1] = red_dt[i]
    end

    red_dt_degx = maximum(sum(mon(red_dt,j)[i] for i in 2+A.npdv:1+2*A.npdv) - sum(mon(red_dt,j)[i] for i in 2:1+A.npdv) for j in 2:length(red_dt))
    g1,g2 = separate(gb,A)
    if length(g2) == 0
        if B 
            return spol, g1, red_dt, OrePoly{T,M}[], TracerDerRedMap{M}(SortedSet{M}(order(A)),0)
        else
            return spol, g1, red_dt, OrePoly{T,M}[]
        end
    end

    spol_degx = sum(mon(red_dt,1)[i] for i in 2+A.npdv:1+2*A.npdv)
    spol = GD_reduction1!(spol,g1,A,geob,tmp_poly)
    rho =  1
    ll = maximum(sum(mon(p,1)[i] for i in 2+A.npdv:1+2*A.npdv) for p in g2)
    s = max(max(red_dt_degx,0) + max(spol_degx,0),ll)
    if B 
        echelon, next_incr, set_trace = GD_prereduction_init(g2, g1, s + rho, A,geob,tmp_poly,tracer = Val(B))
    else
        echelon, next_incr = GD_prereduction_init(g2, g1, s + rho, A,geob,tmp_poly,tracer = Val(B))
    end
    spol = reduce_with_echelon!(echelon,spol,A,geob,tmp_poly)
    while true
        @label letsgo
        done = SortedSet{M}(order(A))
        todo = SortedSet{M}(order(A))
        append!(todo,mons(spol))
        itt = 0 
        while !isempty(todo)
            m = poplast!(todo) 
            m_degx = sum(m[i] for i in 2+A.npdv:1+2*A.npdv)
            if m_degx > s + rho 
                s += 1 
                next_incr = GD_prereduction_increment!(echelon,next_incr,g1,SortedSet{M}(order(A)),A,geob,tmp_poly)
                spol = reduce_with_echelon!(echelon,spol,A,geob,tmp_poly)
                @goto letsgo
            end
            push!(done,m)
            geob = addmul_geobucket!(geob,one(ctx(A)),m,L,A)
            Lm = GD_reduction1!(g1,A,geob,tmp_poly)
            Lm = reduce_with_echelon!(echelon,Lm,A,geob,tmp_poly)
            for mo in mons(Lm)
                if !(mo in done)
                    push!(todo,mo)
                end
            end
            itt += 1 
            # if itt == 5
            #     prettyprint(echelon,A)
            #     error("fin")
            # end
        end
        if B 
            return  spol, g1, red_dt, echelon, TracerDerRedMap{M}(set_trace,s+rho)
        else
            return  spol, g1, red_dt, echelon
        end 
    end

end




        

function der_red_map_precomp(spol :: OrePoly, gb :: Vector{OrePoly{T,M}},A :: OreAlg,geob :: GeoBucket,tmp_poly :: ReuseOrePoly) where {T,M}
    red_dt = find_red_dt(gb,A)
    g1,g2 = separate(gb,A)
    if length(g2) == 0
        return spol, g1, red_dt, OrePoly{T,M}[]
    end

    lm_g1 = [mon(g,1) for g in g1]
    spol = GD_reduction1!(spol,g1,A)
    # d = maxdegx(spol,A)
    # tmp =  maximum([maximum([sum(m[i] for i in 2+A.npdv:1+2*A.npdv) - sum(m[i] for i in 2:1+A.npdv) for m in mons(p)]) for p in g2])
    d = max(d,tmp) + 2
    echelon, next_incr = GD_prereduction_init(g2, g1, d, A)
    lm_ech = SortedSet(order(A),[mon(g,1) for g in echelon])
    spol = reduce_with_echelon!(echelon,spol,A)
    nd = maxdegx(spol,A)
    stairs = SortedSet{M}[]
    # maximal increase in degree when multiplying by dt
    l = maximum([sum(mon(red_dt,j)[i] for i in 2+A.npdv:1+2*A.npdv) - sum(mon(red_dt,j)[i] for i in 2:1+A.npdv) for j in 2:length(red_dt)])
    if l <= 0 
        #is it overkill ? 
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

    t = time()
    for i in nd+1:d
        understair, irred = next_slice(understair, lm_g1, lm_ech, A)
        push!(stairs, SortedSet{M}(order(A),irred))
    end

    itt = 0
    t = time()
    while !termination_criterion_der_red_map_precomp(stairs,nd,l)
        itt += 1
        newlm = SortedSet{M}(order(A))
        next_incr = GD_prereduction_increment!(echelon,next_incr,g1,SortedSet{M}(order(A)),A;newlm = newlm)
        union!(lm_ech, newlm)
        understair, irred = next_slice(understair, lm_g1, lm_ech, A)
        push!(stairs,SortedSet{M}(order(A),irred))
        update_stairs!(stairs,newlm,nd,A)
    end

    spol = reduce_with_echelon!(echelon,spol,A)
    # error("fin")
    return spol, g1, red_dt, echelon
end



function find_red_dt(g :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T,M}
    for (j,p) in enumerate(g) 
        if all(mon(p,1)[i] == 0 for i in 2:nvars(A)) && mon(p,1)[1] == 1
            tmp = p 
            deleteat!(g,j)
            return tmp 
        end
    end
    error("No operator with lm equals to dt found.")
    return zero(A)
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

# version with geobucket
function find_der_red_map(spol :: OrePoly, g1 :: Vector{OrePoly{T,M}}, red_dt :: OrePoly{T,M}, echelon :: Vector{OrePoly{T,M}},A :: OreAlg,geob :: GeoBucket,tmp_poly :: ReuseOrePoly) where {T, M}
    im = Dict{M,OrePoly{T,M}}()
    toadd = SortedSet{M}(order(A),mons(spol))
    cdt = one(ctx(A))
    mdt = makemon(1,A)
    c = opp(inv(coeff(red_dt,1),ctx(A)),ctx(A))

    dt = OrePoly([one(ctx(A))], [makemon(1,A)])


    while !isempty(toadd)
        m = pop!(toadd)
        dtm = OrePoly([cdt],[mdt*m])
        init_GeoBucket!(geob,dtm)
        addmul_geobucket!(geob,c,m,red_dt,A)
        res = GD_reduction1!(g1, A,geob,tmp_poly)
        res = reduce_with_echelon!(echelon, res,A,geob,tmp_poly)

        im[m] = res
        for m in mons(res)
            if !haskey(im,m)
                push!(toadd,m)
            end
        end 
    end
    return im 
end

# If tracer = true, a trace is returned
function der_red_map(A :: OreAlg, spol :: OrePoly, gb :: Vector{OrePoly{T,M}},param :: MCTParam;tracer :: Val{B} = Val(false)) where {T,M,B}
    if use_geobucket(param) 
        geob = GeoBucket(zero(A))
        tmp_poly = ReuseOrePoly(1,A)
        if B 
            spol, g1, red_dt, echelon, trace = der_red_map_precomp2(spol,gb,A,geob,tmp_poly,tracer = Val(B))
            tmp = find_der_red_map(spol,g1,red_dt,echelon,A,geob,tmp_poly)
            return (tmp, spol), trace
        else 
            spol, g1, red_dt, echelon = der_red_map_precomp2(spol,gb,A,geob,tmp_poly)
            return find_der_red_map(spol,g1,red_dt,echelon,A,geob,tmp_poly), spol
        end 
    else
        spol, g1, red_dt, echelon = der_red_map_precomp(spol,gb,A)
        return find_der_red_map(spol,g1,red_dt,echelon,A), spol
    end
end

# subsequent call where the trace is used to speed up computations 
# it assumes geobuckets are used 
function der_red_map(A :: OreAlg, trace :: TracerDerRedMap, spol :: OrePoly, gb :: Vector{OrePoly{T,M}},param :: MCTParam) where {T,M}
    geob = GeoBucket(zero(A))
    tmp_poly = ReuseOrePoly(1,A)
    red_dt = find_red_dt(gb,A)
    g1,g2 = separate(gb,A)
    echelon, _ = GD_prereduction_init_apply(g2, g1, trace.i,trace.s, A,geob,tmp_poly)

    spol = GD_reduction1!(spol,g1,A,geob,tmp_poly)
    spol = reduce_with_echelon!(echelon,spol,A,geob,tmp_poly) 

    debug(param) && @debug "starting find_der_red_map at time $(time()). "

    tmp = find_der_red_map(spol,g1,red_dt,echelon,A,geob,tmp_poly)  
    debug(param) && @debug "der_red_map computed at time $(time()). "
    debug(param) && @debug "the echelon form has length $(length(echelon))."
    return tmp, spol
end