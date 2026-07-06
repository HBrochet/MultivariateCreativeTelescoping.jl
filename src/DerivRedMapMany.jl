function find_red_dts(g :: Vector{OrePoly{T,M}}, A::OreAlg{T,C,M,O}) where {T,C,M,O}
    nrdv = A.nrdv
    red_dts = Vector{OrePoly{T,M}}(undef, nrdv)
    red_mons = [makemon(i,A) for i in 1:nrdv]
    found = falses(nrdv)
    delete_idxs = Int[]
    sizehint!(delete_idxs, nrdv)
    j = 1
    while j <= length(g)
        p = g[j]
        m = mon(p,1)
        idx = 0
        for i in 1:nrdv
            if m == red_mons[i]
                idx = i
                break
            end
        end
        if idx == 0
            j += 1
            continue
        end
        found[idx] && error("Several operators with lm equals to $(A.inp.ratdiffvars[2][idx]) found")
        red_dts[idx] = p
        found[idx] = true
        push!(delete_idxs, j)
        j += 1
    end
    isempty(delete_idxs) || deleteat!(g, delete_idxs)
    all(found) && return red_dts
    missing_ops = join(A.inp.ratdiffvars[2][findall(.!found)], ", ")
    error("No operator with lm equals to $(missing_ops) found")
end

function red_dt_tail(red_dt :: OrePoly{T,M}, A :: OreAlg{T,C,M,O}) where {T,C,M,O}
    length(red_dt) <= 1 && return zero(A)
    L = undefOrePoly(length(red_dt)-1,A)
    for i in 2:length(red_dt)
        L[i-1] = red_dt[i]
    end
    return L
end

function der_red_map_precomp_many2(spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}, A :: OreAlg{T,C,M,O}, geob :: GeoBucket, tmp_poly :: ReuseOrePoly, homogeneous :: Bool = false; rho::Integer = 1, homogeneous_degree::Integer = 1) where {T,C,M,O}
    rho_gap = homogeneous ? Int(rho) * (Int(homogeneous_degree) - 1) : Int(rho)
    red_dts = find_red_dts(gb,A)
    Ls = [red_dt_tail(red_dt, A) for red_dt in red_dts]
    red_dts_degx = maximum(red_dt_xgain(red_dt, A) for red_dt in red_dts)
    g1,g2 = separate(gb,A)
    if length(g2) == 0
        spol = GD_reduction1!(spol,g1,A,geob,tmp_poly)
        return spol, g1, red_dts, OrePoly{T,M}[]
    end

    spol = GD_reduction1!(spol,g1,A,geob,tmp_poly)
    ll = maximum(sum(mon(p,1)[i] for i in xvars_range(A)) for p in g2)
    s = max(max(red_dts_degx,0),ll)
    echelon, next_incr = GD_prereduction_init(g2, g1, s + rho_gap, A,geob,tmp_poly)
    spol = reduce_with_echelon!(echelon,spol,A,geob,tmp_poly)
    while true
        done = SortedSet{M}(order(A))
        todo = SortedSet{M}(order(A))
        append!(todo,mons(spol))
        restart = false
        while !isempty(todo)
            m = poplast!(todo)
            m_degx = sum(m[i] for i in xvars_range(A))
            if m_degx > s + rho_gap
                s += 1
                next_incr = GD_prereduction_increment!(echelon,next_incr,g1,SortedSet{M}(order(A)),A,geob,tmp_poly)
                spol = reduce_with_echelon!(echelon,spol,A,geob,tmp_poly)
                restart = true
                break
            end
            push!(done,m)
            for L in Ls
                isempty(L) && continue
                geob = addmul_geobucket!(geob,one(ctx(A)),m,L,A)
                Lm = GD_reduction1!(g1,A,geob,tmp_poly)
                Lm = reduce_with_echelon!(echelon,Lm,A,geob,tmp_poly)
                for mo in mons(Lm)
                    if !(mo in done)
                        push!(todo,mo)
                    end
                end
            end
        end
        restart && continue
        return spol, g1, red_dts, echelon
    end
end

function find_der_red_map_many(spol :: OrePoly{T,M}, g1 :: Vector{OrePoly{T,M}}, red_dts :: Vector{OrePoly{T,M}}, echelon :: Vector{OrePoly{T,M}}, A :: OreAlg{T,C,M,O}) where {T,C,M,O}
    der_maps = [Dict{M,OrePoly{T,M}}() for i in 1:length(red_dts)]
    basis_set = SortedSet{M}(order(A))
    toadd = SortedSet{M}(order(A),mons(spol))
    dts = [OrePoly([one(ctx(A))], [makemon(i,A)]) for i in 1:length(red_dts)]

    while !isempty(toadd)
        m = popfirst!(toadd)
        push!(basis_set,m)
        for i in 1:length(red_dts)
            dtm = mul(dts[i],OrePoly([one(ctx(A))],[m]),A)
            dtm = div!(dtm,[red_dts[i]],A)
            dtm = GD_reduction1!(dtm, g1, A)
            dtm = reduce_with_echelon!(echelon, dtm, A)
            der_maps[i][m] = dtm
            for mo in mons(dtm)
                if !(mo in basis_set)
                    push!(toadd,mo)
                end
            end
        end
    end
    return der_maps, collect(basis_set)
end

function find_der_red_map_many(spol :: OrePoly{T,M}, g1 :: Vector{OrePoly{T,M}}, red_dts :: Vector{OrePoly{T,M}}, echelon :: Vector{OrePoly{T,M}}, A :: OreAlg{T,C,M,O}, geob :: GeoBucket, tmp_poly :: ReuseOrePoly) where {T,C,M,O}
    der_maps = [Dict{M,OrePoly{T,M}}() for i in 1:length(red_dts)]
    basis_set = SortedSet{M}(order(A))
    toadd = SortedSet{M}(order(A),mons(spol))
    coeffs_red = [opp(inv(coeff(red_dts[i],1),ctx(A)),ctx(A)) for i in 1:length(red_dts)]

    while !isempty(toadd)
        m = popfirst!(toadd)
        push!(basis_set,m)
        for i in 1:length(red_dts)
            dtm = OrePoly([one(ctx(A))],[makemon(i,A)*m])
            init_GeoBucket!(geob,dtm)
            addmul_geobucket!(geob,coeffs_red[i],m,red_dts[i],A)
            res = GD_reduction1!(g1, A,geob,tmp_poly)
            res = reduce_with_echelon!(echelon, res,A,geob,tmp_poly)
            der_maps[i][m] = res
            for mo in mons(res)
                if !(mo in basis_set)
                    push!(toadd,mo)
                end
            end
        end
    end
    return der_maps, collect(basis_set)
end

function der_maps_to_matrices(der_maps :: Vector{Dict{M,OrePoly{T,M}}}, basis :: Vector{M}, A :: OreAlg{T,C,M,O}) where {T,C,M,O}
    S = matrix_space(parent(one(ctx(A))),length(basis),length(basis))
    bij_inv = Dict{M,Int}(m=>i for (i,m) in enumerate(basis))
    res = [S() for i in 1:length(der_maps)]
    for i in 1:length(der_maps)
        for (j,m) in enumerate(basis)
            for (c,mo) in der_maps[i][m]
                res[i][bij_inv[mo],j] = c
            end
        end
    end
    return res
end

"""
    der_red_map_many(A :: OreAlg, spol :: OrePoly, gb :: Vector{OrePoly}, param :: MCTParam = mct_param(), homogeneous::Bool = false; rho = 1, homogeneous_degree = 1)

For each parameter derivation in `A.inp.ratdiffvars[2]`, finds a relation with leading monomial
`d_ti`, computes a common stable monomial basis, then returns the corresponding derivation maps
and their matrices. The `j`th column of `der_mats[i]` is the image of `basis[j]`.
In the homogeneous case, `rho` is scaled to `rho * (homogeneous_degree - 1)`.
"""
function der_red_map_many(A :: OreAlg{T,C,M,O}, spol_ :: OrePoly{T,M}, gb_ :: Vector{OrePoly{T,M}}, param :: MCTParam = mct_param(), homogeneous ::Bool = false; rho::Integer = 1, homogeneous_degree::Integer = 1) where {T,C,M,O}
    gb = deepcopy(gb_)
    spol = deepcopy(spol_)
    geob = GeoBucket(zero(A))
    tmp_poly = ReuseOrePoly(1,A)
    spol, g1, red_dts, echelon = der_red_map_precomp_many2(spol,gb,A,geob,tmp_poly,homogeneous; rho = rho, homogeneous_degree = homogeneous_degree)
    der_maps, basis = find_der_red_map_many(spol,g1,red_dts,echelon,A,geob,tmp_poly)
    return (dops = copy(A.inp.ratdiffvars[2]),
            red_dts = red_dts,
            basis = basis,
            der_maps = der_maps,
            der_mats = der_maps_to_matrices(der_maps, basis, A),
            spol = spol)
end
