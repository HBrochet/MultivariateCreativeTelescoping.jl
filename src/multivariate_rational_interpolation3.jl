# this function assumes that the base field of A contains parameters and that on argument args, 
# the function F returns a vector of OrePoly. 
# It computes F(args...) using multivariate sparse evaluation/interpolation.
# The reference algorithms are Cuyt and Lee's algorithm as well as Zippel's algorithm 

function multivariate_rational_interpolation(F :: Function, A :: OreAlg, args...)
    R = ctx(A).R
    C = coefficient_ring(R)
    n = number_of_variables(R)

    R0, _ = polynomial_ring(C,"_z")

    shift = [C(rand(Int)) for i in 1:n]
    # shift = [C(0) for i in 1:n]
    a = [C(rand(Int)) for i in 1:n] # starting point 

    ### step 1: perform the change of variable F(w_1,...,w_n) -> F(zw_1 + s1,...wx_n + sn),
    ### where w_1 = point, w_i = a[i] and si = shift[i], and reconstruct F as a function of z.
    ### This gives the size of F(x_1,...,x_n), the total degree of each component of F 
    ### and the number of points necessary for each reconstruction, as well as bound, a parameter for *
    ### further Cauchy interpolation

    point = C(rand(Int))
    ev1, nb_points, bounds = mri_first_reconstruction(F,A,point, a,shift,R0,args...)
    
    # evsn contains the numerators of the output of F reconstructed in _z. 
    # more precisely, evsn[i][j][k][l]  is the reconstruction of the numerator the ith coordinate of F 
    # with its input evaluated at pts[j][k][l].
    # the index j ranges from 1 to n, the number of variables, and for each k, the points in  pts[j][k] 
    # are used to reconstruct one skeleton with variables 1,...j-1 in Zippel's algorithm

    pts, evsn, evsd = init_evs(ev1,point,R0,A)

    ### step 2: reconstruct every coefficient of F(zx_1 + s1,...zx_n + sn). 
    ### If not enough points have been computed, compute more

    b = false
    while !b
        b, i,d, cs = try_mpi(pts,evsn,evsd,a,shift,R,A)
        if b 
            return OrePoly(cs,mons(ev1)) 
        else
            add_ev_point!(i,d,pts,evsn,evsd,nb_points,bounds,a,shift,R0,F,A,args...)
        end
    end
end

function mri_first_reconstruction(F :: Function, A :: OreAlg, point :: T, a :: Vector{T}, shift :: Vector{T},R0 :: Ring, args...) where T 
    let res
    R = ctx(A).R
    C = coefficient_ring(R)
    n = number_of_variables(R)


    rps = [C(rand(Int))]
    po = [rps[1]*point + shift[1],Tuple(a[i]*rps[1] + shift[i] for i in 2:n)...]
    nA = evaluate_coeff_algebra(po,A)
    ev_args = evaluate_coeff(po,nA,args...)
    ev_res = [F(nA,ev_args...)]

    succeeded = false 
    ctr = 1 
    while true
        push!(rps,C(rand(Int)))
        ctr +=1
        let ctr = ctr
        po = [rps[ctr]*point + shift[1],Tuple(a[i]*rps[ctr] + shift[i] for i in 2:n)...]
        end
        nA = evaluate_coeff_algebra(po,A)
        ev_args = evaluate_coeff(po,nA,args...)
        push!(ev_res, F(nA,ev_args...))
        if succeeded 
            nres = cauchy_interpolation_mri(ev_res,rps,R0, A)
            if nres == res 
                nb_points = maximum(Nemo.degree(Nemo.numerator(c,false)) + Nemo.degree(Nemo.denominator(c,false)) + 1 for c in coeffs(res))
                bounds = [Nemo.degree(Nemo.numerator(c,false))+1 for c in coeffs(res)]
                globalstats.counters[:mri_nb_points] += ctr
                return res, nb_points, bounds 
            end 
            res = nres
        else
            try 
                res = cauchy_interpolation_mri(ev_res,rps,R0, A)
                succeeded = true
            catch 
                # @debug "failure, trying more points"
            end 
        end
        # if ctr == 10 
        #     println(rps)
        #     println(ev_res)
        #     res = cauchy_interpolation_mri(ev_res,rps,R0, A)
        #     error("fin")
        # end
    end
    end
end

function init_evs(ev1 :: OrePoly, point ::T,R0 :: Ring, A :: OreAlg) where T
    R = ctx(A).R 
    C = coefficient_ring(R)
    C_elem = elem_type(C)
    n = number_of_variables(R)
    R0_elem = elem_type(R0)

    pts = Vector{Vector{Vector{Vector{C_elem}}}}(undef,n)

    evsn = Vector{Vector{Vector{Vector{R0_elem}}}}(undef,length(ev1)) 
    evsd = Vector{Vector{Vector{Vector{R0_elem}}}}(undef,length(ev1)) 

    for i in 1:n 
        pts[i] = Vector{Vector{C_elem}}[]
    end

    for i in 1:length(ev1)
        evsn[i] = Vector{Vector{Vector{R0_elem}}}(undef,n) 
        evsd[i] = Vector{Vector{Vector{R0_elem}}}(undef,n) 
        for j in 1:n 
            evsn[i][j] = Vector{Vector{R0_elem}}[]
            evsd[i][j] = Vector{Vector{R0_elem}}[]
        end
    end

    # add ev1 to evsn and evsd
    push!(pts[1], [[point]])
    for i in 1:length(ev1)
        t1 = Nemo.numerator(coeff(ev1,i),false)
        t2 = Nemo.denominator(coeff(ev1,i),false)

        co = constant_coefficient(t2)
        t1 = t1/co
        t2 = t2/co

        push!(evsn[i][1], [t1]) 
        push!(evsd[i][1], [t2])
    end 
    return pts, evsn, evsd
end

function try_mpi(pts :: Vector{Vector{Vector{Vector{T}}}},evsn :: Vector{Vector{Vector{Vector{TT}}}},evsd :: Vector{Vector{Vector{Vector{TT}}}}, a :: Vector{T},shift :: Vector{T},R :: Ring, A :: OreAlg) where {T, TT}
    C = coefficient_ring(R)
    n = number_of_variables(R)
    x = gens(R)

    # todo: remember partial skeletons
    nskeleton = [[R(0) for j in 1:length(evsn[i][1][1][1])] for i in 1:length(evsn)]
    dskeleton = [[R(0) for j in 1:length(evsd[i][1][1][1])] for i in 1:length(evsd)]
    ivar_nskel = [[1 for j in 1:length(evsn[i][1][1][1])] for i in 1:length(evsn)] # next var to be reconstructed in nskeleton
    ivar_dskel = [[1 for j in 1:length(evsd[i][1][1][1])] for i in 1:length(evsd)]


    Rz, z = polynomial_ring(R,"_z")
    po = [x[i]*z + shift[i] for i in 1:n]
    for i in 1:length(evsn) 
        nlen = length(evsn[i][1][1][1])
        U = [R(0) for l in 1:nlen]
        for l in nlen:-1:1
            if iszero(Nemo.coeff(evsn[i][1][1][1], l-1))
                continue 
            end
            b, j, d, tmp = try_mpi_reconstruction(pts,evsn,i,l-1,a,shift,U,R,nskeleton, ivar_nskel)
            if !b 
                return false, j,d , [ctx(A).F(0)]
            end
            nskeleton[i][l] = tmp 
            tmp2 = Rz(Nemo.evaluate(tmp, po))
            for (j,c) in enumerate(coefficients(tmp2))
                U[j] += c 
            end 
            
        end
    end

    for i in 1:length(evsd) 
        dlen = length(evsd[i][1][1][1])
        V = [R(0) for l in 1:dlen]
        for l in dlen:-1:1
            if iszero(Nemo.coeff(evsd[i][1][1][1], l-1))
                continue 
            end
            b, j, d, tmp = try_mpi_reconstruction(pts,evsd,i,l-1,a,shift,V,R,dskeleton,ivar_dskel)
            if !b 
                return false, j,d, [ctx(A).F(0)]
            end
            dskeleton[i][l] = tmp 
            tmp2 = Rz(Nemo.evaluate(tmp, po))
            for (j,c) in enumerate(coefficients(tmp2))
                V[j] += c 
            end 
        end
    end
    cs = [ctx(A).F(0) for i in 1:length(evsn)]
    for i in 1:length(evsn)
        num = sum(nskeleton[i])
        den = sum(dskeleton[i])
        cs[i] = ctx(A).F(num) // ctx(A).F(den)
    end 
    return true, 0,0, cs
end


function try_mpi_reconstruction(pts :: Vector{Vector{Vector{Vector{T}}}},evs0 :: Vector{Vector{Vector{Vector{TT}}}},h :: Int, d :: Int, a :: Vector{T},shift :: Vector{T},U :: Vector{TTT},R :: Ring, skels :: Vector{Vector{TTT}},ivar :: Vector{Vector{Int}}) where {T, TT,TTT}
    C = coefficient_ring(R)
    C_elem = elem_type(C)

    n = number_of_variables(R)

    R0 = parent(evs0[1][1][1][1])

    skeleton = skels[h][d+1]
    low = ivar[h][d+1]
    for i in low:n # for each variable 
        len = length(pts[i])
        if (i < n && len < 1) || (i == 1 && len < 2)
            return false, i, 1, skeleton 
        end



        if i == 1 
            rs = Vector{C_elem}(undef,len-1)  # evaluation points of the first variable 
            evs = Vector{C_elem}(undef,len-1)
            evp = Vector{C_elem}(undef,n)
            for j in 2:n 
                evp[j] = a[j]
            end
            for j in 1:len-1 
                rs[j] = pts[i][j][1][1]
                evp[1] = rs[j]
                evs[j] = Nemo.coeff(evs0[h][i][j][1],d) - Nemo.evaluate(U[d+1], evp) 
            end
            tmp = interpolate(R0,rs,evs) 

            push!(rs, pts[i][len][1][1])
            push!(evs, Nemo.coeff(evs0[h][i][len][1],d) - Nemo.evaluate(U[d+1], [rs[len], (a[k] for k in 2:n)...]))
            
            tmp2 = interpolate(R0,rs,evs) 
            if tmp != tmp2
                return false, i, 1, univ_pol_to_multiv_pol(R0, R, tmp, 1)
            end
            skeleton = univ_pol_to_multiv_pol(R0, R, tmp, 1)
            skels[h][d+1] = deepcopy(skeleton)
            ivar[h][d+1] = 2
        elseif i <n #  1 < i n 
            lskel = length(skeleton)
            if length(pts[i][1]) < lskel 
                return false, i, lskel, skeleton 
            end
            npts = length(evs0[h][i])

            @assert all(length(pts[i][j]) >= lskel for j in 1:npts)
            rs = Vector{C_elem}(undef,npts)
            xs = Vector{Vector{C_elem}}(undef, npts)

            t = length(skeleton)
            S = matrix_space(C, t, t)
            mat = S() 
            b = Vector{C_elem}(undef,t)
            
            pt = Vector{C_elem}(undef,n)

            for j in 1:npts
                for (k,mon) in enumerate(monomials(skeleton)) 
                    point = [pts[i][j][k][l] for l in 1:i]
                    for l in 1:i 
                        pt[l] = pts[i][j][k][l]
                    end
                    for l in i+1:n 
                        pt[l] = a[l]
                    end
                    # fill mat and b 
                    b[k] = Nemo.coeff(evs0[h][i][j][k],d) - Nemo.evaluate(U[d+1], pt)
                    for (l,m) in enumerate(monomials(skeleton))
                        tmp =  Nemo.evaluate(m,pt)
                        mat[k,l] = tmp # mon should only have variables in x_1,...,x_i-1$
                    end
                end
                rs[j] = pts[i][j][1][i]
                xs[j] = solve(mat,b,side = :right)
            end

            ### construct the new skeleton 
            new_skeleton = R(0)
            for (j,mon) in enumerate(monomials(skeleton))
                evs = [xs[l][j] for l in 1:npts]
                tmp = interpolate(R0,rs,evs) 
                new_skeleton += univ_pol_to_multiv_pol(R0, R, tmp, i)*mon 
            end

            ### try again with one more point 
            # for (k,mon) in enumerate(monomials(skeleton)) 
            #     point = [pts[i][npts][k][l] for l in 1:i]
            #     pt = [point..., Tuple(a[l] for l in i+1:n)... ]

            #     # fill mat and b 
            #     b[k] = Nemo.coeff(evs0[h][i][npts][k],d) - Nemo.evaluate(U[d+1], pt)
            #     for (l,m) in enumerate(monomials(skeleton)) 
            #         tmp =  Nemo.evaluate(m,pt)
            #         mat[k,l] = tmp # mon should only have variables in x_1,...,x_i-1$
            #     end
            # end
            push!(rs, a[i])
            push!(xs, [c for c in coefficients(skels[h][d+1])])

            new_skeleton2 = R(0)
            for (j,mon) in enumerate(monomials(skeleton))
                evs = [xs[l][j] for l in 1:npts+1]
                tmp = interpolate(R0,rs,evs) 
                new_skeleton2 += univ_pol_to_multiv_pol(R0, R, tmp, i)*mon 
            end
            if new_skeleton != new_skeleton2 
                return false, i, lskel , skeleton 
            else 
                skeleton = new_skeleton 
                skels[h][d+1] = deepcopy(skeleton)
                ivar[h][d+1] = i+1
            end
        else # i = n
            # the last variable is easy to reconstruct as the polynomial to reconstruct is homogeneous of degree d 
            p = a[n] 
            new_skeleton = R(0)
            lvar = gens(R)[n] # last variable
            for t in terms(skeleton)
                c = Nemo.coeff(t,1)
                m = monomial(t,1)
                deg = total_degree(m)
                new_skeleton += c/p^(d-deg)*m*lvar^(d-deg)
            end
            skeleton = new_skeleton
            skels[h][d+1] = deepcopy(skeleton)
            ivar[h][d+1] = i+1
        end
        # println("current skeleton $(i)")
        # println(skeleton)
    end
    # println("final skeleton")
    # println(skeleton)
    return true, 0, 0, skeleton 
end

# add d evaluation points that will be used to reconstruct one evaluated skeleton for the ith variable
# if the vectors of evs[_][j][k] have size smaller that d, add points to them instead

function add_ev_point!(j :: Int,d :: Int ,pts, evsn ,evsd, nb_points :: Int,bounds :: Vector{Int},a :: Vector{T}, shift :: Vector{T}, R0 :: Ring,F :: Function ,A :: OreAlg, args...) where T
    R = ctx(A).R
    C = coefficient_ring(R)
    C_elem = elem_type(C)
    n = number_of_variables(R)
    M =eltype_mo(A)
    R0_elem = elem_type(R0)
    
    if length(evsn[1][j]) > 0 
        tmp = length(evsn[1][j][1])
        d = max(d, tmp)
    end
    let boo
    let d=d
        if length(evsn[1][j]) > 0 
            boo = any(length(evsn[1][j][k]) < d for k in length(evsn[1][j]))
        else
            boo = false 
        end
    end
    if length(evsn[1][j]) > 0 && boo
        # then add evaluations to the vectors of evs[_][j][k]  until they all have size d

        rps = Vector{C_elem}(undef,nb_points)
        ev_res = Vector{OrePoly{UInt32,M}}(undef,nb_points)
        po = Vector{C_elem}(undef,n)
        [C(rand(Int)) for i in 1:nb_points]

        for k in 1:length(evsn[1][j])
            for l in length(evsn[1][j][k])+1:d 
                point = Vector{C_elem}(undef,j)
                for i in 1:j-1 
                    point[i] = C(rand(Int))
                end
                point[j] = pts[j][k][1][j]
                for m in 1:nb_points 
                    rps[m] = C(rand(Int))
                    for i in 1:j 
                        po[i] = point[i]*rps[m] + shift[i]
                    end
                    for i in j+1:n 
                        po[i] = a[i]*rps[m] + shift[i]
                    end
                    nA = evaluate_coeff_algebra(po,A)
                    ev_args = evaluate_coeff(po,nA,args...)
                    ev_res[m] = F(nA,ev_args...)    
                    globalstats.counters[:mri_nb_points] += 1
                end
                ev = cauchy_interpolation_mri(ev_res,rps,R0, A,bounds=bounds)
                push!(pts[j][k],point)
                for i in 1:length(ev)
                    t1 = Nemo.numerator(coeff(ev,i),false)
                    t2 = Nemo.denominator(coeff(ev,i),false)
            
                    co = constant_coefficient(t2)
                    t1 = t1/co
                    t2 = t2/co
            
                    push!(evsn[i][j][k], t1) 
                    push!(evsd[i][j][k], t2)
                end 
            end
        end
    else
        # one more point is added to evs[_][j]
        p = C(rand(Int))
        push!(pts[j], Vector{C_elem}[])
        k = length(pts[j])
        for i in 1:length(evsn)
            push!(evsn[i][j], Vector{R0_elem}[])
            push!(evsd[i][j], Vector{R0_elem}[])
        end

        rps = Vector{C_elem}(undef,nb_points)
        ev_res = Vector{OrePoly{UInt32,M}}(undef,nb_points)
        po = Vector{C_elem}(undef,n)
        for l in 1:d 
            point = Vector{C_elem}(undef,j) 
            for i in 1:j-1 
                point[i] = C(rand(Int))
            end
            point[j] = p
            for m in 1:nb_points 
                rps[m] = C(rand(Int))
                for i in 1:j 
                    po[i] = point[i]*rps[m] + shift[i]
                end
                for i in j+1:n 
                    po[i] = a[i]*rps[m] + shift[i]
                end
                nA = evaluate_coeff_algebra(po,A)
                ev_args = evaluate_coeff(po,nA,args...)
                ev_res[m] = F(nA,ev_args...)
                globalstats.counters[:mri_nb_points] += 1
            end
            ev = cauchy_interpolation_mri(ev_res,rps,R0, A,bounds = bounds)

            push!(pts[j][k],point)
            for i in 1:length(ev)
                t1 = Nemo.numerator(coeff(ev,i),false)
                t2 = Nemo.denominator(coeff(ev,i),false)
        
                co = constant_coefficient(t2)
                t1 = t1/co
                t2 = t2/co
        
                push!(evsn[i][j][k], t1) 
                push!(evsd[i][j][k], t2)
            end 
        end
    end
    end
    return nothing
end

function univ_pol_to_multiv_pol(R1 :: PolyRing, R2 :: MPolyRing, p ::RingElement,i ::Int) 
    res = R2(0)
    cs = coefficients(p)
    xi = gen(R2,i)
    for j in 0:length(p)-1 
        res += cs[j]*xi^j 
    end
    return res 
end
