
struct CIParam{A,B} end 

function ci_param(denisone :: Val{A} = Val(false),
                  tracer :: Val{B} = Val(false)) where {A,B}
    return CIParam{A,B}()
end

denisone(p :: CIParam{A,B}) where {A,B} = A
tracer(p :: CIParam{A,B}) where {A,B} = B 


# It assumes that A has only one parameter
# the function f should be called as f(ev(A), ev(arg1...),arg2) where ev evaluates the parameter
function compute_with_cauchy_interpolation(f :: Function, A :: OreAlg, args...;param ::CIParam = ci_param())
    let ev_args;
    randpoints = Int[]
    npoints = 1
    bound = 1 
    succeeded = false
    globalstats.counters[:number_evaluation] += 1

    let prev_cbl

    glen = guess_length(args...)
    nA = evaluate_parameter_algebra(1,A) # first argument is not needed unless there are variables T in A
    vpoints, vargs = evaluate_parameter_many(glen,randpoints,nA,args...;denisone=Val(denisone(param)))
    let vargs = vargs
        ev_args =ntuple(length(args)) do i
            vargs[i][1]
    end; end 
    if tracer(param)
        tmp, trace = f(nA,ev_args...;tracer=Val(true))
        ev_res = [tmp]
    else 
        ev_res = [f(nA,ev_args...)]
    end
    vctr = 2 
    push!(randpoints,vpoints[1])
    # globalstats.counters[:number_evaluation] += 1
    t = ctx(A).vars[1]
    prd = t - vpoints[1] # product of the t-randpoints[i], precomputation for cauchy interpolation


    while true 
        if vctr > glen 
            vpoints,vargs = evaluate_parameter_many(glen,randpoints,nA,args...;denisone=Val(denisone(param)))
            vctr = 1
        end
        push!(randpoints,vpoints[vctr])
        let vargs = vargs; let vctr = vctr
            ev_args =ntuple(length(args)) do i
                vargs[i][vctr]
        end; end; end
        prd = prd*(t-vpoints[vctr])
        vctr +=1
        npoints += 1 

        globalstats.counters[:number_evaluation] += 1

        # globalstats.counters[:number_evaluation] += 1
        # @debug "evaluation of t at a point ($(npoints)th)" 
        if tracer(param)
            re = f(nA,trace,ev_args...)
        else
            re = f(nA,ev_args...)
        end
        push!(ev_res, re)

        # trying interpolation
        if succeeded 
            # @debug "interpolation to reconstruct stable mon set"
            rcbl = random_cbl(ev_res,nA)
            cbl = cauchy_interpolation(rcbl,randpoints, A;prd = prd)
            if Nemo.denominator(cbl,false) == Nemo.denominator(prev_cbl,false)
                ev_den = evaluate_parameter_cbl(cbl,randpoints,A)
                den = Nemo.denominator(cbl,false)
                return cauchy_interpolation_known_den(ev_res, randpoints,ev_den,den, A)
            end
            prev_cbl = cbl
            # @debug "reconstructions don't match, trying another point"
            succeeded = false
        elseif npoints == 2^bound + 1
            try 
                # @debug "interpolation to reconstruct stable mon set"
                rcbl = random_cbl(ev_res,nA)
                prev_cbl = cauchy_interpolation(rcbl,randpoints, A;prd = prd)
                succeeded = true
                # @debug "success, trying one more point"
            catch 
                # @debug "failure, trying more points"
            end 
            bound += 1
        end
        if npoints ==100
            # println("randpoints")
            # println(randpoints)
            # print("evres")
            # prettyprint(ev_res,nA)
            @debug "interpolation to reconstruct stable mon set"
            rcbl = random_cbl(ev_res,nA)
            prev_cbl = cauchy_interpolation(rcbl,randpoints, A;prd = prd)
            println("random cbl")
            println(prev_cbl)
            succeeded = true
            @debug "success, trying one more point"
            error("fin")
        end
    end
    end
end
end


function cauchy_interpolation(v :: Vector{UInt32}, points :: Vector{Int}, A :: OreAlg; prd ::fpPolyRingElem = Nothing)
    R = base_ring(ctx(A).R) 
    ev = [R(Int(c)) for c in v]
    evp = [R(p) for p in points]
    num, den = cauchy_interpolation(ctx(A).R, ctx(A).vars[1], evp, ev, div(length(points), 2),prd = prd)
    # println("degree num $(Nemo.degree(num)) and den $(Nemo.degree(den))")
    return ctx(A).F(num) / ctx(A).F(den)
end

# not used anymore since I use many = true
function cauchy_interpolation(pol :: Vector{OrePoly{K,M}}, points :: Vector{Int}, A :: OreAlg; prd ::fpPolyRingElem = Nothing) where {K,M}
    cs = Vector{eltype_co(A)}(undef,length(pol[1]))
    ms  = mons(pol[1])
    bound = div(length(points), 2)
    S = ctx(A).R
    R = base_ring(S)

    F = ctx(A).F
    var = ctx(A).vars[1]
    for i in 1:length(pol[1])
        ev_pol = [R(Int(coeff(pol[j],i))) for j in 1:length(pol)]
        ev_points = [R(p) for p in points]

        num, den = cauchy_interpolation(S,vary, ev_points, ev_pol, bound,prd = prd)

        num = sum([F(cc.data)*var^i for (i,c) in enumerate(coefficients(num)) for cc in coefficients(c)])
        den = sum([F(cc.data)*var^i for (i,c) in enumerate(coefficients(den)) for cc in coefficients(c)])
        cs[i] = num/den 
    end
    return OrePoly(cs,ms)
end

function cauchy_interpolation(vec :: Vector{Dict{M,OrePoly{T,M}}}, randpoints :: Vector{Int}, A :: OreAlg; prd ::fpPolyRingElem = Nothing) where {T,M}
    res = Dict{M,OrePoly{eltype_co(A),M}}()
    for m in keys(vec[1])
        tmp = [vec[i][m] for i in 1:length(randpoints)]
        res[m] = cauchy_interpolation(tmp, randpoints, A,prd = prd)
    end
    return res 
end

function cauchy_interpolation(vec :: Vector{Tuple{Dict{M,OrePoly{T,M}}, OrePoly{T,M}}}, randpoints :: Vector{Int}, A :: OreAlg; prd ::fpPolyRingElem = Nothing) where {T,M}
    return cauchy_interpolation([vec[i][1] for i in 1:length(vec)], randpoints, A,prd = prd), cauchy_interpolation([vec[i][2] for i in 1:length(vec)], randpoints, A,prd = prd)
end


# x :: points 
# y :: evaluated rational function at points
function cauchy_interpolation(S::PolyRing,var ::TT, x::Vector{T}, y::Vector{T},bound ::Int; prd ::fpPolyRingElem = Nothing) where {TT <: RingElement, T <: RingElement}
    if isnothing(prd)
        prevR = prod([var-xi for xi in x])
    else 
        prevR = prd 
    end   
     R = interpolate(S,x,y)

    #fast approach 
    R0 = prevR
    M = mat_gcd_remainder_with_degree_cond(prevR,R,bound)
    V = M[2,2]
    v = SVector{2}(prevR,R)
    R = (M*v)[2]

    if gcd(R0,V) == S(1)
        return R,V
    end
    error("Cauchy interpolation has no solution")
end

function cauchy_interpolation_naive(S::PolyRing,var ::TT, x::Vector{T}, y::Vector{T},bound ::Int; prd ::fpPolyRingElem = Nothing) where {TT <: RingElement, T <: RingElement}
    if isnothing(prd)
        prevR = prod([var-xi for xi in x])
    else 
        prevR = prd 
    end
    R = interpolate(S,x,y)
    U = S(0)
    V = S(1)
    prevU = S(1)
    prevV = S(0)
    R0 = copy(prevR)
    while Nemo.degree(R) >= bound 
        (Q,R),prevR = divrem(prevR,R),R
        prevU, U = U, prevU - Q*U
        prevV, V = V, prevV - Q*V 
    end
    if gcd(R0,V) == S(1)
        return R,V
    end
    error("Cauchy interpolation has no solution")
end
# see algorithmes efficaces en calcul formel
function half_gcd(A :: fpPolyRingElem, B :: fpPolyRingElem)
    R = parent(A)
    n = Nemo.degree(A)
    n2 = Nemo.degree(B)
    m = div(n,2) + isodd(n)
    if n2 < m 
        return SMatrix{2,2,fpPolyRingElem}(R(1),R(0),R(0),R(1))
    end

    f = R(collect(coefficients(A))[m+1:end])
    g = R(collect(coefficients(B))[m+1:end])
    M = half_gcd(f,g)

    Ap,Bp = M*SVector{2,fpPolyRingElem}(A,B)
    if Nemo.degree(Bp) < m 
        return M 
    end
    Q, Cp = divrem(Ap,Bp)
    Mp = SMatrix{2,2,fpPolyRingElem}(R(0),R(1),R(1),-Q)
    l = 2*m-Nemo.degree(Bp)
    b = R(collect(coefficients(Bp))[l+1:end])
    c = R(collect(coefficients(Cp))[l+1:end])
    Mpp = half_gcd(b,c)
    return Mpp*Mp*M
end

function mat_gcd_remainder_with_degree_cond(p :: fpPolyRingElem, q :: fpPolyRingElem, l :: Int) 
    n = Nemo.degree(p)
    v = SVector{2,fpPolyRingElem}(p,q)
    R = parent(p)

    if Nemo.degree(q) < l  
        return SMatrix{2,2,fpPolyRingElem}([R(1),R(0),R(0),R(1)])
    elseif l <= div(n,2)
        M = half_gcd(p,q)
        (pp,qq) = M*v
        if Nemo.degree(qq) < l
            return M 
        elseif p == pp
            # do one step here
           quo,r = divrem(pp,qq)
           pp, qq = qq, r
           M = SMatrix{2,2,fpPolyRingElem}([R(0),R(1),R(1),-quo])*M
        end
        Mrec= mat_gcd_remainder_with_degree_cond(pp,qq,l)
        return Mrec*M
    else # l> n/2
        m = 2*l-n 
        pp = R(collect(coefficients(p))[m+1:end])
        qq = R(collect(coefficients(q))[m+1:end])
        M = half_gcd(pp,qq)
        return M
    end
end









function cauchy_interpolation_known_den(vec :: Vector{Dict{M,OrePoly{T,M}}}, randpoints :: Vector{Int}, ev_den :: Vector{T},den ::TT, A :: OreAlg) where {T,TT,M}
    res = Dict{M,OrePoly{eltype_co(A),M}}()
    for m in keys(vec[1])
        tmp = [vec[i][m] for i in 1:length(randpoints)]
        res[m] = cauchy_interpolation_known_den(tmp, randpoints, ev_den,den, A)
    end
    return res 
end

function cauchy_interpolation_known_den(vec :: Vector{Tuple{Dict{M,OrePoly{T,M}}, OrePoly{T,M}}}, randpoints :: Vector{Int}, ev_den :: Vector{T},den ::TT, A :: OreAlg) where {T,TT,M}
    return cauchy_interpolation_known_den([vec[i][1] for i in 1:length(vec)], randpoints,ev_den,den, A), cauchy_interpolation_known_den([vec[i][2] for i in 1:length(vec)], randpoints,ev_den,den, A)
end

function cauchy_interpolation_known_den(vec :: Vector{Vector{OrePoly{T,M}}}, randpoints :: Vector{Int}, ev_den :: Vector{T}, den ::TT, A :: OreAlg) where {T,TT,M}
    return [cauchy_interpolation_known_den([vec[i][j] for i in 1:length(vec)], randpoints,ev_den,den, A) for j in 1:length(vec[1])]
end

function cauchy_interpolation_known_den(ev_res :: Vector{OrePoly{T,M}},randpoints ::Vector{Int}, ev_den :: Vector{T},den::TT, A ::OreAlg) where {T,TT,M}
    if length(ev_res[1]) == 0 
        return zero(A)
    end
    # multiply by den 
    nA = evaluate_parameter_algebra(randpoints[1],A)

    for i in 1:length(ev_res)
        mul!(ev_den[i],ev_res[i],nA)
    end
    
    R = ctx(A).R
    S = base_ring(R)
    Rrp = [S(r) for r in randpoints]
    co = Vector{eltype_co(A)}(undef,length(mons(ev_res[1])))
    rden = ctx(A).F(den)
    for i in 1:length(mons(ev_res[1]))
        num = interpolate(R,Rrp,[S(Int(coeff(ev_res[j],i))) for j in 1:length(randpoints)])
        co[i] = ctx(A).F(num) / rden 
    end    
    return OrePoly(co,mons(ev_res[1])) 
end



function random_cbl(ev_res::Vector{Tuple{Dict{M, OrePoly{T, M}}, OrePoly{T,M}}},A :: OreAlg) where {T,M}
    vres = [zero(ctx(A)) for i in 1:length(ev_res)]
    for k in keys(ev_res[1][1])
        for l in 1:length(ev_res[1][1][k])
            lambda = convertn(rand(Int),ctx(A))
            for i in 1:length(ev_res)
                vres[i] = add(vres[i],mul(coeff(ev_res[i][1][k],l),lambda,ctx(A)),ctx(A))
            end
        end
    end
    for l in 1:length(ev_res[1][2])
        lambda = convertn(rand(Int),ctx(A))
        for i in 1:length(ev_res)
            vres[i] = add(vres[i],mul(coeff(ev_res[i][2],l),lambda,ctx(A)),ctx(A))
        end
    end
    return vres 
end

function random_cbl(ev_res::Vector{OrePoly{T,M}},A :: OreAlg) where {T,M}
    vres = [zero(ctx(A)) for i in 1:length(ev_res)]
    for l in 1:length(ev_res[1])
        lambda = convertn(rand(Int),ctx(A))
        for i in 1:length(ev_res)
            vres[i] = add(vres[i],mul(coeff(ev_res[i],l),lambda,ctx(A)),ctx(A))
        end
    end
    return vres 
end

function random_cbl(ev_res::Vector{Vector{OrePoly{T,M}}},A :: OreAlg) where {T,M}
    vres = [zero(ctx(A)) for i in 1:length(ev_res)]
    for l in 1:length(ev_res[1])
        for ll in 1:length(ev_res[1][l])
            lambda = convertn(rand(Int),ctx(A))
            for i in 1:length(ev_res)
                vres[i] = add(vres[i],mul(coeff(ev_res[i][l],ll),lambda,ctx(A)),ctx(A))
            end
        end
    end
    return vres 
end




### Guess the number of points that we will need for Cauchy Interpolation by finding the largest degree+1 in the coefficients
function guess_length(pol :: OrePoly) 
    m = 0 
    for c in coeffs(pol)
        m = max(m,length(Nemo.numerator(c,false)))
        m = max(m,length(Nemo.denominator(c,false)))
    end
    return m
end
    
function guess_length(g :: Vector{OrePoly{T,M}})  where {T, M}
    m = 0 
    for p in g 
        m = max(m,guess_length(p))
    end
    return m
end

function guess_length(d :: Dict{M,OrePoly{T,M}})  where {T, M}
    m = 0 
    for k in keys(d)
        m = max(m,guess_length(d[k]))
    end
    return m
end
    
    
function guess_length(g :: Vector{Vector{OrePoly{T,M}}}) where  {T,M}
    m = 0 
    for p in g 
        m = max(m,guess_length(p))
    end
    return m
end

function guess_length(mat :: Generic.MatSpaceElem{Generic.FracFieldElem{Nemo.fpPolyRingElem}})
    m = 0 
    nc = number_of_columns(mat)
    nr = number_of_rows(mat)
    for j in 1:nc
        for i in 1:nr 
            m = max(m,length(Nemo.numerator(mat[i,j],false)))
            m = max(m,length(Nemo.denominator(mat[i,j],false)))
        end
    end
    return m
end
function guess_length(t :: Tuple)
    return guess_length(t[1])
end
function guess_length(:: Any)
    return 0
end

function guess_length(args...)
    m = 0 
    for arg in args 
        m = max(m,guess_length(arg))
    end
    return m
end

function add_rand_point!(vec :: Vector{Int}, p :: Int)
    while true 
        point = mod(rand(Int),p)
        if !(point in vec) 
            push!(vec,point)
            return point
        end
    end
end
