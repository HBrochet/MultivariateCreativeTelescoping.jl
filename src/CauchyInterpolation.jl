
function add_rand_point!(vec :: Vector{Int}, p :: Int)
    while true 
        point = mod(rand(Int),p)
        if !(point in vec) 
            push!(vec,point)
            return point
        end
    end
end

# It assumes that A has only one parameter
function compute_with_cauchy_interpolation(f :: Function, A :: OreAlg, args...)
    randpoints = Int[]
    npoints = 1
    bound = 1 
    succeeded = false
    let prev_res 

    @debug "evaluation of t at a point ($(npoints)th)" 
    point = add_rand_point!(randpoints,ctx(A).char)
    nA = evaluate_parameter_algebra(point,A)
    ev_res = [f(nA,evaluate_parameter(point,nA, args...)...)]

    while true 
        npoints += 1 
        point = add_rand_point!(randpoints,ctx(A).char)
        
        @debug "evaluation of t at a point ($(npoints)th)" 
        nA = evaluate_parameter_algebra(point,A)
        push!(ev_res, f(nA,evaluate_parameter(point,nA, args...)...))

        # trying interpolation
        if succeeded 
            @debug "interpolation to reconstruct stable mon set"
            res = cauchy_interpolation(ev_res, randpoints, A)

            if prev_res == res
                return res 
            else 
                @debug "reconstructions don't match, trying another point"
                prev_res = res
            end
        elseif npoints == bound^2 + 1
            try 
                @debug "interpolation to reconstruct stable mon set"
                prev_res = cauchy_interpolation(ev_res, randpoints, A)
                succeeded = true
                @debug "success, trying one more point"
            catch 
                @debug "failure, trying more points"
                bound += 1
            end 
        end
        # if npoints == 2
        #     @debug "interpolation to reconstruct stable mon set"
        #     prev_res = cauchy_interpolation(ev_res, randpoints, A)
        #     succeeded = true
        #     @debug "success, trying one more point"
        # end
    end
    end
end



function cauchy_interpolation(pol :: Vector{OrePoly{K,M}}, points :: Vector{Int}, A :: OreAlg) where {K,M}
    cs = Vector{eltype_co(A)}(undef,length(pol[1]))
    ms  = mons(pol[1])
    bound = div(length(points), 2)
    R, _ = polynomial_ring(Native.GF(ctx(A).char), "_x")
    S, vary = polynomial_ring(R, "_y")

    F = ctx(A).F
    var = ctx(A).vars[1]
    for i in 1:length(pol[1])
        ev_pol = [R(Int(coeff(pol[j],i))) for j in 1:length(pol)]
        ev_points = [R(p) for p in points]

        num, den = cauchy_interpolation(S,vary, ev_points, ev_pol, bound)

        num = sum([F(cc.data)*var^i for (i,c) in enumerate(coefficients(num)) for cc in coefficients(c)])
        den = sum([F(cc.data)*var^i for (i,c) in enumerate(coefficients(den)) for cc in coefficients(c)])
        cs[i] = num/den 
    end
    return OrePoly(cs,ms)
end

function cauchy_interpolation(vec :: Vector{Dict{M,OrePoly{T,M}}}, randpoints :: Vector{Int}, A :: OreAlg) where {T,M}
    res = Dict{M,OrePoly{eltype_co(A),M}}()
    for m in keys(vec[1])
        tmp = [vec[i][m] for i in 1:length(randpoints)]
        res[m] = cauchy_interpolation(tmp, [randpoints...], A)
    end
    return res 
end

function cauchy_interpolation(vec :: Vector{Tuple{Dict{M,OrePoly{T,M}}, OrePoly{T,M}}}, randpoints :: Vector{Int}, A :: OreAlg) where {T,M}
    return cauchy_interpolation([vec[i][1] for i in 1:length(vec)], randpoints, A), cauchy_interpolation([vec[i][2] for i in 1:length(vec)], randpoints, A)
end


# x :: points 
# y :: evaluated rational function at points
function cauchy_interpolation(S::PolyRing,var ::TT, x::Vector{T}, y::Vector{T},bound ::Int) where {TT <: RingElement, T <: RingElement}
    R = interpolate(S,x,y)
    U = S(0)
    V = S(1)
    prevR = prod([var-xi for xi in x])
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
    # println("CI failed")
    error("Cauchy interpolation failed")
end


