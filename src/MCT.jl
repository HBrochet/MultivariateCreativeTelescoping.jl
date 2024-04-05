# function wc_and_der_red_map(g :: Vector{OrePoly{T,M}},spol :: OrePoly, init :: WeylClosureInit, A :: OreAlg) where {T,M}
#     gb = weyl_closure(g,A,init.wci)
#     return der_red_map(spol, gb, A)
# end

function MCT_internal(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
    push!(A.nomul,1)
    der_map, spol = compute_with_cauchy_interpolation(der_red_map, A, spol, gb)
    deleteat!(A.nomul, length(A.nomul))
    return find_LDE(der_map, spol, A)
end

#It assumes that A has the good ordering
"""
    MCT(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
"""
function MCT(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
    if char(A) > 0 
        error("characteristic must be zero")
    else
        return compute_with_CRT(MCT_internal,A, spol,gb,A)
    end
end




function find_LDE(map_ :: Dict{M,OrePoly{T,M}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    nrel = spol
    @debug "looking for the LDE"
    echelon_derivatives = OrePoly{T,M}[]
    echelonvect = Vector{T}[]
    ord = 0
    while true  
        @debug "dealing with $(ord)th derivative"
        nrel_red = deepcopy(nrel)
        v = reduce_with_echelon!(echelon_derivatives,nrel_red,A,augmented = true,echelonvect = echelonvect)
        if length(nrel_red) == 0 
            mons = [makemon(1,A)^(i-1) for i in 1:length(v)]
            res =OrePoly(v,mons)
            normalize!(res,A)
            return res
        end
        add_echelon!(echelon_derivatives, nrel_red, A, augmented = true, echelonvect = echelonvect, vect = v)
        nrel = compute_next_rel(map_, nrel,A)
        ord += 1 
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