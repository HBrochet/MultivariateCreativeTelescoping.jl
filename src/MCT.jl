function MCT_internal(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
    push!(A.nomul,1)
    der_map, spol = compute_with_cauchy_interpolation(der_red_map, A, spol, gb; many=true)
    # der_map, spol = der_red_map(A, spol, gb)

    deleteat!(A.nomul, length(A.nomul))
    # return find_LDE(der_map, spol, A)
    return find_LDE_by_interpolation2(der_map, spol, A)
    # return find_LDE_by_interpolation(der_map, spol, A)

end

#It assumes that A has the good ordering
"""
    MCT(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
"""
function MCT(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
    if char(A) > 0 
        return MCT_internal(spol,gb,A)
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
        nrel_red, v = reduce_with_echelon!(echelon_derivatives,nrel_red,A,augmented = true,echelonvect = echelonvect)
        if length(nrel_red) == 0 
            mons = [makemon(1,A)^(i-1) for i in 1:length(v)]
            res =OrePoly(v,mons)
            normalize!(res,A)
            return res
        end
        if ord == 1
            error("fin")
        end
        add_echelon!(echelon_derivatives, nrel_red, A, augmented = true, echelonvect = echelonvect, vect = v)
        nrel = compute_next_rel(map_, nrel,A)
        ord += 1 

    end
end

function find_LDE_by_interpolation(der_map :: Dict{M,OrePoly{T,M}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    rels = find_first_lin_dep_derivatives(der_map,spol,A)
    # return find_linear_relation(A,rels,A)
    # return res = find_linear_relation(A,rels,A)

    res = compute_with_cauchy_interpolation(find_linear_relation,A,rels,A;many=true)

    dgcd = gcd([Nemo.denominator(c,false) for c in coeffs(res)])
    ngcd = gcd([Nemo.numerator(c,false) for c in coeffs(res)])
    println("degree of gcds $(Nemo.degree(ngcd)) $(Nemo.degree(dgcd))")
    mul!(ctx(A).F(dgcd)/ctx(A).F(ngcd), res, A)
    
    return res
end


function find_LDE_by_interpolation2(der_map :: Dict{M,OrePoly{T,M}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    rels = find_first_lin_dep_derivatives(der_map,spol,A)
    mat = mct_op_to_mat(rels, A)
    # return my_kernel(A,mat)
    return  compute_with_cauchy_interpolation(my_kernel,A,mat,many=true)
end

function mct_op_to_mat(rels :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T,M}
    s = SortedSet{eltype_mo(A)}(order(A),m  for p in rels for m in mons(p))
    len = length(s)
    bijection = [pop!(s) for i in 1:len]
    bij_inv = Dict{eltype_mo(A),Int}(m=>i for (i,m) in enumerate(bijection))

    S = matrix_space(parent(coeff(rels[1],1)),length(rels),len)
    mat = S() 
    for i in 1:length(rels)
        for (c,m) in rels[i]
            mat[i,bij_inv[m]] = c
        end
    end
    return mat 
end

function my_kernel(A::OreAlg,mat :: MatElem)
    return line_mat_to_orepoly(kernel(mat),A)
end



function find_first_lin_dep_derivatives(map_ :: Dict{M,OrePoly{T,M}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    point = mod(rand(Int),char(A))
    nA = evaluate_parameter_algebra(point,A)
    TT = eltype_co(nA)
    echelon_derivatives = OrePoly{TT,M}[]
    echelonvect = Vector{TT}[]
    ord = 0
    nrel = spol
    rels = OrePoly{T,M}[]
    while true  
        push!(rels,nrel)
        nrel_red = evaluate_parameter(nrel,point,nA)
        nrel_red, v = reduce_with_echelon!(echelon_derivatives,nrel_red,nA,augmented = true,echelonvect = echelonvect)
        if length(nrel_red) == 0 
            return rels
        end
        add_echelon!(echelon_derivatives, nrel_red, nA, augmented = true, echelonvect = echelonvect, vect = v)
        nrel = compute_next_rel(map_, nrel,A)
        ord += 1 
    end
end

function find_linear_relation(A :: OreAlg, rels :: Vector{OrePoly{T,M}},:: OreAlg) where {T, M}
    echelon_derivatives = OrePoly{T,M}[]
    echelonvect = Vector{T}[]
    for i in 1:length(rels)
        rels[i], v = reduce_with_echelon!(echelon_derivatives,rels[i],A,augmented = true,echelonvect = echelonvect, unitary = true)
        if length(rels[i]) == 0 
            mons = [makemon(1,A)^(i-1) for i in 1:length(v)]
            res =OrePoly(v,mons)
            normalize!(res,A)
            return res
        end
        add_echelon!(echelon_derivatives, rels[i], A, augmented = true, echelonvect = echelonvect, vect = v, unitary = true)
    end
end

function compute_next_rel(map_ :: Dict{M,OrePoly{T,M}}, pol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    cs = deepcopy(coeffs(pol))
    ms = deepcopy(mons(pol))

    for i in 1:length(cs)
        cs[i] = derivative(cs[i])
    end
    res = OrePoly(cs,ms)
    normalize!(res,A)

    for (c,m) in pol 
        res = add!(res,mul(c,map_[m],A) ,A)
    end
    return res 
end

function  line_mat_to_orepoly(line_mat ::Nemo.fpMatrix, A :: OreAlg)
    l = number_of_columns(line_mat)
    res = undefOrePoly(l,A)
    for i in l:-1:1
        res[l-i+1] = (convertn(line_mat[1,i].data,ctx(A)), makemon(1,A)^(i-1))
    end
    normalize!(res,A)
    return res
end

function  line_mat_to_orepoly(line_mat ::Generic.MatSpaceElem{Generic.FracFieldElem{Nemo.fpPolyRingElem}}, A :: OreAlg)
    l = number_of_columns(line_mat)
    res = undefOrePoly(l,A)
    for i in l:-1:1
        res[l-i+1] = (line_mat[1,i], makemon(1,A)^(i-1))
    end
    normalize!(res,A)
    return res
end