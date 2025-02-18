
function MCT_internal(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg,param :: MCTParam) where {T,M}
    push!(A.nomul,1)
    par = CIParam{false,false}()
    der_map, spol = compute_with_cauchy_interpolation(der_red_map, A, spol, gb,param,param=par)
    debug(param) && @debug "der_red_map successfully reconstructed, starting to compute the LDE"
    deleteat!(A.nomul, length(A.nomul))
    return find_LDE_by_interpolation(der_map, spol, A)

end

"""
    MCT(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
"""
function MCT(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg;param :: MCTParam = mct_param()) where {T,M}
    if char(A) > 0 
        return MCT_internal(spol,gb,A,param)
    else
        return compute_with_CRT(MCT_internal,A, spol,gb,A,param)
    end
end





function find_LDE_by_interpolation(der_map :: Dict{M,OrePoly{T,M}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    rels, den = find_first_lin_dep_derivatives(der_map,spol,A)
    mat = mct_op_to_mat(rels, A)
    par = CIParam{true,false}()
    res = compute_with_cauchy_interpolation(my_kernel,A,mat,param=par)
    for i in 1:length(res)
        coeffs(res)[i] = coeffs(res)[i] * ctx(A).F(den)^(mons(res)[i][1]+1)
    end
    clear_denominators!(res,A)
    ngcd = gcd([Nemo.numerator(c,false) for c in coeffs(res)])
    mul!(1 // ctx(A).F(ngcd), res, A)
    return res
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

    den = denominator(spol,A)
    for k in keys(map_)
        den = lcm(den,denominator(map_[k],A))
    end
    Fden = ctx(A).F(den)
    for k in keys(map_)
        mul!(Fden,map_[k],A)
    end
    mul!(Fden,spol,A)
    nrel = spol
    rels = OrePoly{T,M}[]
    while true  
        push!(rels,nrel)
        nrel_red = evaluate_parameter(nrel,point,nA)
        nrel_red, v = reduce_with_echelon_augmented!(echelon_derivatives,nrel_red,nA,echelonvect)
        if length(nrel_red) == 0 
            return rels, Fden
        end
        add_echelon!(echelon_derivatives, nrel_red, nA, augmented = true, echelonvect = echelonvect, vect = v)
        nrel = compute_next_rel(map_, nrel,Fden,ord,A)
        ord += 1 
    end
end


function compute_next_rel(map_ :: Dict{M,OrePoly{T,M}}, pol :: OrePoly{T,M},den ::T,ord::Int, A :: OreAlg) where {T,M}
    cs = deepcopy(coeffs(pol))
    ms = deepcopy(mons(pol))
    for i in 1:length(cs)
        cs[i] = ctx(A).F(derivative(Nemo.numerator(cs[i],false))*den)
    end
    res = OrePoly(cs,ms)
    normalize!(res,A)
    if den != one(ctx(A))
        cpol = deepcopy(pol)
        mul!(ctx(A).F(ord+1)*derivative(den),cpol,A)
        res = sub!(res,cpol,A)
    end
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