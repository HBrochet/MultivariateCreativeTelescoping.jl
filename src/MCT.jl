
function _mct_coeff_derivative(c, A::OreAlg)
    xname = A.inp.ratdiffvars[1][1]
    return coeff_derivative(c, A.drvars_to_int[xname], A)
end

function MCT_internal(spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}, A::OreAlg{T,C,M,O}, param :: MCTParam) where {T,C,M,O}
    push!(A.nomul,1)
    par = ci_param(comp = Val(:fast), same_den = Val(false))
    der_map, spol = compute_with_cauchy_interpolation(der_red_map, A, spol, gb,param,param=par)
    debug(param) && @debug "der_red_map successfully reconstructed, starting to compute the LDE"
    deleteat!(A.nomul, length(A.nomul))
    return find_LDE_by_interpolation(der_map, spol, A)

end

"""
    MCT_direct(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg; param :: MCTParam = mct_param())

Compute MCT directly over the current coefficient field, without CRT or Cauchy interpolation.
"""
function MCT_direct(spol :: OrePoly{T,M}, _gb :: Vector{OrePoly{T,M}}, A::OreAlg{T,C,M,O}; param :: MCTParam = mct_param()) where {T,C,M,O}
    push!(A.nomul,1)
    gb = deepcopy(_gb)
    der_map, spol = der_red_map(A, spol, gb, param)
    debug(param) && @debug "der_red_map computed directly, starting to compute the LDE"
    deleteat!(A.nomul, length(A.nomul))
    return find_LDE_direct(der_map, spol, A)
end

function MCT_direct(A::OreAlg{T,C,M,O}, spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}; param :: MCTParam = mct_param()) where {T,C,M,O}
    return MCT_direct(spol, gb, A; param = param)
end

"""
    MCT_direct_nofrac(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg; param :: MCTParam = mct_param())

Compute MCT directly over the current coefficient field, without fraction-free processing
for finding the LDE.
"""
function MCT_direct_nofrac(spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}, A::OreAlg{T,C,M,O}; param :: MCTParam = mct_param()) where {T,C,M,O}
    push!(A.nomul,1)
    der_map, spol = der_red_map(A, spol, gb, param)
    debug(param) && @debug "der_red_map computed directly (no fraction-free LDE), starting to compute the LDE"
    deleteat!(A.nomul, length(A.nomul))
    return find_LDE_direct_nofrac(der_map, spol, A)
end

function MCT_direct_nofrac(A::OreAlg{T,C,M,O}, spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}; param :: MCTParam = mct_param()) where {T,C,M,O}
    return MCT_direct_nofrac(spol, gb, A; param = param)
end

"""
    MCT(spol :: OrePoly, gb :: Vector{OrePoly{T,M}}, A::OreAlg) where {T,M}
"""
function MCT(spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}, A::OreAlg{T,C,M,O};param :: MCTParam = mct_param()) where {T,C,M,O}
    if char(A) > 0 
        return MCT_internal(spol,gb,A,param)
    else
        return compute_with_CRT(MCT_internal,A, spol,gb,A,param)
    end
end

function MCT(A::OreAlg{T,C,M,O}, spol :: OrePoly{T,M}, gb :: Vector{OrePoly{T,M}}; param :: MCTParam = mct_param()) where {T,C,M,O}
    return MCT(spol,gb,A;param = param)
end





function find_LDE_by_interpolation(der_map :: Dict{M,OrePoly{T,M}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    rels, den = find_first_lin_dep_derivatives(der_map,spol,A)
    mat = mct_op_to_mat(rels, A)
    par = ci_param(;denisone=Val(true), comp = Val(:fast),same_den=Val(false))
    res = compute_with_cauchy_interpolation(my_kernel,A,mat,param=par)
    for i in 1:length(res)
        coeffs(res)[i] = coeffs(res)[i] * ctx(A).F(den)^(mons(res)[i][1]+1)
    end 
    clear_denominators!(res,A)
    ngcd = gcd([Nemo.numerator(c,false) for c in coeffs(res)])
    mul!(1 // ctx(A).F(ngcd), res, A)
    return res
end

function find_LDE_direct(der_map :: Dict{M,OrePoly{T,M}}, spol :: OrePoly{T,M}, A :: OreAlg) where {T,M}
    iszero(spol, A) && return one(A)

    den = denominator(spol, A; normalize = Val(true))
    for v in values(der_map)
        den = lcm(den,denominator(v, A; normalize = Val(true)))
    end
    Fden = ctx(A).F(den)
    for v in values(der_map)
        mul!(Fden,v,A)
    end
    mul!(Fden,spol,A)

    nrels = 1 + length(der_map)
    rels = Vector{OrePoly{T,M}}(undef,nrels)
    rels[1] = spol
    nrel = spol
    ord = 0
    for i in 2:nrels
        nrel = compute_next_rel(der_map, nrel, den, ord, A)
        rels[i] = nrel
        ord += 1
    end

    mat = mct_op_to_mat(rels, A)
    res = my_kernel(A, mat)
    for i in 1:length(res)
        coeffs(res)[i] = coeffs(res)[i] * Fden^(mons(res)[i][1]+1)
    end

    clear_denominators!(res,A)
    ngcd = gcd([Nemo.numerator(c,false) for c in coeffs(res)])
    mul!(1 // ctx(A).F(ngcd), res, A)
    return res
end

function mct_op_to_mat(rels :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T,M}
    s = SortedSet{eltype_mo(A)}(order(A),m  for p in rels for m in mons(p))
    len = length(s)
    bijection = [popfirst!(s) for i in 1:len]
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
        nrel_red = evaluate_parameter(nrel,point,1,nA)
        nrel_red, v = reduce_with_echelon_augmented!(echelon_derivatives,nrel_red,nA,echelonvect)
        if length(nrel_red) == 0 
            return rels, Fden
        end
        add_echelon!(echelon_derivatives, nrel_red, nA, augmented = true, echelonvect = echelonvect, vect = v)
        nrel = compute_next_rel(map_, nrel, den, ord, A)
        ord += 1 
    end
end


function compute_next_rel(map_ :: Dict{M,OrePoly{T,M}}, pol :: OrePoly{T,M}, den :: S, ord::Int, A :: OreAlg) where {T,M,S}
    cs = deepcopy(coeffs(pol))
    ms = deepcopy(mons(pol))
    for i in 1:length(cs)
        cs[i] = ctx(A).F(_mct_coeff_derivative(Nemo.numerator(cs[i], true), A) * den)
    end
    res = OrePoly(cs,ms)
    normalize!(res,A)
    if !isone(den)
        cpol = deepcopy(pol)
        mul!(ctx(A).F(ord + 1) * ctx(A).F(_mct_coeff_derivative(den, A)), cpol, A)
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

function  line_mat_to_orepoly(line_mat ::Generic.MatSpaceElem{<:Generic.FracFieldElem}, A :: OreAlg)
    l = number_of_columns(line_mat)
    res = undefOrePoly(l,A)
    for i in l:-1:1
        res[l-i+1] = (line_mat[1,i], makemon(1,A)^(i-1))
    end
    normalize!(res,A)
    return res
end
