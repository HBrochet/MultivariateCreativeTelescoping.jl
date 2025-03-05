
function multivariate_rational_interpolation_known_support(F :: Function, supp :: OrePoly{K,M}, A :: OreAlg, args...) where {K,M} 
    nb_points = maximum(Nemo.length(Nemo.numerator(c)) + Nemo.length(Nemo.denominator(c)) - 1 for c in coeffs(supp))

    F_elem = elem_type(ctx(A).F)
    R = ctx(A).R
    C = coefficient_ring(R)
    C_elem = elem_type(C)
    n = number_of_variables(R)

    points = Vector{Vector{C_elem}}(undef,nb_points)
    evs = Vector{Vector{UInt32}}(undef,nb_points)
    for i in 1:nb_points 
        points[i] = [C(rand(Int)) for i in 1:n]

        nA = evaluate_coeff_algebra(points[i],A)
        ev_args = evaluate_coeff(points[i],nA,args...)

        evs[i] = coeffs(F(nA,ev_args...))
    end
    
    a = R(1)

    cs = Vector{F_elem}(undef,length(supp))
    for i in 1:length(supp)
        # println("i $i c $(coeff(supp,i))")
        l = length(Nemo.numerator(coeff(supp,i))) + length(Nemo.denominator(coeff(supp,i))) - 1 
        S = matrix_space(C, l, l)
        mat = S()
        b = Vector{C_elem}(undef,l)

        nu = Nemo.numerator(coeff(supp,i))
        ln = length(nu)

        de = Nemo.denominator(coeff(supp,i))
        ld = length(de)
        for j in 1:l 
            for k in 2:ln
                m = exponent_vector_fmpz(nu, k)
                set_exponent_vector!(a,1,m) 
                mat[j,k-1] = Nemo.evaluate(a,points[j])
            end
            c = C(evs[j][i])
            for k in 1:ld 
                m = exponent_vector_fmpz(de, k)
                set_exponent_vector!(a,1,m)
                mat[j,ln-1 + k] = Nemo.evaluate(a,points[j])*(-c) 
            end
            m = exponent_vector_fmpz(nu, 1)
            set_exponent_vector!(a,1,m)
            b[j] = -Nemo.evaluate(a,points[j])
        end
        vec = solve(mat,b,side = :right)
        # reconstruct r.f.
        m = exponent_vector_fmpz(nu, 1)
        set_exponent_vector!(a,1,m)
        num = deepcopy(a)
        for k in 2:ln 
            m = exponent_vector_fmpz(nu, k)
            set_exponent_vector!(a,1,m)
            num += vec[k-1]*a
        end
        den = zero(R)
        for k in 1:ld 
            m = exponent_vector_fmpz(de, k)
            set_exponent_vector!(a,1,m)
            den += vec[ln-1 + k]*a 
        end
        cs[i] = ctx(A).F(num) / ctx(A).F(den)
    end
    return OrePoly(cs,mons(supp))
end


function mri_crt(gens :: Vector{OrePoly{K,M}}, A :: OreAlg;tracer :: Val{B} = Val(false)) where {K,M,B} 
    param = f4_param(stophol = Val(true))
    tmp = multivariate_rational_interpolation(f4_mri,A,gens,param) 
    if B 
        return tmp, tmp
    else 
        return tmp 
    end
end

function mri_crt(trace :: OrePoly, gens :: Vector{OrePoly{K,M}}, A :: OreAlg) where {K,M}
    param = f4_param(stophol = Val(true))
    return multivariate_rational_interpolation_known_support(f4_mri,trace,A,gens,param)
end