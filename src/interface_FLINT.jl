# missing functions in Nemo.

function evaluate(pol :: fpPolyRingElem, v:: Vector{UInt})
    len = UInt(length(v))
    res = Vector{UInt}(undef,len)
    ccall((:nmod_poly_evaluate_nmod_vec_fast, libflint), 
        Nothing,
        (Ptr{UInt} ,Ref{fpPolyRingElem},Ptr{UInt},Int),
        res,pol,v,len)
    return res
end

function evaluate_many(r :: Generic.FracFieldElem{fpPolyRingElem}, v :: Vector{UInt};denisone :: Val{T} = Val(false)) where T
    if T
        num = Nemo.numerator(r,false)
        evnum = evaluate(num,v)
        return evnum 
    else
        den = Nemo.denominator(r,false)
        num = Nemo.numerator(r,false)
        evnum = evaluate(num,v)
        evden = evaluate(den,v)
        return evnum,evden
    end
end


function evaluate(pol :: fpMPolyRingElem, v:: Vector{fpFieldElem})
    z = @ccall libflint.nmod_mpoly_evaluate_all_ui(pol::Ref{fpMPolyRingElem}, v::Ptr{fpFieldElem}, parent(pol)::Ref{fpMPolyRing})::UInt
    return base_ring(parent(pol))(z)
end


function half_gcd_flint(a :: fpPolyRingElem, b :: fpPolyRingElem)
    R = parent(a)
    n = Nemo.degree(a)
    m = div(n,2) + isodd(n)
    if iszero(b) || Nemo.degree(b) < m
        return SMatrix{2,2,fpPolyRingElem}(R(1), R(0), R(0), R(1))
    end

    m11 = R()
    m12 = R()
    m21 = R()
    m22 = R()
    A = R()
    B = R()

    @ccall libflint.nmod_poly_hgcd(
        m11::Ref{fpPolyRingElem},
        m12::Ref{fpPolyRingElem},
        m21::Ref{fpPolyRingElem},
        m22::Ref{fpPolyRingElem},
        A::Ref{fpPolyRingElem},
        B::Ref{fpPolyRingElem},
        a::Ref{fpPolyRingElem},
        b::Ref{fpPolyRingElem}
    )::Int

    # FLINT returns M, A, B such that [a; b] = M * [A; B].
    # We return the inverse so half_gcd matches the Julia convention.
    det = m11*m22 - m12*m21
    if iszero(det)
        error("nmod_poly_hgcd returned a singular matrix")
    end
    invdet = inv(det)
    return SMatrix{2,2,fpPolyRingElem}(invdet*m22, invdet*(-m21), invdet*(-m12), invdet*m11)
end
