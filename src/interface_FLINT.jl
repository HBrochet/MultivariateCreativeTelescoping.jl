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


# slong _nmod_poly_hgcd(mp_ptr *M, slong *lenM, mp_ptr A, slong *lenA, mp_ptr B, slong *lenB, mp_srcptr a, slong lena, mp_srcptr b, slong lenb, nmod_t mod)
# function half_gcd(a :: fpPolyRingElem, b :: fpPolyRingElem)
#     p = parent(a)
#     A = p()
#     B = p()

#     lenM = UInt(4)
#     M = Vector{fpPolyRingElem}(undef,4) 
#     for i in 1:4 
#         M[i] = p()
#     end

#     ccall((:_nmod_poly_hgcd, libflint), 
#     UInt,
#     (Ref{typeof(M)},Ptr{UInt},Ref{fpPolyRingElem},Ptr{UInt},Ref{fpPolyRingElem},Ptr{UInt},Ref{fpPolyRingElem},Ptr{UInt},Ref{fpPolyRingElem},Ptr{UInt},Ptr{UInt} ),
#     M  ,lenM,A,v,len)

