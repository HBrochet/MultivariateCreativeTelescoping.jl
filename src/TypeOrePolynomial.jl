
# Defining it as a subtype of AbstractVector permit to easily sort them
abstract type AbsOrePolynomial{K, M} <:AbstractVector{Tuple{K, M}}end


# represent a polynomial by a vector of monomial and coefficient 
struct OrePoly{K, M <: AbsOreMonomial} <: AbsOrePolynomial{K, M} 
    coeffs :: Vector{K}
    mons :: Vector{M}
end
Base.show(io::IO, p :: OrePoly) = print(io,"Ore polynomial")
Base.show(io::IO, ::MIME"text/plain", p :: OrePoly) = print(io,"Ore polynomial")

Base.show(io::IO, v :: Vector{OrePoly}) = print(io,"Vector of Ore polynomials")
Base.show(io::IO, ::MIME"text/plain", v :: Vector{OrePoly{T,M}}) where {T,M} = print(io,"Vector of Ore polynomials")



mons(P :: OrePoly) = P.mons 
Base.@propagate_inbounds mon(P :: OrePoly, i :: Integer) = P.mons[i]
Base.@propagate_inbounds lm(P :: OrePoly) = P.mons[1]
Base.@propagate_inbounds lc(P :: OrePoly) = P.coeffs[1]

coeffs(P :: OrePoly) = P.coeffs
Base.@propagate_inbounds coeff(P :: OrePoly, i :: Integer) = P.coeffs[i]

"""
    Base.length(p :: OrePoly)

Return the number of terms in the polynomial p.
"""
Base.length(P :: OrePoly) = length(coeffs(P))
Base.size(P :: OrePoly) = size(coeffs(P))
Base.copy(P::OrePoly) = OrePoly(copy(coeffs(P)),copy(mons(P)))
Base.iszero(p :: OrePoly) = length(p) == 0 

Base.@propagate_inbounds Base.getindex(P :: OrePoly, i :: Integer) = (getindex(P.coeffs, i), getindex(P.mons, i))
# Base.popfirst!(P :: OrePoly) = (popfirst!(P.coeffs),pôpfirst!(P.mons))
Base.zero(:: OrePoly{K,M})  where {K,M} = OrePoly(K[],M[])

function Base.resize!(P :: OrePoly, r) 
    resize!(P.coeffs, r) 
    resize!(P.mons,r)
end

Base.@propagate_inbounds function Base.setindex!(P::OrePoly{K, M}, t::Tuple{K, M}, i) where {M, K}
    setindex!(P.coeffs, t[1], i)
    setindex!(P.mons, t[2], i)
end

function degree(p :: OrePoly)
    return maximum([degree(m[2]) for m in p])
end


### Reusable OrePoly 

mutable struct ReuseOrePoly{K,M}
    op :: OrePoly{K,M}
    ind :: Int 
end


Base.length(p :: ReuseOrePoly) = length(p.op)
Base.iszero(p :: ReuseOrePoly) = p.ind == 0 

Base.@propagate_inbounds function Base.getindex(P :: ReuseOrePoly, i :: Integer)
    op = P.op 
    return     (getindex(op.coeffs, i), getindex(op.mons, i))
end
    
    

Base.@propagate_inbounds function Base.setindex!(P::ReuseOrePoly{K, M}, t::Tuple{K, M}, i) where {M, K}
    op = P.op 
    (c,m) = t
    setindex!(op.coeffs, c, i)
    setindex!(op.mons, m, i)
end

function grow_to!(p ::ReuseOrePoly, l :: Int)
    op = p.op
    resize!(op.coeffs,l)
    resize!(op.mons,l)
end

function Base.push!(p :: ReuseOrePoly,c :: K, m :: M) where {K,M}
    ind = p.ind
    if length(p) == ind 
        grow_to!(p,2*ind)
    end
    p[ind+1] = (c,m) 
    p.ind += 1
    nothing
end 

