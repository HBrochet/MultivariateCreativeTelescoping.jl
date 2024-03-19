
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
Base.show(io::IO, ::MIME"text/plain", v :: Vector{OrePoly}) = print(io,"Vector of Ore polynomials")



mons(P :: OrePoly) = P.mons 
Base.@propagate_inbounds mon(P :: OrePoly, i :: Integer) = P.mons[i]
coeffs(P :: OrePoly) = P.coeffs
Base.@propagate_inbounds coeff(P :: OrePoly, i :: Integer) = P.coeffs[i]

Base.length(P :: OrePoly) = length(coeffs(P))
Base.size(P :: OrePoly) = size(coeffs(P))
Base.copy(P::OrePoly) = OrePoly(copy(coeffs(P)),copy(mons(P)))

Base.@propagate_inbounds Base.getindex(P :: OrePoly, i :: Integer) = (getindex(P.coeffs, i), getindex(P.mons, i))
 
function Base.resize!(P :: OrePoly, r) 
    resize!(P.coeffs, r) 
    resize!(P.mons,r)
end

Base.@propagate_inbounds function Base.setindex!(P::OrePoly{K, M}, t::Tuple{K, M}, i) where {M, K}
    setindex!(P.coeffs, t[1], i)
    setindex!(P.mons, t[2], i)
end

