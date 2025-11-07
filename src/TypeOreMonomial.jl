abstract type AbsOreMonomial{E} end 


# represent a monomial by a vector of exponent
struct OreMonVE{N,E <: Integer} <: AbsOreMonomial{E}
    exp :: SVector{N,E}
end 
Base.show(io::IO, m :: OreMonVE) = println(io,"Ore monomial")


Base.exp(P :: OreMonVE) = P.exp
getexp(m :: OreMonVE, i :: Int)  = m.exp[i]
makemon(v ::Vector{E}) where E = OreMonVE(SVector{length(v)}(v))
degree(m ::OreMonVE) = sum(m.exp)

Base.length(:: OreMonVE{N,E}) where {N, E} = N

nvars(:: OreMonVE{N,E}) where {N, E} = N
nvars(::Type{OreMonVE{N,E}}) where {N,E} = N
exptype(::Type{OreMonVE{N,E}}) where {N,E} = E
Base.getindex(M :: OreMonVE, i :: Integer) = getindex(M.exp,i)

Base.:(*)(a :: OreMonVE, b :: OreMonVE) = OreMonVE(a.exp + b.exp) 
Base.:(/)(a :: OreMonVE, b :: OreMonVE) = OreMonVE(a.exp - b.exp) 
Base.:(^)(a :: OreMonVE{N,E}, i :: Integer) where {N,E} = OreMonVE(E(i) * a.exp)

@inline Base.lcm(a ::OreMonVE, b ::OreMonVE{N,E}) where {N,E} = OreMonVE(SVector{N,E}(max(a[i],b[i]) for i in 1:N))
Base.sum(a ::OreMonVE) = sum(a.exp)
