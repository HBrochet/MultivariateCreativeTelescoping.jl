# partially adapted from a file by Pierre Lairez

# File defining the coefficents types that can be used in this package.  
# The coefficients are stored in a type T, except for the output of addmul and submul operations 
# which are stored in a type Tbuf. A context is defined to store necessary informations in order to 
# perform operations on coefficients.  


### Abstract types 

abstract type AbsCoefficient end

#  Expect:
#     - Base.isone(::T) :: Bool
#     - Base.iszero(::T) :: Bool
#     - add(::T, ::T, ::NmodLikeΓ{T}) :: T
#     - mul(::T, ::T, ::NmodLikeΓ{T}) :: T
#     - inv(::T, ::NmodLikeΓ{T}) :: T
#     - opp(::T, ::NmodLikeΓ{T}) :: T
#     - submul(::Tbuf, ::T, ::T, ::NmodLikeΓ{T, Tbuf}) :: Tbuf
#     - normal(::Tbuf, ::NmodLikeΓ{T, Tbuf}) :: Tbuf
#     - Base.convert(::Integer, ::NomodLikeΓ{T, Tbuf}) :: T
#     - convertn(::Integer, ::NomodLikeΓ{T, Tbuf}) :: T
#     - inflate(::T, ::NmodLikeΓ{T, Tbuf}) :: Tbuf
#     - deflate(::Tbuf, ::NmodLikeΓ{T, Tbuf}) :: T


abstract type AbsContextCoeff{T, Tbuf} end

eltype1_ctx(ctx :: AbsContextCoeff{T, Tbuf}) where {T, Tbuf} = T
eltype2_ctx(ctx :: AbsContextCoeff{T, Tbuf}) where {T, Tbuf} = Tbuf

abstract type NmodLikeΓ{T, Tbuf} <: AbsContextCoeff{T, Tbuf} end

inflate(x :: T,::NmodLikeΓ{T, Tbuf}) where {T, Tbuf} = Tbuf(x)
deflate(x :: Tbuf,::NmodLikeΓ{T, Tbuf}) where {T, Tbuf} = x%T


#. Z/pZ

#.. Abstract type

abstract type NmodΓ{T<:Unsigned, Tbuf<:Unsigned}<: NmodLikeΓ{T,Tbuf} end

Base.zero(ctx :: NmodΓ{T,Tbuf}) where {T, Tbuf} = T(0) 
Base.one(ctx :: NmodΓ{T,Tbuf}) where {T, Tbuf} = T(1) 

Base.zero(T :: DataType,ctx :: NmodΓ) = T(0) 
Base.one(T :: DataType,ctx :: NmodΓ) = T(1) 

Base.iszero(x :: T, ctx :: NmodΓ{TT, Tbuf}) where {T,TT, Tbuf} = x == T(0)
Base.isone(x :: T, ctx :: NmodΓ{T, Tbuf}) where {T, Tbuf} = x == T(1)

#.. 32-bit modulus and 64-bit buffer 

struct Nmod32Γ<:NmodΓ{UInt32, UInt64}
    char::UInt32
    maxshift::UInt64

    function Nmod32Γ(char)
        uchar = UInt64(char)
        nbits = 64 - leading_zeros(uchar)
        @assert 2*nbits <= 64

        new(uchar, uchar << leading_zeros(uchar))
    end
end

Base.convert(T :: Integer,ctx :: Nmod32Γ) = UInt32(T)
#... arithmetic operations
# Must implement inv, mul, addmul and normal

@inline function add(a::UInt32, b::UInt32, ctx::Nmod32Γ)
    c0 = a+b
    if c0 >= ctx.char
        return c0 - ctx.char
    else 
        return c0
    end
end

# assume that we do not take the opposite of 0
function opp(a::UInt32, ctx::Nmod32Γ)
    return ctx.char - a
end

@inline function sub(a::UInt32, b::UInt32, ctx::Nmod32Γ)
    return add(a,opp(b,ctx),ctx)
end
@inline function addmul(a::UInt64, b::UInt32, c::UInt32,ctx::Nmod32Γ)
    z0 = a + UInt64(b)*UInt64(c)
    z1 = z0 - ctx.maxshift
    z = (z0 < a) ? z1 : z0
    return z
end

@inline function submul(a::UInt64, b::UInt32, c::UInt32, ctx::Nmod32Γ)
    z0 = a - UInt64(b)*UInt64(c)
    z1 = z0 + ctx.maxshift
    z = (z0 > a) ? z1 : z0
    #if (z%ctx.char + (UInt64(b)%ctx.char)*(UInt64(c)%ctx.char)%ctx.char)%ctx.char != a%ctx.char
    #    error((ctx, a, b, c))
    #end
    return z
end



function Base.inv(a::UInt32, ctx::Nmod32Γ)
    invmod(UInt64(a), UInt64(ctx.char)) % UInt32
end

function mul(a, b, ctx::Nmod32Γ)
    (UInt64(a)*UInt64(b) % ctx.char) % UInt32
end

function normal(a::UInt64, ctx::Nmod32Γ)
    a % ctx.char
end

function convertn(x :: Integer,ctx :: Nmod32Γ)
    return convert(mod(x,Int(ctx.char)),ctx)
end


### Rational functions

#. Univariate rational functions
abstract type RatFunCtx{T,TT} <: NmodLikeΓ{T, TT} end
abstract type UnivRatFunCtx{T,TT} <: RatFunCtx{T,TT} end

const RatFun = Union{
    Generic.FracFieldElem{fpPolyRingElem},
    Generic.FracFieldElem{FpPolyRingElem},
    Generic.FracFieldElem{ZZPolyRingElem},
    Generic.FracFieldElem{fpMPolyRingElem},
    Generic.FracFieldElem{FpMPolyRingElem},
    Generic.FracFieldElem{ZZMPolyRingElem},
}


# shared for all rational-function contexts
opp(a :: RatFun, ::RatFunCtx) = -a
add(a :: RatFun, b :: RatFun, ::RatFunCtx) = a + b
sub(a :: RatFun, b :: RatFun, ::RatFunCtx) = a - b
mul(a :: RatFun, b :: RatFun, ::RatFunCtx) = a * b
Base.inv(a :: RatFun, ::RatFunCtx) = inv(a)
submul(a :: RatFun, b :: RatFun, c :: RatFun, ::RatFunCtx) = a - b * c
normal(a :: RatFun, ::RatFunCtx) = a
inflate(a :: RatFun, ::RatFunCtx) = a
deflate(a :: RatFun, ::RatFunCtx) = a
Base.convert(a, ctx::RatFunCtx) = ctx.F(a)
convertn(a, ctx::RatFunCtx) = ctx.F(a)

Base.one(ctx::RatFunCtx) = ctx.F(1)
Base.zero(::Type{T}, ctx::RatFunCtx) where {T} = ctx.F(0)
Base.zero(ctx::RatFunCtx) = ctx.F(0)
Base.iszero(x :: RatFun, ctx::RatFunCtx) = x == ctx.F(0)
Base.isone(x :: RatFun, ctx::RatFunCtx) = x == ctx.F(1)

function evaluate(a::RatFun, p::Vector{S}) where {S}
    num = Nemo.numerator(a, false)
    den = Nemo.denominator(a, false)

    if a isa UnivRatFunModp || a isa UnivRatFunModP || a isa UnivRatFunQQ
        length(p) == 1 || error("the vector p has size $(length(p))")
        return Nemo.evaluate(a, p[1])
    end

    return Nemo.evaluate(num, p) // Nemo.evaluate(den, p)
end

function evaluate(a::RatFun, p)
    num = Nemo.numerator(a, false)
    den = Nemo.denominator(a, false)
    return Nemo.evaluate(num, p) // Nemo.evaluate(den, p)
end

#. ratfun mod small p 

const UnivRatFunModp = Generic.FracFieldElem{fpPolyRingElem}

struct UnivRatFunModpCtx <: UnivRatFunCtx{UnivRatFunModp, UnivRatFunModp}
    F :: Generic.FracField{Nemo.fpPolyRingElem}
    R :: fpPolyRing
    vars :: Vector{fpPolyRingElem}
    char :: Int  
end


# ratfun mod large p 

const UnivRatFunModP = Generic.FracFieldElem{FpPolyRingElem}

struct UnivRatFunModPCtx <: UnivRatFunCtx{UnivRatFunModP, UnivRatFunModP}
    F :: Generic.FracField{Nemo.FpPolyRingElem}
    R :: FpPolyRing
    vars :: Vector{FpPolyRingElem}
    char :: ZZRingElem 
end

# rafun over Q
const UnivRatFunQQ = Generic.FracFieldElem{ZZPolyRingElem}

struct UnivRatFunQQCtx <: UnivRatFunCtx{UnivRatFunQQ,UnivRatFunQQ}
    F :: Generic.FracField{ZZPolyRingElem}
    R :: ZZPolyRing
    vars :: Vector{ZZPolyRingElem}
    char :: Int
    UnivRatFunQQCtx(F,R,v) = new(F,R,v,0)
end

#. Multivariate rational functions

abstract type MRatFunCtx{T,TT} <: RatFunCtx{T,TT}  end

#. ratfun mod small p 

const RatFunModp = Generic.FracFieldElem{fpMPolyRingElem}

struct RatFunModpCtx <: MRatFunCtx{RatFunModp, RatFunModp}
    F :: Generic.FracField{Nemo.fpMPolyRingElem}
    R :: fpMPolyRing
    vars :: Vector{fpMPolyRingElem}
    char :: Int  
end


# ratfun mod large p 

const RatFunModP = Generic.FracFieldElem{FpMPolyRingElem}

struct RatFunModPCtx <: MRatFunCtx{RatFunModP, RatFunModP}
    F :: Generic.FracField{Nemo.FpMPolyRingElem}
    R :: FpMPolyRing
    vars :: Vector{FpMPolyRingElem}
    char :: ZZRingElem  
end

##
const RatFunQQ = Generic.FracFieldElem{ZZMPolyRingElem}

struct RatFunQQCtx <: MRatFunCtx{RatFunQQ,RatFunQQ}
    F :: Generic.FracField{ZZMPolyRingElem}
    R :: ZZMPolyRing
    vars :: Vector{ZZMPolyRingElem}
    char :: Int
    RatFunQQCtx(F,R,v) = new(F,R,v,0)
end

## defining rational coefficient 
struct QQCtx <: NmodLikeΓ{QQFieldElem,QQFieldElem} 
    char :: Int 
    QQCtx() = new(0)
end

opp(a :: QQFieldElem, :: QQCtx) = - a
add(a :: QQFieldElem, b :: QQFieldElem, ::QQCtx) = a + b 
sub(a :: QQFieldElem, b :: QQFieldElem, ::QQCtx) = a - b 
mul(a :: QQFieldElem, b :: QQFieldElem, ::QQCtx) = a * b 
Base.inv(a :: QQFieldElem, :: QQCtx) = 1/a 
submul(a :: QQFieldElem, b :: QQFieldElem, c :: QQFieldElem, ::QQCtx) = a - b * c 
normal(a :: QQFieldElem, :: QQCtx) = a 
inflate(a :: QQFieldElem, :: QQCtx) = a 
deflate(a :: QQFieldElem, :: QQCtx) = a 
Base.convert(a :: Int, :: QQCtx) = QQ(a)
convertn(a :: Integer, :: QQCtx) = QQ(a)

Base.one(:: QQCtx) = QQ(1)
Base.zero(:: T,:: QQCtx) where T = QQ(0)
Base.zero(:: QQCtx) = QQ(0)

Base.iszero(x :: QQFieldElem, :: QQCtx) = x == QQ(0)
Base.isone(x :: QQFieldElem, :: QQCtx) = x == QQ(1)
