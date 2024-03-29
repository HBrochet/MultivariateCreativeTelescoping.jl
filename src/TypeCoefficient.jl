# partially adapted from a file of Pierre Lairez

abstract type AbsCoefficient end

#= Expect:
    - Base.isone(::T) :: Bool
    - Base.iszero(::T) :: Bool
    - add(::NmodLikeΓ{T}, ::T, ::T) :: T
    - mul(::NmodLikeΓ{T}, ::T, ::T) :: T
    - inv(::NmodLikeΓ{T}, ::T) :: T
    - submul(::NmodLikeΓ{T, Tbuf}, ::Tbuf, ::T, ::T) :: Tbuf
    - normal(::NmodLikeΓ{T, Tbuf}, ::Tbuf) :: Tbuf
    - inflate(::NmodLikeΓ{T, Tbuf}, ::T) :: Tbuf
    - deflate(::NmodLikeΓ{T, Tbuf}, ::Tbuf) :: T

    - convert(::NomodLikeΓ{T, Tbuf}, ::Integer) :: T
    - convertn
    - opp 
    - normal =#


abstract type AbsContextCoeff{T, Tbuf} end

eltype1_ctx(ctx :: AbsContextCoeff{T, Tbuf}) where {T, Tbuf} = T
eltype2_ctx(ctx :: AbsContextCoeff{T, Tbuf}) where {T, Tbuf} = Tbuf

abstract type NmodLikeΓ{T, Tbuf} <: AbsContextCoeff{T, Tbuf} end

inflate(x :: T,::NmodLikeΓ{T, Tbuf}) where {T, Tbuf} = Tbuf(x)
deflate(x :: Tbuf,::NmodLikeΓ{T, Tbuf}) where {T, Tbuf} = x%T


#. Z/pZ

#.. Abstract type

abstract type NmodΓ{T<:Unsigned, Tbuf<:Unsigned}<: NmodLikeΓ{T,Tbuf} end

# (ctx::NmodΓ)(a::AA.ResElem) = ctx(a.data)
# abstractalgebra(ctx :: NmodΓ) = Nemo.NmodRing(UInt64(ctx.char))

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

convert(T :: Integer,ctx :: Nmod32Γ) = UInt32(T)
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

function opp(a::UInt32, ctx::Nmod32Γ)
    return ctx.char - a
end

function inv(a::UInt32, ctx::Nmod32Γ)
    invmod(UInt64(a), UInt64(ctx.char)) % UInt32
end

function mul(a, b, ctx::Nmod32Γ)
    (UInt64(a)*UInt64(b) % ctx.char) % UInt32
end

function normal(a::UInt64, ctx::Nmod32Γ)
    a % ctx.char
end

function convertn(x :: Integer,ctx :: Nmod32Γ)
    return convert(x,ctx) % ctx.char
end


#.. 32-bit modulus with 128-bits buffer
# currently not used

mutable struct Nmod32xΓ<:NmodΓ{UInt32, UInt128}
    char::UInt32
    twotothe64::UInt32        # 2^64 mod char

    function Nmod32xΓ(char)
        uchar = UInt32(char)
        p = (UInt64(2^32) % uchar)^2 % uchar % UInt32

        new(char, p)
    end
end


#... arithmetic operations
# Must implement inv, mul, submul and normal

function add(ctx::Nmod32xΓ, a::UInt32, b::UInt32)
    c0 = a+b
    c1 = c0 - ctx.char
    d = max(c0, c1)
    #@assert d%ctx.char == (a%ctx.char + b%ctx.char)%ctx.char
    return d
end

@inline function addmul(ctx::Nmod32xΓ, a::UInt128, b::UInt32, c::UInt32)
    z0 = a + UInt128(UInt64(b)*UInt64(c))
    return z0
end


@inline function submul(ctx::Nmod32xΓ, a::UInt128, b::UInt32, c::UInt32)
    z0 = a - UInt128(UInt64(b)*UInt64(c))
    return z0
end

function inv(ctx::Nmod32xΓ, a::UInt32)
    UInt32(invmod(UInt64(a), UInt64(ctx.char)))
end

function mul(ctx::Nmod32xΓ, a::UInt32, b::UInt32)
    UInt32((UInt64(a)*UInt64(b) % ctx.char) % UInt32)
end

function normal(ctx::Nmod32xΓ, a::UInt128)
    return a % ctx.char
    # hi = (a >> 64) % UInt64
    # lo = a % UInt64
    # add(ctx,
    #     ctx.twotothe64 * hi % ctx.char % UInt32,
    #     lo % ctx.char % UInt32)
end



#. Multivariate rational functions
abstract type RatFun{T,TT} <: NmodLikeΓ{T, TT} end

#. Multivariate rational functions

const RatFunModp = Generic.FracFieldElem{fpMPolyRingElem}

struct RatFunCtx <: RatFun{RatFunModp, RatFunModp}
    F :: Generic.FracField{Nemo.fpMPolyRingElem}
    R :: fpMPolyRing
    vars :: Vector{fpMPolyRingElem}
    char :: Int  
end

opp(a :: RatFunModp, ctx :: RatFunCtx) = - a
add(a :: RatFunModp, b :: RatFunModp, ctx ::RatFunCtx) = a + b 
mul(a :: RatFunModp, b :: RatFunModp, ctx ::RatFunCtx) = a * b 
inv(a :: RatFunModp, ctx :: RatFunCtx) = 1/a 
submul(a :: RatFunModp, b :: RatFunModp, c :: RatFunModp, ctx ::RatFunCtx) = a - b * c 
normal(a :: RatFunModp, ctx :: RatFunCtx) = a 
inflate(a :: RatFunModp, ctx :: RatFunCtx) = a 
deflate(a :: RatFunModp, ctx :: RatFunCtx) = a 
convert(a :: Integer, ctx :: RatFunCtx) = ctx.F(a)
convertn(a :: Integer, ctx :: RatFunCtx) = ctx.F(a)

Base.one(ctx :: RatFunCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: RatFunCtx) where T = ctx.F(0)
Base.zero(ctx :: RatFunCtx) = ctx.F(0)

Base.iszero(x :: RatFunModp, ctx :: RatFunCtx) = x == ctx.F(0)
Base.isone(x :: RatFunModp, ctx :: RatFunCtx) = x == ctx.F(1)



##
const RatFunQQ = Generic.FracFieldElem{ZZMPolyRingElem}

struct RatFunQQCtx <: RatFun{RatFunQQ,RatFunQQ}
    F :: Generic.FracField{ZZMPolyRingElem}
    R :: ZZMPolyRing
    vars :: Vector{ZZMPolyRingElem}
end

opp(a :: RatFunQQ, ctx :: RatFunQQCtx) = - a
add(a :: RatFunQQ, b :: RatFunQQ, ctx ::RatFunQQCtx) = a + b 
mul(a :: RatFunQQ, b :: RatFunQQ, ctx ::RatFunQQCtx) = a * b 
inv(a :: RatFunQQ, ctx :: RatFunQQCtx) = 1/a 
submul(a :: RatFunQQ, b :: RatFunQQ, c :: RatFunQQ, ctx ::RatFunQQCtx) = a - b * c 
normal(a :: RatFunQQ, ctx :: RatFunQQCtx) = a 
inflate(a :: RatFunQQ, ctx :: RatFunQQCtx) = a 
deflate(a :: RatFunQQ, ctx :: RatFunQQCtx) = a 
convert(a :: Int, ctx :: RatFunQQCtx) = ctx.F(a)
convertn(a :: Integer, ctx :: RatFunQQCtx) = ctx.F(a)

Base.one(ctx :: RatFunQQCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: RatFunQQCtx) where T = ctx.F(0)
Base.zero(ctx :: RatFunQQCtx) = ctx.F(0)

Base.iszero(x :: RatFunQQ, ctx :: RatFunQQCtx) = x == ctx.F(0)
Base.isone(x :: RatFunQQ, ctx :: RatFunQQCtx) = x == ctx.F(1)



## defining rational coefficient 
struct QQCtx <: NmodLikeΓ{QQFieldElem,QQFieldElem} end

opp(a :: QQFieldElem, ctx :: QQCtx) = - a
add(a :: QQFieldElem, b :: QQFieldElem, ctx ::QQCtx) = a + b 
mul(a :: QQFieldElem, b :: QQFieldElem, ctx ::QQCtx) = a * b 
inv(a :: QQFieldElem, ctx :: QQCtx) = 1/a 
submul(a :: QQFieldElem, b :: QQFieldElem, c :: QQFieldElem, ctx ::QQCtx) = a - b * c 
normal(a :: QQFieldElem, ctx :: QQCtx) = a 
inflate(a :: QQFieldElem, ctx :: QQCtx) = a 
deflate(a :: QQFieldElem, ctx :: QQCtx) = a 
convert(a :: Int, ctx :: QQCtx) = QQ(a)
convertn(a :: Integer, ctx :: QQCtx) = QQ(a)

Base.one(ctx :: QQCtx) = QQ(1)
Base.zero(x :: T, ctx :: QQCtx) where T = QQ(0)
Base.zero(ctx :: QQCtx) = QQ(0)

Base.iszero(x :: QQFieldElem, ctx :: QQCtx) = x == QQ(0)
Base.isone(x :: QQFieldElem, ctx :: QQCtx) = x == QQ(1)

