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
    return convert(mod(x,Int(ctx.char)),ctx)
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


#. Univariate rational functions
abstract type RatFunCtx{T,TT} <: NmodLikeΓ{T, TT} end
abstract type UnivRatFunCtx{T,TT} <: RatFunCtx{T,TT} end

#. ratfun mod small p 

const UnivRatFunModp = Generic.FracFieldElem{fpPolyRingElem}

struct UnivRatFunModpCtx <: UnivRatFunCtx{UnivRatFunModp, UnivRatFunModp}
    F :: Generic.FracField{Nemo.fpPolyRingElem}
    R :: fpPolyRing
    vars :: Vector{fpPolyRingElem}
    char :: Int  
end

opp(a :: UnivRatFunModp, ctx :: UnivRatFunModpCtx) = - a
add(a :: UnivRatFunModp, b :: UnivRatFunModp, ctx ::UnivRatFunModpCtx) = a + b 
sub(a :: UnivRatFunModp, b :: UnivRatFunModp, ctx ::UnivRatFunModpCtx) = a - b 
mul(a :: UnivRatFunModp, b :: UnivRatFunModp, ctx ::UnivRatFunModpCtx) = a * b 
inv(a :: UnivRatFunModp, ctx :: UnivRatFunModpCtx) = 1/a 
submul(a :: UnivRatFunModp, b :: UnivRatFunModp, c :: UnivRatFunModp, ctx ::UnivRatFunModpCtx) = a - b * c 
normal(a :: UnivRatFunModp, ctx :: UnivRatFunModpCtx) = a 
inflate(a :: UnivRatFunModp, ctx :: UnivRatFunModpCtx) = a 
deflate(a :: UnivRatFunModp, ctx :: UnivRatFunModpCtx) = a 
convert(a :: Integer, ctx :: UnivRatFunModpCtx) = ctx.F(a)
convertn(a :: Integer, ctx :: UnivRatFunModpCtx) = ctx.F(a)

Base.one(ctx :: UnivRatFunModpCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: UnivRatFunModpCtx) where T = ctx.F(0)
Base.zero(ctx :: UnivRatFunModpCtx) = ctx.F(0)

Base.iszero(x :: UnivRatFunModp, ctx :: UnivRatFunModpCtx) = x == ctx.F(0)
Base.isone(x :: UnivRatFunModp, ctx :: UnivRatFunModpCtx) = x == ctx.F(1)

function evaluate(a :: UnivRatFunModp, p :: T) where T 
    return Nemo.evaluate(Nemo.numerator(a,false),p)//Nemo.evaluate(Nemo.denominator(a,false),p)
end

function evaluate(a :: UnivRatFunModp, p :: Vector{T}) where T 
    if length(p) > 1 
        error("the vector p has size > 1")
    end
    return Nemo.evaluate(Nemo.numerator(a,false),p[1])//Nemo.evaluate(Nemo.denominator(a,false),p[1])
end

# ratfun mod large p 

const UnivRatFunModP = Generic.FracFieldElem{FpPolyRingElem}

struct UnivRatFunModPCtx <: UnivRatFunCtx{UnivRatFunModP, UnivRatFunModP}
    F :: Generic.FracField{Nemo.FpPolyRingElem}
    R :: FpPolyRing
    vars :: Vector{FpPolyRingElem}
    char :: Int  
end

opp(a :: UnivRatFunModP, ctx :: UnivRatFunModPCtx) = - a
add(a :: UnivRatFunModP, b :: UnivRatFunModP, ctx ::UnivRatFunModPCtx) = a + b 
sub(a :: UnivRatFunModP, b :: UnivRatFunModP, ctx ::UnivRatFunModPCtx) = a - b 
mul(a :: UnivRatFunModP, b :: UnivRatFunModP, ctx ::UnivRatFunModPCtx) = a * b 
inv(a :: UnivRatFunModP, ctx :: UnivRatFunModPCtx) = 1/a 
submul(a :: UnivRatFunModP, b :: UnivRatFunModP, c :: UnivRatFunModP, ctx ::UnivRatFunModPCtx) = a - b * c 
normal(a :: UnivRatFunModP, ctx :: UnivRatFunModPCtx) = a 
inflate(a :: UnivRatFunModP, ctx :: UnivRatFunModPCtx) = a 
deflate(a :: UnivRatFunModP, ctx :: UnivRatFunModPCtx) = a 
convert(a :: Integer, ctx :: UnivRatFunModPCtx) = ctx.F(a)
convertn(a :: Integer, ctx :: UnivRatFunModPCtx) = ctx.F(a)

Base.one(ctx :: UnivRatFunModPCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: UnivRatFunModPCtx) where T = ctx.F(0)
Base.zero(ctx :: UnivRatFunModPCtx) = ctx.F(0)

Base.iszero(x :: UnivRatFunModP, ctx :: UnivRatFunModPCtx) = x == ctx.F(0)
Base.isone(x :: UnivRatFunModP, ctx :: UnivRatFunModPCtx) = x == ctx.F(1)

##
const UnivRatFunQQ = Generic.FracFieldElem{ZZPolyRingElem}

struct UnivRatFunQQCtx <: UnivRatFunCtx{UnivRatFunQQ,UnivRatFunQQ}
    F :: Generic.FracField{ZZPolyRingElem}
    R :: ZZPolyRing
    vars :: Vector{ZZPolyRingElem}
    char :: Int
    UnivRatFunQQCtx(F,R,v) = new(F,R,v,0)
end

opp(a :: UnivRatFunQQ, ctx :: UnivRatFunQQCtx) = - a
add(a :: UnivRatFunQQ, b :: UnivRatFunQQ, ctx ::UnivRatFunQQCtx) = a + b 
sub(a :: UnivRatFunQQ, b :: UnivRatFunQQ, ctx ::UnivRatFunQQCtx) = a - b 
mul(a :: UnivRatFunQQ, b :: UnivRatFunQQ, ctx ::UnivRatFunQQCtx) = a * b 
inv(a :: UnivRatFunQQ, ctx :: UnivRatFunQQCtx) = 1/a 
submul(a :: UnivRatFunQQ, b :: UnivRatFunQQ, c :: UnivRatFunQQ, ctx ::UnivRatFunQQCtx) = a - b * c 
normal(a :: UnivRatFunQQ, ctx :: UnivRatFunQQCtx) = a 
inflate(a :: UnivRatFunQQ, ctx :: UnivRatFunQQCtx) = a 
deflate(a :: UnivRatFunQQ, ctx :: UnivRatFunQQCtx) = a 
convert(a :: Int, ctx :: UnivRatFunQQCtx) = ctx.F(a)
convertn(a :: Integer, ctx :: UnivRatFunQQCtx) = ctx.F(a)

Base.one(ctx :: UnivRatFunQQCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: UnivRatFunQQCtx) where T = ctx.F(0)
Base.zero(ctx :: UnivRatFunQQCtx) = ctx.F(0)

Base.iszero(x :: UnivRatFunQQ, ctx :: UnivRatFunQQCtx) = x == ctx.F(0)
Base.isone(x :: UnivRatFunQQ, ctx :: UnivRatFunQQCtx) = x == ctx.F(1)



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

opp(a :: RatFunModp, ctx :: RatFunModpCtx) = - a
add(a :: RatFunModp, b :: RatFunModp, ctx ::RatFunModpCtx) = a + b 
sub(a :: RatFunModp, b :: RatFunModp, ctx ::RatFunModpCtx) = a - b 
mul(a :: RatFunModp, b :: RatFunModp, ctx ::RatFunModpCtx) = a * b 
inv(a :: RatFunModp, ctx :: RatFunModpCtx) = 1/a 
submul(a :: RatFunModp, b :: RatFunModp, c :: RatFunModp, ctx ::RatFunModpCtx) = a - b * c 
normal(a :: RatFunModp, ctx :: RatFunModpCtx) = a 
inflate(a :: RatFunModp, ctx :: RatFunModpCtx) = a 
deflate(a :: RatFunModp, ctx :: RatFunModpCtx) = a 
convert(a :: Integer, ctx :: RatFunModpCtx) = ctx.F(a)
convertn(a :: Integer, ctx :: RatFunModpCtx) = ctx.F(a)

Base.one(ctx :: RatFunModpCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: RatFunModpCtx) where T = ctx.F(0)
Base.zero(ctx :: RatFunModpCtx) = ctx.F(0)

Base.iszero(x :: RatFunModp, ctx :: RatFunModpCtx) = x == ctx.F(0)
Base.isone(x :: RatFunModp, ctx :: RatFunModpCtx) = x == ctx.F(1)

function evaluate(a :: RatFunModp, p :: Vector{T}) where T
    return Nemo.evaluate(Nemo.numerator(a,false),p)//Nemo.evaluate(Nemo.denominator(a,false),p)
end

# ratfun mod large p 

const RatFunModP = Generic.FracFieldElem{FpMPolyRingElem}

struct RatFunModPCtx <: MRatFunCtx{RatFunModP, RatFunModP}
    F :: Generic.FracField{Nemo.FpMPolyRingElem}
    R :: FpMPolyRing
    vars :: Vector{FpMPolyRingElem}
    char :: Int  
end

opp(a :: RatFunModP, ctx :: RatFunModPCtx) = - a
add(a :: RatFunModP, b :: RatFunModP, ctx ::RatFunModPCtx) = a + b 
sub(a :: RatFunModP, b :: RatFunModP, ctx ::RatFunModPCtx) = a - b 
mul(a :: RatFunModP, b :: RatFunModP, ctx ::RatFunModPCtx) = a * b 
inv(a :: RatFunModP, ctx :: RatFunModPCtx) = 1/a 
submul(a :: RatFunModP, b :: RatFunModP, c :: RatFunModP, ctx ::RatFunModPCtx) = a - b * c 
normal(a :: RatFunModP, ctx :: RatFunModPCtx) = a 
inflate(a :: RatFunModP, ctx :: RatFunModPCtx) = a 
deflate(a :: RatFunModP, ctx :: RatFunModPCtx) = a 
convert(a :: Integer, ctx :: RatFunModPCtx) = ctx.F(a)
convertn(a :: Integer, ctx :: RatFunModPCtx) = ctx.F(a)

Base.one(ctx :: RatFunModPCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: RatFunModPCtx) where T = ctx.F(0)
Base.zero(ctx :: RatFunModPCtx) = ctx.F(0)

Base.iszero(x :: RatFunModP, ctx :: RatFunModPCtx) = x == ctx.F(0)
Base.isone(x :: RatFunModP, ctx :: RatFunModPCtx) = x == ctx.F(1)

##
const RatFunQQ = Generic.FracFieldElem{ZZMPolyRingElem}

struct RatFunQQCtx <: MRatFunCtx{RatFunQQ,RatFunQQ}
    F :: Generic.FracField{ZZMPolyRingElem}
    R :: ZZMPolyRing
    vars :: Vector{ZZMPolyRingElem}
    char :: Int
    RatFunQQCtx(F,R,v) = new(F,R,v,0)
end

opp(a :: RatFunQQ, ctx :: RatFunQQCtx) = - a
add(a :: RatFunQQ, b :: RatFunQQ, ctx ::RatFunQQCtx) = a + b 
sub(a :: RatFunQQ, b :: RatFunQQ, ctx ::RatFunQQCtx) = a - b 
mul(a :: RatFunQQ, b :: RatFunQQ, ctx ::RatFunQQCtx) = a * b 
inv(a :: RatFunQQ, ctx :: RatFunQQCtx) = 1/a 
submul(a :: RatFunQQ, b :: RatFunQQ, c :: RatFunQQ, ctx ::RatFunQQCtx) = a - b * c 
normal(a :: RatFunQQ, ctx :: RatFunQQCtx) = a 
inflate(a :: RatFunQQ, ctx :: RatFunQQCtx) = a 
deflate(a :: RatFunQQ, ctx :: RatFunQQCtx) = a 
convert(a, ctx :: RatFunQQCtx) = ctx.F(a)
convertn(a, ctx :: RatFunQQCtx) = ctx.F(a)

Base.one(ctx :: RatFunQQCtx) = ctx.F(1)
Base.zero(x :: T, ctx :: RatFunQQCtx) where T = ctx.F(0)
Base.zero(ctx :: RatFunQQCtx) = ctx.F(0)

Base.iszero(x :: RatFunQQ, ctx :: RatFunQQCtx) = x == ctx.F(0)
Base.isone(x :: RatFunQQ, ctx :: RatFunQQCtx) = x == ctx.F(1)

function evaluate(a :: RatFunQQ, p :: Vector{T}) where T
    return Nemo.evaluate(Nemo.numerator(a,false),p)//Nemo.evaluate(Nemo.denominator(a,false),p)
end

## defining rational coefficient 
struct QQCtx <: NmodLikeΓ{QQFieldElem,QQFieldElem} 
    char :: Int 
    QQCtx() = new(0)
end

opp(a :: QQFieldElem, ctx :: QQCtx) = - a
add(a :: QQFieldElem, b :: QQFieldElem, ctx ::QQCtx) = a + b 
sub(a :: QQFieldElem, b :: QQFieldElem, ctx ::QQCtx) = a - b 
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


# struct SLPCtx{N,T,P} <: AbsContextCoeff{SLP{N,T}, SLP{N,T}}
#     R :: P
# end

# SLPCtx(N :: Int,R :: ZZMPolyRing) = SLPCtx{N,elem_type(R),typeof(R)}(R)


# add(a :: SLP, b :: SLP,:: SLPCtx) = SLP(:add,a,b)
# sub(a :: SLP, b :: SLP,:: SLPCtx) = SLP(:sub,a,b)
# mul(a :: SLP, b :: SLP,:: SLPCtx) = SLP(:mul,a,b)
# opp(a :: SLP,:: SLPCtx) = SLP(:opp,a)
# shift(a :: SLP, m::SVector,:: SLPCtx) = SLP(m,:shift,nothing,a,nothing)
# # mulshift(a :: SLP, b :: SLP,m::SVector,:: SLPCtx) = SLP(m,:mulshift,nothing,a,b)

# # creates a leaf
# SLP(a :: K, ctx :: SLPCtx{N,T,P}) where {K,N,T,P} =  SLP{N,T}(SVector{N,Int16}(0 for i in 1:N),:pol,ctx.R(a),nothing,nothing)

# convertn(a :: T,ctx::SLPCtx) where T = SLP(ctx.R(a),ctx)
# Base.one(ctx :: SLPCtx) = convertn(1,ctx)
# Base.zero(ctx :: SLPCtx) = convertn(0,ctx)

            


