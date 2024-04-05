abstract type AbsMonomialOrder{N} <: Base.Order.Ordering end

nvars(o :: AbsMonomialOrder{N}) where N = N
struct EmptyMORD{N} <: AbsMonomialOrder{N} end

function  Base.Order.lt(ord :: AbsMonomialOrder, a::Tuple{K, M}, b :: Tuple{K, M}) where {K, M <: AbsOreMonomial}
    return lt(ord, a[2], b[2])
end

function  Base.Order.lt(ord :: AbsMonomialOrder, a :: M, b :: M) where M <: AbsOreMonomial
    return lt(ord, a, b)
end








#.. Lexicographic order
struct Lex{N} <: AbsMonomialOrder{N} end
Lex(N) = Lex{N}()

ordervector(::Lex{N}, e::SVector{N}) where N = e

#.. Graded reverse lexicographic order
struct Grevlex{N} <: AbsMonomialOrder{N} end
Grevlex(N) = Grevlex{N}()

# should be compiled statically
ordervector(::Grevlex{N}, e::SVector{N}) where N = insert(deleteat(reverse(-e), N), 1, sum(e))


""" Lexicographic on two blocks """
struct Block{N1, N2, N, O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} <: AbsMonomialOrder{N}
    o1 :: O1
    ind1 :: SVector{N1,Int}
    o2 :: O2
    ind2 :: SVector{N2,Int}
end
Block(N1, N2, N, O1, O2) = Block{N1, N2, N, O1, O2}()
length_b1(b :: Block{N1, N2, N, O1, O2}) where {N,N1,N2,O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} = N1
length_b2(b :: Block{N1, N2, N, O1, O2}) where {N,N1,N2,O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} = N2


Block(o1,ind1, o2,ind2) = Block{nvars(o1),nvars(o2),nvars(o1)+nvars(o2), typeof(o1), typeof(o2)}(o1,ind1, o2,ind2)


Base.split(b::Block, e::SVector)  =
    (SVector{length_b1(b)}(e[i] for i in b.ind1), SVector{length_b2(b)}(e[i] for i in b.ind2))


@generated function _concat(::Size{p}, va :: StaticVector, ::Size{q}, vb :: StaticVector) where {p,q}
    :(SVector{p[1]+q[1]}($([:(va[$i]) for i in 1:p[1]]...), $([:(vb[$i]) for i in 1:q[1]]...)))
end
concat(va, vb) = _concat(Size(va), va, Size(vb), vb)

function ordervector(b::Block, e::SVector)
    spl = split(b, e)
    return concat(ordervector(b.o1, spl[1]), ordervector(b.o2, spl[2]))
end



function make_order(s ::String,strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    blocs = split(s, ">")
    filter!(b -> length(split(b)) > 1, blocs)
    ai = Int[] # list of already added exponent index
    ord, _ = make_order_rec(blocs,strvar_to_indexp,ai,Val(M))
    return ord
end

function make_order_rec(blocs ::Vector{SubString{String}} ,strvar_to_indexp :: Dict{String,E},ai :: Vector{Int}, ::Val{M}) where {E, M<:AbsOreMonomial}
    str = split(strip(blocs[1])," ")
    ord1, ind1 = make_order_aux(str,strvar_to_indexp,ai,Val(M))
    if length(blocs) == 1 
        return ord1, ind1
    end
    ord2, _ = make_order_rec(blocs[2:end],strvar_to_indexp,ai,Val(M))

    ind2 = setdiff(Set(i for i in 1:nvars(ord1) + nvars(ord2)),Set(ind1))
    ind2 = [i for i in ind2]
    sort!(ind2)
    ind2 = SVector{length(ind2)}(ind2)

    return Block(ord1, ind1, ord2, ind2), concat(ind1,ind2)
end

function make_order_aux(str :: Vector{SubString{String}},strvar_to_indexp :: Dict{String,E},ai :: Vector{Int}, ::Val{M}) where {E, M<:AbsOreMonomial}
    len = length(str)
    if str[1] == "grevlex" 
        ord = Grevlex(len-1)
    elseif str[1] == "lex"
        ord = Lex(len-1)
    else
        error("order not recognised")
    end
    v= Vector(undef,length(str)-1)
    for i in 1:length(str)-1
        l = strvar_to_indexp[str[i+1]] 
        v[i] = l - count(x -> x<l,ai)

    end
    append!(ai,[strvar_to_indexp[str[i]] for i in 2:length(str)])
    return (ord, SVector{len-1}(v))
end



@generated function lt(t::AbsMonomialOrder{N}, e::SVector{N}, f::SVector{N}) where N
    quote
        oe = ordervector(t, e)
        of = ordervector(t, f)

        $([:(oe[$i] < of[$i] && return true ; oe[$i] > of[$i] && return false) for i in 1:N]...)

        return false
    end
end

function lt(t::AbsMonomialOrder, m1 :: OreMonVE, m2 :: OreMonVE)
    return lt(t,m1.exp,m2.exp)
end


function Base.Order.lt(ord :: AbsMonomialOrder, a :: AbsOrePolynomial, b :: AbsOrePolynomial)
    for i in 1:min(length(a),length(b))
        if lt(ord,mon(a,i),mon(b,i))
            return true 
        elseif lt(ord,mon(a,i),mon(b,i))
            return false
        end
    end
    if length(a) > length(b) 
        return false
    elseif length(a) < length(b)
        return true 
    else
        return false
    end
end

