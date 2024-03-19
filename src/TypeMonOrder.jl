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

length_b1(b :: Block{N, O1, O2}) where {N,N1,N2,O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} = N1
length_b2(b :: Block{N, O1, O2}) where {N,N1,N2,O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} = N2


Block(o1,ind1, o2,ind2) = Block{nvars(o1),nvars(o2),nvars(o1)+nvars(o2), typeof(o1), typeof(o2)}(o1,ind1, o2,ind2)


Base.split(b::Block{N}, e::SVector{N}) where {N} =
    (SVector{length_b1(b)}(e[i] for i in b.ind1), SVector{length_b2(b)}(e[i] for i in b.ind2))

function ordervector(b::Block{N}, e::SVector{N}) where N 
    spl = split(b, e)
    return concat(ordervector(b.o1, spl[1]), ordervector(b.o2, spl[2]))
end



function make_order(s ::String,strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    blocs = split(s, ">")

    str = split(strip(blocs[1])," ")
    ord1, ind1 = make_order_aux(str,strvar_to_indexp,Val(M))
    if length(blocs) == 1 
        return ord1 
    end

    for i in 2:length(blocs) 
        str =  split(strip(blocs[i])," ")
        ord2, ind2 = make_order_aux(str,strvar_to_indexp,Val(M))
        ord1 = Block(ord1, ind1, ord2, ind2)
        ind1 = symdiff(ind1,ind2)
    end
    return ord1
end

function make_order_aux(str :: Vector{SubString{String}},strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    len = length(str)
    vec = Vector{Vector{E}}
    if str[1] == "grevlex" 
        ord = Grevlex(len-1)
    elseif str[1] == "lex"
        ord = Lex(len-1)
    else
        error("order not recognised")
    end
    return (ord, SVector{len-1}(strvar_to_indexp[str[i]] for i in 2:length(str)))
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

