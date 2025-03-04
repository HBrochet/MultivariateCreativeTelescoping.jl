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

#.. Weight order

struct Weight{N} <: AbsMonomialOrder{N}
    weights :: SVector{N, Int}
end
Weight(N ::Int, w:: SVector) = Weight{N}(w)

ordervector(o::Weight{N}, e::SVector{N}) where N = SVector{1}(sum(o.weights .* e))

""" Lexicographic on two blocks """
struct Block{N1, N2, N, O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} <: AbsMonomialOrder{N}
    o1 :: O1
    ind1 :: SVector{N1,Int}
    o2 :: O2
    ind2 :: SVector{N2,Int}
end
# Block(N1, N2, N, O1, O2) = Block{N1, N2, N, O1, O2}()
length_b1(b :: Block{N1, N2, N, O1, O2}) where {N,N1,N2,O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} = N1
length_b2(b :: Block{N1, N2, N, O1, O2}) where {N,N1,N2,O1 <: AbsMonomialOrder{N1}, O2 <: AbsMonomialOrder{N2}} = N2


Block(o1,ind1, o2,ind2) = Block{nvars(o1),nvars(o2),nvars(o1)+nvars(o2), typeof(o1), typeof(o2)}(o1,ind1, o2,ind2)
Block(o1,ind1, o2,ind2,ninds) = Block{nvars(o1),nvars(o2),ninds, typeof(o1), typeof(o2)}(o1,ind1, o2,ind2)



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



function make_order2(s ::String,strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    blocs = split(s, ">")
    filter!(b -> length(split(b)) > 1, blocs)
    ord, _ = make_order_rec(blocs,strvar_to_indexp,Val(M))
    return ord
end

function make_order_rec(blocs ::Vector{SubString{String}} ,strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    str = split(strip(blocs[1])," ")
    ord1, ind1 = make_order_aux(str,strvar_to_indexp,Val(M))

    if length(blocs) == 1 
        return ord1, ind1
    end
    ord2, ind2 = make_order_rec(blocs[2:end],strvar_to_indexp,Val(M))

    inds = Vector(ind1)
    append!(inds,ind2)
    inds = union(inds) #get every unique index in ind1 union ind2
    sort!(inds)

    ind1 = SVector{length(ind1)}(findfirst(x->x==i,inds) for i in ind1)
    ind2 = SVector{length(ind2)}(findfirst(x->x==i,inds) for i in ind2)

    return Block(ord1, ind1, ord2, ind2,length(inds)), SVector{length(inds)}(inds)
end

function make_order_aux(str :: Vector{SubString{String}},strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    let len
    if str[1] == "grevlex" 
        len = length(str)-1
        ord = Grevlex(len)
    elseif str[1] == "lex"
        len = length(str)-1
        ord = Lex(len)
    elseif str[1] == "weight"
        len = div(length(str)-1,2)
    else
        error("order not recognised")
    end
    ind= Int[]
    re = r"^-?\d+$"
    k=1
    w = Int[]
    for i in 2:length(str)
        if occursin(re,str[i])
            push!(w,parse(Int,str[i]))
            continue
        end
        l = strvar_to_indexp[str[i]]
        push!(ind,l)
        k = k+1 
    end
    if str[1] == "weight"
        ord = Weight(len,SVector{length(w),Int}(w))
    end
    return (ord, SVector{length(ind)}(ind))
end
end



# @generated function lt(t::AbsMonomialOrder{N}, e::SVector{N,E}, f::SVector{N,E}) where {N,E}
#     quote
#         oe = ordervector(t, e)
#         of = ordervector(t, f)

#         $([:(oe[$i] < of[$i] && return true ; oe[$i] > of[$i] && return false) for i in 1:N]...)

#         return false
#     end
# end

# function lt(t::AbsMonomialOrder, m1 :: OreMonVE, m2 :: OreMonVE)
#     return lt(t,m1.exp,m2.exp)
# end


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



##############################


function make_order(order ::String,dic :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    A = parse_order(order,dic,Val(M))
    return create_order(A,Val(M))
end

function parse_order(s ::String,strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    blocs = split(s, ">")
    vec = Vector{E}[]
    for l in blocs 
        str =  split(strip(l)," ")
        if length(str) == 1 
            continue 
        end
        append!(vec,ordervec(str,strvar_to_indexp,Val(M)))
    end
    return vec
end


function create_order(A :: Vector{Vector{I}},::Val{M}) where {M <: AbsOreMonomial,I <: Integer}
    N = nvars(M)
    namestruct = Symbol("Order",ord_ctr)
    B = adjoint(SMatrix{nvars(M),length(A)}( collect(Iterators.flatten(A))))

    eval(Meta.parse("struct $(namestruct) <: AbsMonomialOrder{$(N)} end"))
    prog = "function lt(order :: $(namestruct), a :: OreMonVE, b :: OreMonVE)\n"

    for i in 1:length(A)
        prog *= "t1 = Int16(0)"
        prog *= prod("+ a[$(j)]*Int16($(B[i,j]))" for j in 1:N)
        prog *= "\nt2 = Int16(0)"
        prog *= prod("+ b[$(j)]*Int16($(B[i,j]))" for j in 1:N)
        prog *= "\nif t1 < t2 
        return true
    elseif t1 > t2 
        return false
    end\n"
    end
    prog *= "return false\n end"

    prog2 = "function max_deg_block(order :: $(namestruct), a :: OreMonVE)\n"
    prog2 *= "v = SVector{$(N),Int16}("*prod("Int16($(B[1,i])), " for i in 1:N-1)*"Int16($(B[1,N])))\n"
    prog2 *= "return Int(sum(a.exp.*v))\nend"
    eval(Meta.parse(prog))
    eval(Meta.parse(prog2))

    global ord_ctr += 1
    return eval(Meta.parse("$(namestruct)()"))
end




function makeexp(::Val{M}, i :: Integer) where {N, E , M <: OreMonVE{N,E}}
    return E[j == i ? E(1) : E(0) for j in 1:N]
end

function ordervecgrevlex(str :: Vector{SubString{String}},strvar_to_indexp :: Dict{String,E},::Val{M}) where {E, M<:AbsOreMonomial}
    vec = Vector{E}[]
    push!(vec, sum([makeexp(Val(M), strvar_to_indexp[str[i]]) for i in 2:length(str)]))
    append!(vec,[-makeexp(Val(M), strvar_to_indexp[str[i]]) for i in length(str):-1:2])
    return vec
end

function orderveclex(str :: Vector{SubString{String}},strvar_to_indexp :: Dict{String,E},::Val{M}) where {E, M<:AbsOreMonomial}
    return [makeexp(Val(M), strvar_to_indexp[str[i]]) for i in 2:length(str)]
end






function ordervec(str :: Vector{SubString{String}},strvar_to_indexp :: Dict{String,E}, ::Val{M}) where {E, M<:AbsOreMonomial}
    len = length(str)
    vec = Vector{Vector{E}}
    if str[1] == "grevlex" 
        return ordervecgrevlex(str, strvar_to_indexp, Val(M))
    elseif str[1] == "lex"
        return orderveclex(str, strvar_to_indexp, Val(M))
    else
        error("order not recognised")
    end
end