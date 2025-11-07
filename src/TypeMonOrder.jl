abstract type AbsMonomialOrder{N} <: Base.Order.Ordering end

nvars(:: AbsMonomialOrder{N}) where N = N
struct EmptyMORD{N} <: AbsMonomialOrder{N} end

function  Base.Order.lt(ord :: AbsMonomialOrder, a::Tuple{K, M}, b :: Tuple{K, M}) where {K, M <: AbsOreMonomial}
    return lt(ord, a[2], b[2])
end

function  Base.Order.lt(ord :: AbsMonomialOrder, a :: M, b :: M) where M <: AbsOreMonomial
    return lt(ord, a, b)
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