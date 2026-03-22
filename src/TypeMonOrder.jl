abstract type AbsMonomialOrder{N} <: Base.Order.Ordering end

struct GeneratedMonomialOrder{N,R,C,B} <: AbsMonomialOrder{N} end

@inline function _order_weight_at(B::Tuple, R::Int, i::Int, j::Int)
    return B[(j - 1) * R + i]
end

function order_weight_matrix(A :: Vector{Vector{I}}, ::Val{M}) where {I <: Integer, M <: AbsOreMonomial}
    R = length(A)
    C = nvars(M)
    data = ntuple(k -> begin
        i = (k - 1) % R + 1
        j = (k - 1) ÷ R + 1
        Int16(A[i][j])
    end, R * C)
    return SMatrix{R,C,Int16,R*C}(data)
end

nvars(:: AbsMonomialOrder{N}) where N = N

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
    B = order_weight_matrix(A, Val(M))
    return GeneratedMonomialOrder{N,size(B, 1),size(B, 2),Tuple(B)}()
end

@generated function lt(::GeneratedMonomialOrder{N,R,C,B}, a :: OreMonVE{N,E}, b :: OreMonVE{N,E}) where {N,R,C,B,E <: Integer}
    rows = Any[]
    for i in 1:R
        terms_a = [:(t1 += a[$j] * Int16($(_order_weight_at(B, R, i, j)))) for j in 1:C if _order_weight_at(B, R, i, j) != 0]
        terms_b = [:(t2 += b[$j] * Int16($(_order_weight_at(B, R, i, j)))) for j in 1:C if _order_weight_at(B, R, i, j) != 0]
        push!(rows, quote
            t1 = Int16(0)
            $(terms_a...)
            t2 = Int16(0)
            $(terms_b...)
            if t1 < t2
                return true
            elseif t1 > t2
                return false
            end
        end)
    end
    return Expr(:block, rows..., :(return false))
end

@generated function max_deg_block(::GeneratedMonomialOrder{N,R,C,B}, a :: OreMonVE{N,E}) where {N,R,C,B,E <: Integer}
    terms = [:(acc += a[$j] * Int16($(_order_weight_at(B, R, 1, j)))) for j in 1:C if _order_weight_at(B, R, 1, j) != 0]
    return quote
        acc = Int16(0)
        $(terms...)
        return Int(acc)
    end
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
