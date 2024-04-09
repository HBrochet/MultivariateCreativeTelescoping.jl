"""
    parse_vector_OrePoly(s :: String, A :: OreAlg)

Return a vector of OrePoly corresponding to the parseable string s in the algebra A.
"""
function parse_vector_OrePoly(s :: String, A :: OreAlg)
    expr = Meta.parse(s)
    if expr.head != :vect 
        error("input is not a parsable vector")
    end
    return [expr_to_OreAlg(i,A) for i in expr.args]
end

"""
    parse_OrePoly(s :: String, A :: OreAlg)

Return an OrePoly corresponding to the parseable string s in the algebra.
"""
function parse_OrePoly(s :: String, A :: OreAlg)
    expr = Meta.parse(s)
    r = expr_to_OreAlg(expr, A)
    return normalize!(r,A)
end

function expr_to_OreAlg(expr :: Symbol, A :: OreAlg)
    if haskey(A.strvar_to_indexp, "$expr")
        i = A.strvar_to_indexp["$expr"]
        return makepoly(one(ctx(A)),makemon(i,A))
    end
    if haskey(A.ratvars, "$expr")
        c = A.ratvars["$expr"]
        return makepoly(c, makemon(0,A))
    end
    error("pb while parsing symbol $(expr) not recognized")
end

function expr_to_OreAlg(expr :: Number, A :: OreAlg)
    return makepoly(convertn(expr,ctx(A)),makemon(-1,A))
end


function expr_to_OreAlg(expr :: Expr, A :: OreAlg)
    if expr.head == :call 
        return  expr_to_OreAlg_call(expr, A)
    elseif expr.head == :macrocall
        return expr_to_OreAlg_macrocall(expr,A)
    end
end

function expr_to_OreAlg_macrocall(expr :: Expr, A :: OreAlg)
    if (expr.args[1] == GlobalRef(Core, Symbol(Core.var"@big_str"))) || (expr.args[1] ==  GlobalRef(Core, Symbol(Core.var"@int128_str"))) ||  (expr.args[1] ==  GlobalRef(Core, Symbol(Core.var"@uint128_str")))
        return makepoly(convertn(parse(BigInt,expr.args[3]),ctx(A)),makemon(-1,A))
    end
    error("error while parsing arg $(expr.args[1]) not recognised")
end


function expr_to_OreAlg_call(expr :: Expr, A :: OreAlg)
    op = expr.args[1]
    a = expr_to_OreAlg(expr.args[2],A)
    len = length(expr.args)
    if op == :- && len == 2
        return mul(opp(one(ctx(A)),ctx(A)),a,A)
    elseif op == :+
        for i in 3:len
            a = add!(a,expr_to_OreAlg(expr.args[i],A),A)
        end
        return a
    elseif op == :*
        for i in 3:len
            a = mul(a,expr_to_OreAlg(expr.args[i],A),A)
        end
        return a
    end

    # there are only 3 arguments for the latter cases
    b = expr_to_OreAlg(expr.args[3],A)
    if op == :- 
        sub!(a,b,A)
        return a
    elseif op == :/ || op == ://
        #b should be a constant
        coeffs(b)[1] = inv(coeff(b,1), ctx(A))
        return mul(a, b, A)
    elseif op == :^ 
        #b should be an integer
        powa = deepcopy(a)
        bound = coeff(b,1)
        if bound isa RatFunQQ || bound isa RatFunModp || bound isa RatFunModP
            if is_zero(bound)
                bound = 0
            else 
                bound = Nemo.coeff(numerator(bound),1)
            end
        end 
        if bound isa ZZRingElem
            bound = bound.d
        elseif bound isa QQFieldElem
            bound = numerator(bound)
        elseif bound isa fpFieldElem
            bound = bound.data
        elseif !(bound isa Integer) # à changer
            bound = Int(evaluate(bound, [1 for i in 1:A.nrdv]).num)
        end
        if bound == 0
            return one(A)
        end
        for i in 2:bound # fastexp would be better but w/e
            powa = mul(a,powa,A)
        end
        return powa
    else 
        error("operator $(op) not recognised/supported")
    end
end

function printmon(m :: OreMonVE{N,E}, A :: OreAlg) where {N,E}
    for i in 1:N 
        if m[i] == E(1) 
            print(A.indexp_to_strvar[i])
        elseif m[i] != E(0) 
            print(A.indexp_to_strvar[i],"^",m[i])
        end
    end
end

"""
    prettyprint(p :: OrePoly,A ::OreAlg)

Print the OrePoly p.
"""
function prettyprint(p :: OrePoly,A ::OreAlg)
    if length(p) == 0 
        println("(0)")
        return
    end
    print("(",coeff(p,1),")")
    printmon(mon(p,1),A)
    for i in 2:length(p)
        print(" + (", coeff(p,i),")")
        printmon(mon(p,i),A)
    end
    println()
end

"""
    prettyprint(v :: Vector{OrePoly{K,M}},A :: OreAlg)

Print the vector of OrePoly v.
"""
function prettyprint(v :: Vector{OrePoly{K,M}},A :: OreAlg) where {K,M}
    println("vector of ", length(v), " OrePoly")
    for i in 1:length(v)
        prettyprint(v[i],A)
    end
end

function prettyprint(v :: Vector{Vector{OrePoly{K,M}}},A :: OreAlg) where {K,M}
    println("vector of ", length(v), " of vector OrePoly")
    for i in 1:length(v)
        prettyprint(v[i],A)
    end

end


function mystring(m :: OreMonVE{N,E}, A :: OreAlg) where {N,E}
    s = "" 
    for i in 1:N 
        if m[i] == E(1) 
            s *= string(A.indexp_to_strvar[i])*"*"
        elseif m[i] != E(0) 
            s *= string(A.indexp_to_strvar[i])*"^"*string(m[i])*"*"
        end
    end
    return s*"1"
end

"""
    mystring(p :: OrePoly, A:: OreAlg)

Returns a string representing the Ore polynomial p.
"""
function mystring(P :: OrePoly, A:: OreAlg)
    s= "" 
    if length(P) > 0
        s *= "("*string(coeff(P,1))*")*"*mystring(mon(P,1),A)
    end
    for i in 2:length(P)
        s *= " + ("*string(coeff(P,i))*")*"*mystring(mon(P,i),A)
    end
    return s
end

"""
    mystring(v :: Vector{OrePoly{K,M}}, A:: OreAlg)

Returns a string representing the vector of Ore polynomials v.
"""
function mystring(P :: Vector{OrePoly{K,M}}, A:: OreAlg) where {K,M}
    s= "[" 
    if length(P) > 0
        s *= mystring(P[1],A)
    end
    for i in 2:length(P)
        s *= ", "*mystring(P[i],A)
    end
    return s*"]"
end