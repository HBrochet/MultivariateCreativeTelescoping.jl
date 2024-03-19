function parse_vector_OrePoly(s :: String, A::OreAlg)
    expr = Meta.parse(s)
    if expr.head != :vect 
        error("input is not a parsable vector")
    end
    return [expr_to_OreAlg(i,A) for i in expr.args]
end

function parse_OrePoly(s :: String, A :: OreAlg)
    expr = Meta.parse(s)
    r = expr_to_OreAlg(expr, A)
    return normalize!(r,A)
end

function expr_to_OreAlg(expr :: Symbol, A :: OreAlg{})
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
        if bound isa RatFunQQ || bound isa RatFunModp
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
        error("operator not recognised/supported")
    end
end

function printmon(m :: OreMonVE{N,E}, A :: AbsOreAlgebra) where {N,E}
    for i in 1:N 
        if m[i] == E(1) 
            print(A.indexp_to_strvar[i])
        elseif m[i] != E(0) 
            print(A.indexp_to_strvar[i],"^",m[i])
        end
    end
end


function prettyprint(P :: OrePoly,A ::AbsOreAlgebra)
    if length(P) > 0
        print("(",coeff(P,1),")")
        printmon(mon(P,1),A)
    end
    for i in 2:length(P)
        print(" + (", coeff(P,i),")")
        printmon(mon(P,i),A)
    end
    println()
end

function prettyprint(v :: Vector{OrePoly{K,M}},A :: AbsOreAlgebra) where {K,M}
    println("vector of ", length(v), " OrePoly")
    for i in 1:length(v)
        prettyprint(v[i],A)
    end
end


function mystring(m :: OreMonVE{N,E}, A :: AbsOreAlgebra) where {N,E}
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