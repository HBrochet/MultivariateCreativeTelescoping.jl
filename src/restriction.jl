function fourrier(p :: OrePoly{C,M}, A::Alg; inv ::Bool = false)
    res = zero(A)
    E = eltype_mo(A)
    for i in 1:length(p)
        m = [mon(p,i)[j] for j in 1:A.nrdv]
        append!(m, [mon(p,i)[j] for j in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv])
        append!(m,[E(0) for j in 1:A.npdv+A.npv])
        append!(m,[mon(p,i)[j] for j in A.ndv+2*A.npdv+1:N])
        m2 = [E(0) for j in 1:A.nrdv+A.npdv]
        append!(m2,[mon(p,i)[j] for j in A.nrdv+1:A.nrdv+A.npdv])
        append!(m2,[E(0) for j in A.nrdv+A.npdv+1:N])
        p1 = makepoly(coeff(p,i),makemon(m))
        if inv 
            s = sum([mon(p,i)[j] for j in A.ndv+A.npdv+1:A.ndv+2*A.npdv])
        else
            s = sum([mon(p,i)[j] for j in A.ndv+1:A.ndv+A.npdv])
        end
        if s & E(1) == E(1)
            p2 = makepoly(opp(one(ctx(A)),ctx(A)), makemon(m2))
        else
            p2 = makepoly(one(ctx(A)), makemon(m2))
        end

        res = add!(res, mul(p1,p2,A),A)
    end
    return res
end

function fourrier(g :: Vector{OrePoly{C,M}}, A:: Alg; inv = false) where {C,M,Alg <: AbsOreAlgebra}
    res = OrePoly{C,M}[]
    for p in g 
        push!(res,fourrier(p,A,inv=inv))
    end
    return res 
end

function change_algebra_invpoleps(A :: OreAlg)
    inp = deepcopy(A.inp)
    inp.order = "lex _Deps > "*inp.order
    push!(A.ratdiffvars,"_eps")
    push!(v,A.ndv+A.npdv+1,"_Deps")
    return OreAlg(inp)
end

function initialize_ann(pol :: OrePoly{T,M},A :: AbsOreAlgebra) where {T,M}
    ann = OrePoly{T,M}[]
    eps = makepoly(one(ctx(A)),makemon(A.nrdv + A.npdv*2,A))
    Deps = makepoly(one(ctx(A)),makemon(A.nrdv + A.npdv,A))
    push!(ann, mul(Deps,add(eps,pol,A),A))

    for i in 1:A.npdv-1
        tmp = makepoly(one(ctx(A)),makemon(A.ndv + i,A))
        sub!(tmp,mul(diff(pol,i,A),Deps,A),A)
        push!(ann,tmp)
    end
    return ann
end