struct Signature{N} 
    mon :: OreMonVE{N,Int16}
    ind :: Int64
end


# Signature(v :: Vector,i :: Int) = Signature{length(v)}(OreMonVE(SVector{length(v)}(map(x -> Int16(x),v))),i)

function Base.lcm(a :: Signature{N}, b :: Signature{N}) where N
    return SVector{N,Int16}([max(a.mon[i],b.mon[i]) for i in 1:N])
end

struct SigPair{N,T, M <: OreMonVE{N,Int16}} 
    op :: OrePoly{T, M}
    sig :: Signature{N}
end
Base.show(io ::IO, ::SigPair) = print(io,"SigPair")
Base.show(io ::IO, ::Vector{SigPair{N,T,M}}) where {N,T,M} = print(io,"Vector of SigPairs")

Base.zero(s::SigPair) = zero(s.op)
Base.iszero(s::SigPair) = length(s.op) == 0
Base.iszero(s::SigPair, :: OreAlg) = length(s.op) == 0

lm(s ::SigPair) =  lm(s.op)
lc(s::SigPair) = lc(s.op)

function mon(s :: SigPair, i :: Int)
    return mon(s.op,i)
end

function coeff(s :: SigPair, i :: Int)
    return coeff(s.op, i)
end

function Base.length(s :: SigPair)
    return length(s.op)
end

function zeroSigPair(A ::OreAlg)
    N = nvars(A)
    m = OreMonVE(SVector{N,Int16}(Int16(0) for i in 1:N))
    return SigPair(zero(A),Signature{N}(m,0))
end




# function iscompatible(a :: SigPair, b :: SigPair, A:: OreAlg)
#     return iscompatible(mon(a,1),mon(b,1),A)
# end

function iscompatible(a :: Signature, b :: Signature, A:: OreAlg)
    return iscompatible(a.mon,b.mon,A)
end


abstract type SigOrder{O <: AbsMonomialOrder} <: Base.Order.Ordering end

struct POT{O} <: SigOrder{O} 
    ord :: O
end
order(so :: POT) = so.ord

struct TOP{O} <: SigOrder{O} 
    ord :: O
end
order(so :: TOP) = so.ord


function lt(so :: TOP{O}, a :: Signature{N}, b :: Signature{N}) where {N, O <: AbsMonomialOrder}
    if lt(order(so),a.mon, b.mon) 
        return true
    elseif lt(order(so),b.mon,a.mon)
        return false
    else
        return a.ind < b.ind
    end
end

function lt(so :: POT{O}, a :: Signature{N}, b :: Signature{N}
                      ) where {N, O <: AbsMonomialOrder}
    if a.ind < b.ind 
        return true
    elseif b.ind < a.ind 
        return false
    else
        return Base.order.lt(order(so),a.mon, b.mon) 
    end
end

function Base.Order.lt(so :: SigOrder{O}, a :: Signature{N}, b :: Signature{N}) where {N, O <: AbsMonomialOrder}
    return lt(so,a,b)
end



function divide(s :: Signature{N}, ss :: Signature{N},A :: OreAlg) where N
    if s.ind != ss.ind
        return false
    end
    return divide(s.mon,ss.mon, A)
end


function Base.:(==)(s :: Signature{N}, ss :: Signature{N}) where N
    return s.ind == ss.ind && s.mon == ss.mon
end

function Base.:(==)(s :: SigPair, ss :: SigPair)
    return  s.sig == ss.sig && s.op == ss.op
end

function Base.copy(s :: SigPair)
    return SigPair(copy(s.op),s.sig)
end

function emptysigpair(A :: OreAlg)
    N = nvars(A)
    sig = Signature([0 for i in 1:N],0)
    wp = emptypol(A)
    return SigPair(wp, sig)
end

function mul(v :: OreMonVE,f :: SigPair{N,K,M}, A :: OreAlg)  where {M,N,K}
    return SigPair(mul(v,f.op,A), Signature(v*f.sig.mon,f.sig.ind))
end


function mul(a :: OreMonVE, sig :: Signature{N})   where N
    return Signature(sig.mon*a, sig.ind)
end

function mul(c :: K,f :: SigPair{N,K,M}, A :: OreAlg)  where {M,N,K}
    return SigPair(mul(c,f.op,A), f.sig)
end

function mul!(v :: OreMonVE,f :: SigPair{N,K,M}, A :: OreAlg)  where {M,N,K}
    return SigPair(mul!(v,f.op,A), Signature(v*f.sig.mon,f.sig.ind))
end


function mul!(c :: K,f :: SigPair{N,K,M}, A :: OreAlg)  where {M,N,K}
    return SigPair(mul!(c,f.op,A), f.sig)
end

function prettyprint(s :: Signature,A :: OreAlg)
    print(s.ind,", ")
    printmon(s.mon,A)
    println()
end    
function prettyprint(g :: Vector{SigPair{N,K,M}}, A :: OreAlg) where {M, N, K} 
    println("Vector of WeylPol of size $(length(g))")
    for i in 1:length(g)
        prettyprint(g[i].sig,A)
        prettyprint(g[i].op,A)
    end
end

function isholonomic(gb :: Vector{SigPair{N,C,M}},A :: OreAlg) where {N,C,M}
    #println("testing holonomy")
    set = [i for i in 1:A.nrdv+2*A.npdv]
    n = A.npdv
    # returns an iterator over every subset of set of size n+1
    it = powerset(set,n+1,n+1) 

    for subset in it 
        gfound = false
        for i in 1:length(gb)
            if A.nlv > 0 && mon(gb[i],1)[A.nrdv+2*A.npdv+A.npv+1] > 0
                continue
            end 
            # test if gb[i] verifies the criterion
            b = true  
            for j in 1:N
                if !(j in subset) && mon(gb[i],1)[j] != 0
                    b = false
                    break
                end
            end
            if b 
                gfound = true
                break
            end
        end
        if !gfound 
            # println(subset)
            return false
        end
    end
    return true
end