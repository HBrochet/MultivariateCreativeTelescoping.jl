


function add_echelon!(echelon :: Vector{OrePoly{T,M}}, P :: OrePoly{T,M}, A :: OreAlg; augmented = false, echelonvect = Vector{T}[], vect = T[],unitary = true ) where {T,M}
    if length(P) == 0
        globalstats.counters[:mct_red_zero] += 1
        return 
    end
    if augmented
        if unitary
            for i in 1:length(vect)
                vect[i] = mul(vect[i], inv(coeff(P,1),ctx(A)),ctx(A))
            end
        end

        for v in echelonvect 
            push!(v, zero(ctx(A)))
        end
    end
    if unitary
        makemonic!(P,A)
    end
    for i in length(echelon):-1:1 
        if lt(order(A),mon(P,1),mon(echelon[i],1))
            insert!(echelon,i+1,P)
            if augmented 
                insert!(echelonvect,i+1,vect)
            end
            return
        end 
    end 
    insert!(echelon,1,P)
    if augmented
        insert!(echelonvect,1,vect)
    end
end

function reduce_with_echelon!(echelon :: Vector{OrePoly{T,M}}, P :: OrePoly{T,M}, A :: OreAlg, geob :: GeoBucket,tmp_poly ::ReuseOrePoly) where {T,M}
    init_GeoBucket!(geob,P)
    return reduce_with_echelon!(echelon,A,geob,tmp_poly)
end

function reduce_with_echelon!(echelon :: Vector{OrePoly{T,M}}, A :: OreAlg, geob :: GeoBucket,tmp_poly ::ReuseOrePoly) where {T,M}
    e = 1 
    while true
        c,m = lt(geob,A)
        if iszero(c)
            break
        end
        div = false
        for i in e:length(echelon)
            me = mon(echelon[i],1)
            if lt(order(A), me, m)
                e = i 
                break
            elseif m == me
                div = true
                geob = add!(geob, opp(c,ctx(A)),echelon[i],A)
                e = i + 1
                break
            end
        end
        if !div
            push!(tmp_poly,c,m)
            rem_lt!(geob,m)
        end
    end
    return copy_to_OrePoly!(tmp_poly,A)
end




# function reduce_with_echelon_augmented!(echelon :: Vector{OrePoly{T,M}}, P :: OrePoly{T,M}, A :: OreAlg, echelonvect :: Vector{Vector{T}}) where {T,M}
#     res = P
#     vect = T[zero(ctx(A)) for i in 1:length(echelon)]
#     push!(vect,one(ctx(A)))

#     r = 1 
#     e = 1 
#     while r <= length(res)
#         c,m = res[r]
#         div = false
#         for i in e:length(echelon)
#             if lt(order(A), mon(echelon[i],1), m)
#                 e = i 
#                 break
#             elseif m == mon(echelon[i],1)
#                 div = true
#                 res = sub!(res,mul(c,echelon[i],A),A)
#                 for l in 1:length(echelonvect[i])
#                     vect[l] = sub(vect[l],mul(c,echelonvect[i][l], ctx(A)),ctx(A))
#                 end
#                 e = i + 1
#                 break
#             end
#         end
#         if !div
#             r += 1
#         end
#     end
#     return res, vect 
# end