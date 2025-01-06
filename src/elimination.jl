const NmodF4Matrix{I, T, Tbuf} = F4Matrix{M,I,T,Alg} where {T,Tbuf,C <: AbsContextCoeff{T, Tbuf},M,O,Alg <:OreAlg{T,C,M,O}}




@inline function criticalloop!(ctx :: NmodLikeΓ{T, Tbuf},
                               buffer :: Vector{Tbuf},
                               red :: OrePoly,
                               mult :: T) where {T, Tbuf}
    maxl = length(red)-3
    l = 2
    @inbounds while l <= maxl
        m = mon(red, l)
        c = coeff(red, l)
        buffer[m] = submul(buffer[m], mult, c, ctx)
        m = mon(red, l+1)
        c = coeff(red, l+1)
        buffer[m] = submul(buffer[m], mult, c, ctx)
        m = mon(red, l+2)
        c = coeff(red, l+2)
        buffer[m] = submul(buffer[m], mult, c, ctx)
        m = mon(red, l+3)
        c = coeff(red, l+3)
        buffer[m] = submul(buffer[m], mult, c, ctx)
        l += 4
    end

    @inbounds for l in l:length(red)
        m = mon(red, l)
        c = coeff(red, l)
        buffer[m] = submul(buffer[m], mult, c, ctx)
    end
end


function fillbuffer!(buffer :: Vector{Tbuf}, row :: OrePoly{K, I},ctx ::NmodLikeΓ{T, Tbuf}) where {K,I,T,Tbuf}
    fill!(buffer, zero(Tbuf, ctx))
    @inbounds for (c,j) in row
        buffer[j] = c
    end
end


function savebuffer!(ctx :: NmodLikeΓ{T, Tbuf},
                    buffer :: Vector{Tbuf},  # coefficients must be normalized
                    r :: OrePoly{T,I},
                    firstterm,
                    normalize :: Bool) where {T, Tbuf,I}


    row = OrePoly(T[],I[])

    if normalize
        @inbounds mult = inv(deflate(normal(buffer[firstterm], ctx), ctx), ctx)
        @inbounds for j in firstterm:length(buffer)
            @assert normal(buffer[j], ctx) == buffer[j]
            c = deflate(buffer[j], ctx)
            iszero(c,ctx) && continue
            push!(row.mons, I(j))
            push!(row.coeffs, mul(mult, c,ctx)) # here, the buffer MUST be normalized
            # @assert ismonic(row)
            #if (BigInt(row.co[end])*buffer[firstterm]-c)%ctx.char != 0
            #    error((buffer[firstterm], mult, c, ctx.char))
            #end
        end
    else
        @inbounds for j in firstterm:length(buffer)
            @assert normal(buffer[j], ctx) == buffer[j]
            c = buffer[j]
            iszero(c,ctx) && continue
            push!(row.mons, I(j))
            push!(row.coeffs, deflate(c, ctx))
        end
    end
    return row
end


function elementary_reduction_in_buffer(mx::NmodF4Matrix{I, T, Tbuf},
                                        buffer::Vector{Tbuf},
                                        j
                                        ) where {I, T, Tbuf}

    buffer[j] = normal(buffer[j],ctx(mx.A))
    if iszero(buffer[j],ctx(mx.A))
        return true
    end

    piv = mx.pivots[j]
    if piv == 0
        return false
    end

    reducer = mx.rows[piv]

    #@assert ismonic(reducer)
    #@assert mon(reducer,1) == j
    globalstats.counters[:f4_line_reductions] += 1

    globalstats.counters[:f4_field_operations] += length(reducer)

    mult = deflate(buffer[j],ctx(mx.A))
    buffer[j] = zero(T,ctx(mx.A))
    # Critical loop
    # globalstats.timings[:critical_loop] +=
    #     @elapsed 
    criticalloop!(ctx(mx.A), buffer, reducer, mult)
    return true
end

function reduce!(mx::NmodF4Matrix{I, T, Tbuf}) where {T, Tbuf, I}
    #@assert all([i == 0 ? true : ismonic(mx.rows[i]) for i in mx.pivots])

    buffer = Vector{Tbuf}(undef, mx.nbcolumns)
    # iter = 0 
    for (i, row) in enumerate(mx.rows)
        if isempty(row) || ispivot(mx, i)
            continue
        end

        fillbuffer!(buffer, row, ctx(mx.A))

        leadingterm = 0

        for j in mon(row,1):mx.nbcolumns
            if !elementary_reduction_in_buffer(mx, buffer, j)
                leadingterm = j
                break
            end
        end

        # Do it as savebuffer! expects it
        for j in leadingterm+1:mx.nbcolumns
            buffer[j] = normal(buffer[j], ctx(mx.A))
        end

        if leadingterm > 0
            # nonzero reduction
            if !(i ∈ mx.donotpivot)  # && mx.εcolumn > leadingterm
                # new pivot
                mx.rows[i] = savebuffer!(ctx(mx.A), buffer, row, leadingterm, true)
                mx.pivots[leadingterm] = i
                push!(mx.newpivots, i)       
            else
                mx.rows[i] = savebuffer!(ctx(mx.A), buffer, row, leadingterm, false)      
            end
        else
            mx.rows[i] = OrePoly(T[],I[]) # pourquoi a-t-on besoin de ça ? 
        end
    end

    return
end


function reducepivots!(mx :: NmodF4Matrix{I, T, Tbuf}) where {T, Tbuf, I}

    buffer = Vector{Tbuf}(undef, mx.nbcolumns)
    # Going bottom up is important for performance
    @inbounds for p in mx.nbcolumns:-1:1
        mx.pivots[p] == 0 && continue
        row = mx.rows[mx.pivots[p]]
        #@assert ismonic(row)
        length(row) <= 1 && continue

        fillbuffer!(buffer, row,ctx(mx.A))
        for j in mon(row, 2):mx.nbcolumns
            elementary_reduction_in_buffer(mx, buffer, j)
        end

        mx.rows[mx.pivots[p]] = savebuffer!(ctx(mx.A), buffer, row, p, false)
    end
end



# #. Interreduction of Gröbner bases

function interreduce(alg :: OreAlg,
                      pols :: Vector{OrePoly{I, T}},
                      fullreduction :: Bool = true
                      ) where {I, T}

    isempty(pols) && return pols
    for p in pols
        if length(p) == 0 
            error("zero vector in basis in interreduce")
        end
    end

    mx = interreductionmx(alg, pols)
    if fullreduction
        reducepivots!(mx)
    end
    reduce!(mx)


    # It is important here that a polynomial in `pols` is not chosen as pivots
    # if its leading term is reducible.
    pols = [row(mx, i) for i in mx.inputrows]

    filter!(p -> !isempty(p), pols)
    sort!(pols, order=order(alg))

    return pols
end