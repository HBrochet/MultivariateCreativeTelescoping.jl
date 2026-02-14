const NmodF4Matrix{I, T, Tbuf} =
    F4Matrix{M,I,T,Alg} where {T,Tbuf,C <: NmodLikeΓ{T, Tbuf},M,O,Alg <:OreAlg{T,C,M,O}}

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

    @inbounds for k in l:length(red)
        m = mon(red, k)
        c = coeff(red, k)
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
            c = deflate(buffer[j], ctx)
            iszero(c,ctx) && continue
            push!(row.mons, I(j))
            push!(row.coeffs, mul(mult, c,ctx))
        end
    else
        @inbounds for j in firstterm:length(buffer)
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
                                        j :: Int,
                                        param :: F4Param
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
    if stat(param)
        globalstats.counters[:f4_line_reductions] += 1
        globalstats.counters[:f4_field_operations] += length(reducer)
    end

    mult = deflate(buffer[j],ctx(mx.A))
    buffer[j] = zero(T,ctx(mx.A))
    criticalloop!(ctx(mx.A), buffer, reducer, mult)
    return true
end

function reduce!(mx::NmodF4Matrix{I, T, Tbuf},param :: F4Param) where {T, Tbuf, I}
    buffer = Vector{Tbuf}(undef, mx.nbcolumns)
    for i in mx.inputrows
        row = mx.rows[i]
        if isempty(row) || ispivot(mx,i)
            continue
        end

        fillbuffer!(buffer, row, ctx(mx.A))

        leadingterm = 0

        for j in mon(row,1):mx.nbcolumns
            if !elementary_reduction_in_buffer(mx, buffer, j,param)
                leadingterm = j
                break
            end
        end

        if leadingterm > 0
            for k in leadingterm+1:mx.nbcolumns
                elementary_reduction_in_buffer(mx, buffer, k,param)
            end
            mx.rows[i] = savebuffer!(ctx(mx.A), buffer, row, leadingterm, true)
            mx.pivots[leadingterm] = i
            push!(mx.newpivots, i)
        else
            mx.rows[i] = OrePoly(T[],I[])
        end
    end
    return
end


function reducepivots!(mx :: NmodF4Matrix{I, T, Tbuf},param :: F4Param) where {T, Tbuf, I}

    buffer = Vector{Tbuf}(undef, mx.nbcolumns)
    @inbounds for p in mx.nbcolumns:-1:1
        mx.pivots[p] == 0 && continue
        row = mx.rows[mx.pivots[p]]
        length(row) <= 1 && continue

        fillbuffer!(buffer, row,ctx(mx.A))
        for j in mon(row, 2):mx.nbcolumns
            elementary_reduction_in_buffer(mx, buffer, j,param)
        end

        mx.rows[mx.pivots[p]] = savebuffer!(ctx(mx.A), buffer, row, p, false)
    end
end

function reducenewpivots!(mx :: NmodF4Matrix{I, T, Tbuf},param :: F4Param) where {T, Tbuf, I}
    buffer = Vector{Tbuf}(undef, mx.nbcolumns)

    piv = [mon(mx.rows[p],1) for p in mx.newpivots]
    sort!(piv,rev=true)

    @inbounds for p in piv
        row = mx.rows[mx.pivots[p]]
        if length(row) < 2
            continue
        end
        fillbuffer!(buffer, row, ctx(mx.A))
        for j in mon(row,2):mx.nbcolumns
            elementary_reduction_in_buffer(mx, buffer, j,param)
        end
        mx.rows[mx.pivots[p]] = savebuffer!(ctx(mx.A), buffer, row, p, false)
    end
end

