const RingF4Matrix{I, T, Tbuf} =
    F4Matrix{M,I,T,Alg} where {T,Tbuf,C <: RingCtx{T, Tbuf},M,O,Alg <:OreAlg{T,C,M,O}}

function fillbuffer!(buffer::Vector{T}, row::OrePoly{T,I}, ctx::RingCtx{T,Tbuf}) where {T,Tbuf,I}
    @inbounds for j in eachindex(buffer)
        buffer[j] = zero(T, ctx)
    end
    @inbounds for (c, j) in row
        buffer[j] = c
    end
end

function savebuffer!(ctx::RingCtx{T,Tbuf},
                    buffer::Vector{T},
                    r::OrePoly{T,I},
                    firstterm,
                    normalize::Bool) where {T,Tbuf,I}
    row = OrePoly(T[], I[])
    @inbounds for j in firstterm:length(buffer)
        c = buffer[j]
        iszero(c, ctx) && continue
        push!(row.mons, I(j))
        push!(row.coeffs, c)
    end
    return row
end

function scale_buffer!(ctx::RingCtx{T,Tbuf}, buffer::Vector{T}, mult::T) where {T,Tbuf}
    isone(mult, ctx) && return buffer
    @inbounds for j in eachindex(buffer)
        mul!(buffer[j], mult, ctx)
    end
    return buffer
end

function elementary_reduction_in_buffer(mx::RingF4Matrix{I, T, Tbuf},
                                        buffer::Vector{T},
                                        j::Int,
                                        param::F4Param) where {I,T,Tbuf}
    ctxA = ctx(mx.A)
    c = buffer[j]
    iszero(c, ctxA) && return true

    piv = mx.pivots[j]
    piv == 0 && return false

    reducer = mx.rows[piv]
    pivco = coeff(reducer, 1)

    gcd_ = gcd(pivco, c) 
    co = divexact(c,gcd_,ctxA)

    scale_buffer!(ctxA, buffer, divexact(pivco,gcd_,ctxA))
    buffer[j] = zero(T, ctxA)

    if stat(param)
        globalstats.counters[:f4_line_reductions] += 1
        globalstats.counters[:f4_field_operations] += length(reducer)
    end


    tmp = parent(co)()
    @inbounds for k in 2:length(reducer)
        m = mon(reducer, k)
        ck = coeff(reducer, k)
        Nemo.mul!(tmp,co,ck)
        sub!(buffer[m], tmp, ctxA)
    end
    return true
end

function reduce!(mx::RingF4Matrix{I, T, Tbuf}, param::F4Param) where {T,Tbuf,I}
    ctxA = ctx(mx.A)
    buffer = Vector{T}(undef, mx.nbcolumns)

    for i in mx.inputrows
        row = mx.rows[i]
        if isempty(row) || ispivot(mx, i)
            continue
        end

        fillbuffer!(buffer, row, ctxA)

        leadingterm = 0
        for j in mon(row, 1):mx.nbcolumns
            if !elementary_reduction_in_buffer(mx, buffer, j, param)
                leadingterm = j
                break
            end
        end

        if leadingterm > 0
            for k in leadingterm+1:mx.nbcolumns
                elementary_reduction_in_buffer(mx, buffer, k, param)
            end

            newrow = savebuffer!(ctxA, buffer, row, leadingterm, false)
            primitive_part!(newrow, mx.A)
            mx.rows[i] = newrow
            mx.pivots[leadingterm] = i
            push!(mx.newpivots, i)
        else
            mx.rows[i] = OrePoly(T[], I[])
        end
    end
    return nothing
end

function reducepivots!(mx::RingF4Matrix{I, T, Tbuf}, param::F4Param) where {T,Tbuf,I}
    ctxA = ctx(mx.A)
    buffer = Vector{T}(undef, mx.nbcolumns)

    @inbounds for p in mx.nbcolumns:-1:1
        mx.pivots[p] == 0 && continue
        row = mx.rows[mx.pivots[p]]
        length(row) <= 1 && continue

        fillbuffer!(buffer, row, ctxA)
        for j in mon(row, 2):mx.nbcolumns
            elementary_reduction_in_buffer(mx, buffer, j, param)
        end

        newrow = savebuffer!(ctxA, buffer, row, p, false)
        primitive_part!(newrow, mx.A)
        mx.rows[mx.pivots[p]] = newrow
    end
end

function reducenewpivots!(mx::RingF4Matrix{I, T, Tbuf}, param::F4Param) where {T,Tbuf,I}
    ctxA = ctx(mx.A)
    buffer = Vector{T}(undef, mx.nbcolumns)

    piv = [mon(mx.rows[p], 1) for p in mx.newpivots]
    sort!(piv, rev = true)

    @inbounds for p in piv
        row = mx.rows[mx.pivots[p]]
        length(row) < 2 && continue

        fillbuffer!(buffer, row, ctxA)
        for j in mon(row, 2):mx.nbcolumns
            elementary_reduction_in_buffer(mx, buffer, j, param)
        end

        newrow = savebuffer!(ctxA, buffer, row, p, false)
        primitive_part!(newrow, mx.A)
        mx.rows[mx.pivots[p]] = newrow
    end
    return nothing
end
