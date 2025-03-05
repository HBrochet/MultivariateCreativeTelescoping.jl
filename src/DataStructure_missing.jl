function Base.append!(s :: SortedSet{M}, v :: Vector{M}) where M
    for m in v 
        push!(s,m)
    end
    return nothing 
end

@inline function Base.popfirst!(m::SortedSet)
    i = DataStructures.beginloc(m.bt)
    i == 2 && throw(BoundsError())
    k = m.bt.data[i].k
    delete!(m.bt, i)
    return k
end


@inline function poplast!(m::SortedSet)
    i = DataStructures.endloc(m.bt)
    i == 2 && throw(BoundsError())
    k = m.bt.data[i].k
    delete!(m.bt, i)
    return k
end