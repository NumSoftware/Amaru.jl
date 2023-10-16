# Simple ordered set with getindex and no collision support
struct OSet{T}
    vals::Vector{T}
    dict::Dict{UInt, T}
end

function OSet()
    return OSet{Any}(Vector{Any}(), Dict{UInt, Any}())
end

function OSet{T}() where T
    return OSet{T}(Vector{T}(), Dict{UInt, T}())
end

Base.length(os::OSet{T}) where T = length(os.vals)

# Implement the iteration interface
function Base.iterate(os::OSet{T}, state=1) where T
    if state > length(os.vals)
        return nothing
    else
        return (os.vals[state], state + 1)
    end
end

function Base.push!(os::OSet{T}, val::T) where T
    hash_value = hash(val)
    haskey(os.dict, hash_value) && return val
    push!(os.vals, val)
    os.dict[hash_value] = val
end

function Base.getindex(os::OSet{T}, index::Int) where T
    if 1 <= index <= length(os.vals)
        return os.vals[index]
    else
        throw(BoundsError(os, index))
    end
end

function Base.delete!(os::OSet{T}, element::T) where T
    hash_value = hash(element)
    
    # Check if the hash value is in the dictionary
    if haskey(os.dict, hash_value)
        # Remove the element from the vector
        filter!(x -> x != element, os.vals)
        
        # Remove the element from the dictionary
        delete!(os.dict, hash_value)
        return os
    else
        error("Element not found in the ordered set.")
    end
end

# Implement the get function
function Base.get(os::OSet{T}, key::T, default) where T
    get(os.dict, hash(key), default)
end
