mutable struct ArrayList
    arrays::Array{AbstractArray,1}
    inuse::Array{Bool, 1}

    function ArrayList(T::Type, dims::Tuple)
        arrays = Array{T, length(dims)}[]
        inuse  = Bool[]
        return new(arrays, inuse)
    end
end

mutable struct ThreadPool  # Pool per thread
    dict::Dict{Tuple, ArrayList}
    function ThreadPool()
        return new(Dict{Tuple, ArrayList}())
    end
end

mutable struct ArrayPool
    tpools::Array{ThreadPool,1}
    function ArrayPool()
        tpools = ThreadPool[ ThreadPool() for i in 1:Threads.nthreads() ]
        return new(tpools)
    end
end

# Gets a new array
function getarray(pool::ArrayPool, T::DataType, dims::Int...)
    tid = Threads.threadid()
    tpool = pool.tpools[tid]
    arraylist = get!(tpool.dict, (T, dims), ArrayList(T, dims))

    idx = findfirst(!, arraylist.inuse) # find first non used array
    if idx===nothing # empty arraylist or all arrays are in use
        array = zeros(T, dims)
        push!(arraylist.arrays, array)
        push!(arraylist.inuse, true)
    else
        array = arraylist.arrays[idx]
        arraylist.inuse[idx] = true
    end

    return array

end

# Gets an array of type Float64
function getarray(pool::ArrayPool, dims::Int...)
    return getarray(pool, Float64, dims...)
end

function Base.zeros(pool::ArrayPool, dims::Int...)
    A = getarray(pool, Float64, dims...)
    A .*= 0.0
    return A
end


function free(pool, arrays...)
    tid = Threads.threadid()
    tpool = pool.tpools[tid]
    
    for A in arrays
        id    = objectid(A)
        T     = eltype(A)
        dims  = size(A)
        arraylist = get!(tpool.dict, (T, dims), nothing)
        if arraylist===nothing
            error()
        end

        idx = findfirst(x->id==objectid(x), arraylist.arrays) # find array position
        if idx===nothing
            error()
        end

        arraylist.inuse[idx] = false
    end
end

function freeall(pool)
    for tid in 1:Threads.nthreads()
        tpool = pool[tid]
        for (_, arraylist) in tpool.dict
            arraylist.inuse .= false
        end
    end
end


# pool = ArrayPool()

# A = getarray(pool, Float64, 2,2)
# B = getarray(pool, Float64, 3,3)

# A[1,1] = 55
# free(pool, A)
# C = getarray(pool, Float64, 2,2)
# free(pool, A, B, C)