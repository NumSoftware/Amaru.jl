# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Types
const Vect = AbstractArray{Float64, 1}
const Matx = AbstractArray{Float64, 2}

# Static arrays
const Vec2 = SVector{2, Float64}
const Vec3 = SVector{3, Float64}
const Vec4 = SVector{4, Float64}
const Vec6 = SVector{6, Float64}
const Mat6x6 = SMatrix{6, 6, Float64, 36} 


function Base.convert(::Type{Vec3}, A::Array{Float64,1})
    n = length(A)
    @assert n in (1,2,3)
    if n==3
        return SVector(A[1], A[2], A[3])
    elseif n==2
        return SVector(A[1], A[2], 0.0)
    else
        return SVector(A[1], 0.0, 0.0)
    end
end



eye(n) = Array{Float64,2}(I, n, n)

# Fancy matrix printing
function print_matrix(M::Array{Float64,2})
    n, m = size(M)
    for i in 1:n
        for j in 1:m
            @printf( "%23.11e", M[i,j] )
        end
        println()
    end
end

function extend!(V::AbstractArray{Float64,1}, n::Int)
    while length(V)<n
        push!(V, 0.0)
    end
    return V
end

# Pseudo determinant of non-square matrices
function norm2(J)

    if ndims(J)==1; return norm(J) end

    r, c = size(J)
    r==c && return det(J)
    (r==1 || c==1) && return norm(J)

    if r==2 && c==3
        j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        j2 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
        j3 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
        return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian norm
    end
    if r==3 && c==2
        j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        j2 = J[1,1]*J[3,2] - J[1,2]*J[3,1]
        j3 = J[2,1]*J[3,2] - J[2,2]*J[3,1]
        return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian norm
    end

    error("No rule to calculate norm2 of a $r x $c matrix")
end

# Y  = α*A*X
# Y  = α*A'*X
# Y += α*A'*X
# Y -= α*A'*X
macro mul(expr)
    β = 0.0
    s = 1.0
    if expr.head == :(=)
    elseif expr.head == :(+=)
        β = 1.0
    elseif expr.head == :(-=)
        s = -1.0
        β =  1.0
    else
        error("@mul: =, +=, -= operator expected, found $(expr.head)")
    end

    Y = expr.args[1]
    rhs = expr.args[2]

    if rhs.args[1] !=  :(*)
        error("@mul: * operator expected, found $(rhs.args[1])")
    end

    if length(rhs.args) == 4
        α = rhs.args[2]
        A = rhs.args[3]
        X = rhs.args[4]
    else
        α = 1.0
        A = rhs.args[2]
        X = rhs.args[3]
    end

    return esc(Expr(:call, :mul!, Y, A, X, :($α*$s), β))
end


# Y  = α*X
# Y += α*X
# Y -= α*X
macro axpy(expr)
    α = 1.0
    β = 0.0
    s = 1.0
    if expr.head == :(=)
    elseif expr.head == :(+=)
        β = 1.0
    elseif expr.head == :(-=)
        s = -1.0
        β =  1.0
    else
        error("@axpy: =, +=, -= operator expected, found $(expr.head)")
    end

    Y   = expr.args[1]
    rhs = expr.args[2]

    if typeof(rhs)==Expr
        if rhs.args[1] !=  :(*)
            error("@axpy: * operator expected, found $(rhs.args[1])")
        end

        if length(rhs.args) == 3
            α = rhs.args[2]
            X = rhs.args[3]
        else
            X = rhs.args[2]
        end
    else
        X = rhs
    end

    if β == 0.0
        return quote
            $(esc(Y)) .= 0.0
            BLAS.axpy!( $(esc(α))*$s, $(esc(X)), $(esc(Y)) )
        end
    else
        return :( BLAS.axpy!( $(esc(α))*$s, $(esc(X)), $(esc(Y)) ) )
    end
end


function Base.split(A::AbstractArray, knots::AbstractArray)
    # @s A
    # @s knots

    # R = Array{eltype(A),1}[]
    R = Array{Int,1}[]
    n = length(A)

    knots = unique(knots)

    if knots[end] == A[end]
        knots = knots[1:end-1]
    end
    
    if knots[1] != A[1]
        knots = [A[1]; knots]
    end

    idx = 0
    for i in 1:n
        for j in idx+1:length(knots)
            if A[i]>=knots[j]
                push!(R, Int[])
                idx += 1
            end
        end
        # push!(R[idx], A[i])
        push!(R[idx], i)
    end
    return R
end

# dumps the content of a sparse matrix to a file in coordinate format
function dump_matrix(A::SparseMatrixCSC, filename::String)
    n, m = size(A)
    rows, cols, vals = findnz(A)
    
    open(filename, "w") do file
        println(file, "%%MatrixMarket matrix coordinate real general")
        println(file, "$n $m $(length(vals))")
        for i in 1:length(vals)
            @printf(file, "%4d %4d %f\n", rows[i], cols[i], vals[i])
        end
    end
end