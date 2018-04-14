# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Types
const Vect = Array{Float64, 1}
const Matx = Array{Float64, 2}

export @showm

# Fancy matrix printing
function print_matrix(M::Array{Float64,2})
    n, m = size(M)
    for i=1:n
        for j=1:m
            @printf( "%23.11e", M[i,j] )
        end
        println()
    end
end

# Pseudo determinant of non-square matrices
function norm2(J)

    if ndims(J)==1; return norm(J) end

    r, c = size(J)
    if r==1; return norm(J) end
    if r==2 && c==3
        j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        j2 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
        j3 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
        return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian determinant
    end
    if r==c; return det(J) end
    error("No rule to calculate norm2 of a $r x $c matrix")
end

macro showm(M)
    return quote
        println($(string(M)), "=")
        Base.showarray(STDOUT, $(esc(M)), false)
        println()
        #=
        dim = size($M)
        if length(dim)==1
            print(" [ ")
            for i=1:dim[1]
                if abs($M[i])<1e-3
                    @printf("%13.6e ", $M[i])
                else
                    @printf("%13.6f ", $M[i])
                end
            end
            println(" ]")
        else
            print(" [ ")
            for i=1:dim[1]
                for j=1:dim[2]
                    if abs($M[i,j])<1e-3
                        @printf("%13.6e ", $M[i,j])
                    else
                        @printf("%13.6f ", $M[i,j])
                    end
                end
                if i<dim[1]
                    print("\n   ")
                end
            end
            println(" ]")
        end
        =#
    end
end


# C  = α*A*B
# C  = α*A'*B'
# C += α*A'*B'
# C -= α*A'*B'
macro gemm(expr)
    β = 0.0
    s = 1.0
    if expr.head == :(=)
    elseif expr.head == :(+=)
        β = 1.0
    elseif expr.head == :(-=)
        s = -1.0
        β =  1.0
    else
        error("@gemm: =, +=, -= operator expected, found $(expr.head)")
    end

    C = expr.args[1]
    rhs = expr.args[2]

    if rhs.args[1] !=  :(*)
        error("@gemm: * operator expected, found $(rhs.args[1])")
    end

    if length(rhs.args) == 4
        α = rhs.args[2]
        A = rhs.args[3]
        B = rhs.args[4]
    else
        α = 1.0
        A = rhs.args[2]
        B = rhs.args[3]
    end
    
    tA = 'N'
    if typeof(A) == Expr 
        if A.head == Symbol("'"); 
            tA = 'T' 
            A  = A.args[1]
        end
    end

    tB = 'N'
    if typeof(B) == Expr
        if B.head == Symbol("'"); 
            tB = 'T' 
            B  = B.args[1]
        end
    end

    return :( BLAS.gemm!($tA, $tB, $(esc(α))*$s, $(esc(A)), $(esc(B)), $β, $(esc(C)) ) )
end


# Y  = α*A*X
# Y  = α*A'*X 
# Y += α*A'*X 
# Y -= α*A'*X 
macro gemv(expr)
    β = 0.0
    s = 1.0
    if expr.head == :(=)
    elseif expr.head == :(+=)
        β = 1.0
    elseif expr.head == :(-=)
        s = -1.0
        β =  1.0
    else
        error("@gemv: =, +=, -= operator expected, found $(expr.head)")
    end

    Y = expr.args[1]
    rhs = expr.args[2]

    if rhs.args[1] !=  :(*)
        error("@gemv: * operator expected, found $(rhs.args[1])")
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
    
    tA = 'N'
    if typeof(A) == Expr 
        if A.head == Symbol("'"); 
            tA = 'T' 
            A  = A.args[1]
        end
    end

    return :( BLAS.gemv!($tA, $(esc(α))*$s, $(esc(A)), $(esc(X)), $(β), $(esc(Y)) ) )
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
            $(esc(Y))[:] = 0.0
            BLAS.axpy!( $(esc(α))*$s, $(esc(X)), $(esc(Y)) )
        end
    else
        return :( BLAS.axpy!( $(esc(α))*$s, $(esc(X)), $(esc(Y)) ) )
    end
end


