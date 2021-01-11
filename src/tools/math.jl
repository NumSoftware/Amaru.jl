# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
Returns max(x,0)
"""
pos(x) = (abs(x)+x)/2.0

"""
Returns min(x,0)
"""
neg(x) = (-abs(x)+x)/2.0


"""
Returns a vector with the real roots of the quadratic polynomial ax^2 + bx + c
"""
function square_roots(a,b,c)
    D = b^2 - 4*a*c # discriminant

    if D<=0.0 # no real roots
        return Float64[]
    else
        x1 = (-b+sqrt(D))/(2.0*a)
        x2 = (-b-sqrt(D))/(2.0*a)
        X  = Float64[ x1, x2 ]
        return sort(X)
    end
end


"""
Returns a vector with the real roots of the cubic polynomial ax^3 + bx^2 + cx + d
"""
function cubic_roots(a,b,c,d)
    A = b/a
    B = c/a
    C = d/a

    Q = (3*B-A^2)/9
    R = (9*A*B-27*C-2*A^3)/54
    D = Q^3 + R^2 # discriminant

    ftol = 1e-4

    if D<=0 # 3 real roots
        th = acos(R/sqrt(-Q^3))
        x1 = 2*sqrt(-Q)*cos(th/3) - A/3
        x2 = 2*sqrt(-Q)*cos(th/3 + 2*π/3) - A/3
        x3 = 2*sqrt(-Q)*cos(th/3 + 4*π/3) - A/3
        X = Float64[x1, x2, x3]

        F = a*X.^3 .+ b.*X.^2 .+ c.*X .+ d
        f = maximum(abs, F)
        f > ftol && @warn "cubic_roots: residue ($f) greather than ftol"
        return sort(X)
    else # one real root: use Newton method
        maxits = 40
        tol = 1e-12
        x   = 0.0
        x0  = 0.0
        der = 3*a*x0^2 + 2*b*x0 + c
        der==0.0 && (x0=π)

        for i=1:maxits
            f   = a*x0^3 + b*x0^2 + c*x0 + d
            der = 3*a*x0^2 + 2*b*x0 + c
            x   = x0 - f/der
            err = abs(x - x0)
            err<tol && break
            x0 = x
            i == maxits && @warn "cubic_roots: max number of iterations reached."
        end
        f = a*x^3 + b*x^2 + c*x + d
        f > ftol && @warn "cubic_roots: residue ($f) greather than ftol"
        return Float64[ x ]
    end
end



function spline_coefficients(Y::Array{Float64,1})
    n   = length(Y)
    neq = n
    M   = zeros(neq, neq)
    B   = zeros(neq)
    
    M[1,1:2]         = [2, 1]
    M[end,end-1:end] = [1, 2]
    B[1]             = 3*(Y[2]-Y[1])
    B[end]           = 3*(Y[n]-Y[n-1])

    v141 = [1, 4, 1]
    for i=2:n-1
        M[i, i-1:i+1] .= v141
        B[i]           = 3*(Y[i+1]-Y[i-1])
    end

    return M\B
end


# Parametric spline
struct Spline
    X::Array{Float64,1}  # y values
    Y::Array{Float64,1}  # y values
    D::Array{Float64,1}  # parametric spline coefficients
    function Spline(X, Y)
        D = spline_coefficients(Y)
        return new(X, Y, D)
    end
end


# Overloading () operator to be used with the ParamSpline structure
function (spline::Spline)(x::Real)
    X = spline.X
    Y = spline.Y
    D = spline.D

    @assert X[1] <= x <= X[end]
    
    # Find polynomial function
    n = length(Y)
    i = max(1, findfirst(xi->x<=xi, X) - 1)

    # Parametric coordinate
    t = (x-X[i])/(X[i+1]-X[i])
    
    # Spline coefficients
    a = Y[i]
    b = D[i]
    c = 3*(Y[i+1]-Y[i]) - 2*D[i] - D[i+1]
    d = 2*(Y[i]-Y[i+1]) + D[i] + D[i+1]

    return a + b*t + c*t^2 + d*t^3
end