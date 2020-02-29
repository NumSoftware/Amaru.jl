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
