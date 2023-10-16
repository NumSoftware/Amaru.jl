# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
Returns max(x,0)
"""
pos(x) = (abs(x)+x)/2.0

"""
Returns min(x,0)
"""
neg(x) = (-abs(x)+x)/2.0


@inline function norm_rms(V::Vector{Float64})
    return norm(V)/√length(V)
end

# eye function
eye(n::Int64) = Array{Float64}(I,n,n)

function signedmaxabs(x, y)
    if abs(x)>abs(y)
        return x
    else
        return y
    end
end


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

        for i in 1:maxits
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


function findrootinterval(f::Function, x1, Δx)
    x2 = x1 + Δx
    f1 = f(x1)
    f2 = f(x2)
    maxits = 50

    # search for a valid interval
    if f1*f2>0
        for i in 1:maxits
            x2 += Δx*1.6^i
            f2  = f(x2)
            f1*f2<0 && break

            i==maxits && return 0.0, 0.0, failure("findrootinterval: Could not find a valid interval")
        end
    end

    # if !(f1*f2<0)
    #     println()
        
    #     @show f1
    #     @show f2
    # end

    return x1, x2, success()
end

function findroot(f::Function, a, b; tol=(b-a)*0.001, ftol=Inf, method=:default)
    if method==:default
        return findroot_default(f, a, b, tol, ftol)
    elseif method==:bisection
        return findroot_bisection(f, a, b, tol, ftol)
    end
end

function findroot_bisection(f::Function, a, b, tol, ftol)
    fa = f(a)
    fb = f(b)
    # if fa*fb>0
    #     @show a
    #     @show b
    # end

    fa*fb>0 && return 0.0, failure("findroot_bisection: function must have opposite sings at endpoints")

    maxits = 50
    x  = (a+b)/2
    for i in 1:maxits
        fx = f(x)
        if fa*fx<0.0
            b  = x
        else
            a  = x
            fa = fx
        end
        x  = (a+b)/2

        (b-a)/2<tol && abs(fx-fa)<ftol && break
        i==maxits && return 0.0, failure("findroot_bisection: maxits reached")
    end

    isnan(x) && return 0.0, failure("findroot_bisection: NaN result")
    return x, success()
end


function findroot_default(f::Function, a, b, tol, ftol)
    # bracketed method that combines bissection, quadratic interpolation and regula falsi methods
    # for points are used at each itearation
    # x1 and x4 are the endpoints
    # each iteration starts with x1, x2 and x4 points
    # x3 is the newest approximation

    maxits = 50

    x1 = a
    x4 = b
    y1 = f(x1)
    y4 = f(x4)

    y1*y4>0 && return 0.0, failure("findroot: function must have opposite sings at endpoints")

    x2 = (y4*x1 - x4*y1)/(y4-y1) # regula falsi first approximation
    y2 = f(x2)
    x = 0
    local n, err, x2, x3, y2, y3
    n, err = 0, 0.0

    for i in 1:maxits

        # find x3
        if (y2-y1)*(y2-y4)>0 
            # bissection
            x3 = (x1+x4)/2
        else
            # quadratic interpolation
            Y1 = y1/((x1-x2)*(x1-x4))
            Y2 = y2/((x2-x1)*(x2-x4))
            Y4 = y4/((x4-x1)*(x4-x2))
            A  = Y1 + Y2 + Y4
            B  = -(x2+x4)*Y1 - (x1+x4)*Y2 - (x1+x2)*Y4
            C  = x2*x4*Y1 + x1*x4*Y2 + x1*x2*Y4
            D  = B^2 - 4*A*C
            quadratic = D>0

            if quadratic
                R  = D^0.5
                x31 = (-B - R)/(2*A)
                x32 = (-B + R)/(2*A)

                # choose the rigth root
                if x1<x31<x4
                    x3 = x31
                elseif x1<x32<x4
                    x3 = x32
                else
                    quadratic = false
                end
            end

            # regula falsi
            if !quadratic
                x3 = (y4*x1 - x4*y1)/(y4-y1)
            end
        end
        
        y3 = f(x3)

        # order x2 and x3
        if x3<x2
            x2, x3 = x3, x2
            y2, y3 = y3, y2
        end

        err = x3-x2

        # update x1, x2 and x4
        if y1*y2<0 || (y2*y3<0 && abs(y1)<abs(y4))
            # keep x1 and x2 and update x4
            x4 = x3
            y4 = y3
        else
            # update x1 and x2 and keep x4
            x1 = x2
            x2 = x3
            y1 = y2
            y2 = y3
        end

        if err<tol && abs(y3-y2)<ftol
            x = x2 # x3 is discarded
            break
        end

        n=i

        i==maxits && return 0.0, failure("findroot: maxits reached")
    end
    # @show err
    # @show tol
    # @show y2
    # @show y3
    # @show n
    # error()

    return x, success()
end


function brent(f::Function, a, b, tol; maxits::Int=50)
    ftol = 2*eps()
    fa = f(a)
    fb = f(b)

    if abs(fa) < abs(fb)
        # swap bounds
        a, b = b, a
        fa, fb = fb, fa
    end
    c  = a
    fc = fa
    d  = c
    x0 = NaN
    
    bisect = true
    for i in 1:maxits
        abs(b-a) < tol && return b

        # use inverse quadratic interpolation if f(a)!=f(b)!=f(c)
        if abs(fa-fc)>ftol && abs(fb-fc)>ftol
            x = a*fb*fc/((fa-fb)*(fa-fc)) +
                b*fa*fc/((fb-fa)*(fb-fc)) +
                c*fa*fb/((fc-fa)*(fc-fb))
        else
            x = b - fb * (b-a)/(fb-fa) # regula falsi
        end

        # Use bisection method if satisfies the conditions.
        delta = 2*eps()*abs(b)
        Dxb = abs(x-b)
        Dcb = abs(b-c)
        Dcd = abs(c-d)
        if x<(3*a+b)/4 && x>b || bisect && Dxb>=Dcb/2 || !bisect && Dxb>=Dcd/2 || bisect && Dcb<delta || !bisect && Dcd<delta
            x = (a+b)/2
            bisect = true
        else
            bisect = false
        end

        abs(x-x0)<tol && return x
        x0 = x

        fx = f(x)
        abs(fx) < ftol && return x

        d = c
        c = b

        if fa*fx<0
            b  = x
            fb = fx
        else
            a  = x
            fa = fx
        end

        if abs(fa) < abs(fb)
            a, b = b, a
            fa, fb = fb, fa
        end
    end
    error("Max iteration exceeded")
end


function spline_coefficients(X,Y)
    n = length(X)
    neq = 4*(n-1) # number of equations
    # System M A = B
    M = zeros(neq, neq)
    B = zeros(neq)
    k = 0 # equation number

    # condition 1: matching y values
    for i in 1:n-1
        k += 1
        M[k, 4*i-3] = 1
        M[k, 4*i-2] = X[i]
        M[k, 4*i-1] = X[i]^2
        M[k, 4*i-0] = X[i]^3
        B[k] = Y[i]
        k +=1
        M[k, 4*i-3] = 1
        M[k, 4*i-2] = X[i+1]
        M[k, 4*i-1] = X[i+1]^2
        M[k, 4*i-0] = X[i+1]^3
        B[k] = Y[i+1]
    end

    # condition 2: same internal first derivatives
    for i=2:n-1
        k += 1
        M[k, 4*(i-1)-2] = 1
        M[k, 4*(i-1)-1] = 2*X[i]
        M[k, 4*(i-1)-0] = 3*X[i]^2
        M[k, 4*i-2] = -1
        M[k, 4*i-1] = -2*X[i]
        M[k, 4*i-0] = -3*X[i]^2
    end

    # condition 3: same internal second derivatives
    for i=2:n-1
        k += 1
        M[k, 4*(i-1)-1] = 2
        M[k, 4*(i-1)-0] = 6*X[i]
        M[k, 4*i-1] = -2
        M[k, 4*i-0] = -6*X[i]
    end

    # condition 4: natural spline
    k += 1
    M[k, 3] = 2
    M[k, 4] = 6*X[1]
    k += 1
    M[k, end-1] = 2
    M[k, end-0] = 6*X[end]
    
    return inv(M)*B
end


# Spline
struct Spline
    X::Array{Float64,1}  # x values
    Y::Array{Float64,1}  # y values
    C::Array{Float64,1}
    function Spline(X, Y)
        C = spline_coefficients(X, Y)
        return new(X, Y, C)
    end
end

# Overloading the () operator
function (spline::Spline)(x::Real)
    n = length(spline.X)
    for i in 1:n-1
        if x <= spline.X[i+1]
            a, b, c, d = spline.C[4*i-3 : 4*i]
            return a + b*x + c*x^2 + d*x^3
        end
    end
end


# # Parametric spline
# struct Spline
#     X::Array{Float64,1}  # x values
#     Y::Array{Float64,1}  # y values
#     # S::Array{Float64,1}  # s values
#     D::Array{Float64,1}  # parametric spline coefficients
#     function Spline(X, Y)
#         # a, b = X[1], X[end]
#         # S = [ (x-a)/(b-a) for x in X ]
#         D = spline_coefficients(Y)
#         return new(X, Y, D)
#         # return new(X, Y, S, D)
#     end
# end


# # Overloading () operator to be used with the Spline structure
# function (spline::Spline)(x::Real)
#     X = spline.X
#     Y = spline.Y
#     # S = spline.S
#     D = spline.D

#     a, b = X[1], X[end]
#     @assert a <= x <= b

#     # s = (x-a)/(b-a)
    
#     # Find polynomial function
#     i = max(1, findfirst(xi->x<=xi, X) - 1)

#     # Parametric coordinate
#     # t = (s-S[i])/(S[i+1]-S[i])
#     t = (x-X[i])/(X[i+1]-X[i])

    
#     # Spline coefficients
#     a = Y[i]
#     b = D[i]
#     c = 3*(Y[i+1]-Y[i]) - 2*D[i] - D[i+1]
#     d = 2*(Y[i]-Y[i+1]) + D[i] + D[i+1]
    
#     return a + b*t + c*t^2 + d*t^3
# end


function deriv(spline::Spline, x::Real)
    X = spline.X
    Y = spline.Y
    # S = spline.S
    D = spline.D

    a, b = X[1], X[end]
    @assert a <= x <= b

    s = (x-a)/(b-a)

    # Find polynomial function
    # i = max(1, findfirst(si->s<=si, S) - 1)
    i = max(1, findfirst(xi->x<=xi, X) - 1)

    # Parametric coordinate
    # t = (s-S[i])/(S[i+1]-S[i])
    t = (x-X[i])/(X[i+1]-X[i])
    
    # Spline coefficients
    a = Y[i]
    b = D[i]
    c = 3*(Y[i+1]-Y[i]) - 2*D[i] - D[i+1]
    d = 2*(Y[i]-Y[i+1]) + D[i] + D[i+1]

    dydt = b + 2*c*t + 3*d*t^2
    dtdx = 1/(X[i+1]-X[i])
    
    return dydt*dtdx
end