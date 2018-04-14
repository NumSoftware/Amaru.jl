# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Tensor definitions using Mandel notation

export Tensor2, Tensor4

const Tensor2 = Array{Float64,1}
const Tensor4 = Array{Float64,2}

tensor2() = zeros(6)
tensor4() = zeros(6,6)

const tI  = [1., 1., 1., 0., 0., 0.]

import Base.trace
J1(T::Tensor2) = sum(T[1:3])
const trace = J1

@inline function J2(T::Tensor2)
    t11, t22, t33, t12, t23, t13 = T
    return 0.5*t11*t11 + 0.5*t22*t22 + 0.5*t33*t33 + t12*t12 + t23*t23 + t13*t13
end


const Psd = [  
    2/3. -1/3. -1/3. 0. 0. 0.
   -1/3.  2/3. -1/3. 0. 0. 0.
   -1/3. -1/3.  2/3. 0. 0. 0.
      0.    0.    0. 1. 0. 0.
      0.    0.    0. 0. 1. 0.
      0.    0.    0. 0. 0. 1. ]

const Isym = eye(6)

@inline function dev(T::Tensor2) # deviatoric tensor
    return Psd*T
end

function J2D(T::Tensor2)
    #return J2(Psd*T)
    t11, t22, t33, t12, t23, t13 = T
    return 1/6*( (t11-t22)^2 + (t22-t33)^2 + (t33-t11)^2 ) + t12*t12 + t23*t23 + t13*t13
end


function J3D(T::Tensor2)
    return J3(Psd*T)
end


function principal(T::Tensor2)::Vect
    t11, t22, t33, t12, t23, t13 = T
    i1 = t11 + t22 + t33
    i2 = t11*t22 + t22*t33 + t11*t33 - t12*t12 - t23*t23 - t13*t13

    if i1==0.0 && i2==0.0
        return zeros(3)
    end

    i3 = t11*(t22*t33 - t23*t23) - t12*(t12*t33 - t23*t13) + t13*(t12*t23 - t22*t13)
    val = round( (2*i1^3 - 9*i1*i2 + 27*i3 )/( 2*(i1^2 - 3*i2)^(3/2) ), 14 )
    θ = 1/3*acos( val )

    r = 2/3*√(i1^2-3*i2)

    s1 = i1/3 + r*cos(θ)
    s2 = i1/3 + r*cos(θ - 2*π/3)
    s3 = i1/3 + r*cos(θ - 4*π/3)

    # sorting
    if s1<s2 s1,s2 = s2,s1 end
    if s2<s3 s2,s3 = s3,s2 end
    if s1<s2 s1,s2 = s2,s1 end

    P = [ s1, s2, s3 ]

    return P
end

function principal_dir(T::Tensor2)
    # full notation
    F = [ T[1]  T[4]  T[6] ;
          T[4]  T[2]  T[5] ;
          T[6]  T[5]  T[3] ]
    L, V = eig(F, permute=false, scale=false)

    # force a clockwise system
    if norm( cross(V[:,2], V[:,3]) - V[:,1] ) > 1e-5
        L[2], L[3] = L[3], L[2]
        V[:,2], V[:,3] = V[:,3], V[:,2]
    end

    # find max value
    val, idx = findmax(L)
    shift = 1 - idx
    L = circshift(L, shift)
    V = circshift(V, (0,shift))
    return L, V
end


const V2M = [ 1., 1., 1., √2., √2., √2. ]
const M2V = [ 1., 1., 1., √.5, √.5, √.5 ] # Use .* operator

@inline function tfull(T::Tensor2)
    t1, t2, t3, t4, t5, t6 = T
    return [ 
        t1 t4 t6
        t4 t2 t5
        t6 t5 t3 ]
end

function matrix2Mandel(M::Array{Float64,2})
    t11 = M[1,1]
    t22 = M[2,2]
    t33 = M[3,3]
    t12 = M[1,2]
    t23 = M[2,3]
    t13 = M[1,3]
    return [t11, t22, t33, t12, t23, t13] .* V2M
end

import Base.norm
function norm(T::Tensor2)
    t11, t22, t33, t12, t23, t13 = T
    return (t11*t11 + t22*t22 + t33*t33 + 2.0*t12*t12 + 2.0*t23*t23 + 2.0*t13*t13)^0.5
end

import Base.dot 
function dot(T1::Tensor2, T2::Tensor2)
    t11, t22, t33, t12, t23, t13 = T1
    s11, s22, s33, s12, s23, s13 = T2
    return t11*s11 + t22*s22 + t33*s33 + 2.0*t12*s12 + 2.0*t23*s23 + 2.0*t13*s13
end

function dyad(T1::Tensor2, T2::Tensor2)
    # TODO: convert to mandel first and fix later
    return T1 * T2'
end

⊗ = dyad

function inner(T1::Tensor4, T2::Tensor4)
    return sum(T1 .* T2)
end

function inner(T1::Tensor4, T2::Tensor2)
    # TODO: convert to mandel first and fix later
    return T1 * T2
end

function inner(T1::Tensor2, T2::Tensor4)
    # TODO: convert to mandel first and fix later
    return T2*T1
end

function inner(T1::Tensor2, T2::Tensor4, T3::Tensor2)
    # TODO: convert to mandel first and fix later
    return dot(T2*T1, T3)
end

∷ = inner

function tensor_rot!(V::Tensor2, T::Array{Float64,2})
    l1, m1, n1 = V[:,1]
    l2, m2, n2 = V[:,2]
    l3, m3, n3 = V[:,3]

    T[1,1] =     l1*l1;  T[1,2] =     m1*m1;  T[1,3] =     n1*n1;  T[1,4] =       l1*m1;  T[1,5] =       m1*n1;  T[1,6] =       n1*l1     
    T[2,1] =     l2*l2;  T[2,2] =     m2*m2;  T[2,3] =     n2*n2;  T[2,4] =       l2*m2;  T[2,5] =       m2*n2;  T[2,6] =       n2*l2     
    T[3,1] =     l3*l3;  T[3,2] =     m3*m3;  T[3,3] =     n3*n3;  T[3,4] =       l3*m3;  T[3,5] =       m3*n3;  T[3,6] =       n3*l3     
    T[4,1] =     l1*l2;  T[4,2] =     m1*m2;  T[4,3] =     n1*n2;  T[4,4] = l1*m2+l2*m1;  T[4,5] = m1*n2+m2*n1;  T[4,6] = l1*n2+l2*n1     
    T[5,1] =     l2*l3;  T[5,2] =     m2*m3;  T[5,3] =     n2*n3;  T[5,4] = l2*m3+l3*m2;  T[5,5] = m2*n3+m3*n2;  T[5,6] = l2*n3+l3*n2     
    T[6,1] =     l3*l1;  T[6,2] =     m3*m1;  T[6,3] =     n3*n1;  T[6,4] = l3*m1+l1*m3;  T[6,5] = m3*n1+m1*n3;  T[6,6] = l3*n1+l1*n3 
end

function tensor_rot(V::Tensor2)
    T = zeros(6,6)
    transformer!(V, T)
    return T
end
