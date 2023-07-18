# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Tensor definitions using Mandel notation

import DataStructures.OrderedDict
export Vec6, Tensor4

# const Tensor2 = Array{Float64,1}
# const Tensor4 = Array{Float64,2}



const V2M = [ 1., 1., 1., SR2, SR2, SR2 ]
const M2V = [ 1., 1., 1., 1.0/SR2, 1.0/SR2, 1.0/SR2 ] # Use .* operator

const tI  = @Vec6 [1., 1., 1., 0., 0., 0.]
const Isym = eye(6)

const Psd = [
    2/3. -1/3. -1/3. 0. 0. 0.
   -1/3.  2/3. -1/3. 0. 0. 0.
   -1/3. -1/3.  2/3. 0. 0. 0.
      0.    0.    0. 1. 0. 0.
      0.    0.    0. 0. 1. 0.
      0.    0.    0. 0. 0. 1. ]

dev(T::Vec6) = Psd*T # deviatoric tensor

# Tensor invariants
LinearAlgebra.tr(T::Vec6) = sum(T[1:3])
J1(T::Vec6) = sum(T[1:3])
J2(T::Vec6) = 0.5*dot(T,T)


# Deviatoric tensor invariants
function J2D(T::Vec6)
    #return J2(Psd*T)
    t11, t22, t33, t12, t23, t13 = T
    t12 /= SR2
    t23 /= SR2
    t13 /= SR2
    return 1/6*( (t11-t22)^2 + (t22-t33)^2 + (t33-t11)^2 ) + t12*t12 + t23*t23 + t13*t13
end


function J3D(T::Vec6)
    return J3(Psd*T)
end


"""
Computes the eigenvalues of a second order tensor written in Mandel notation.
The eigenvalues are sorted from highest to lowest
"""
function eigvals(T::Vec6)
    @assert length(T) == 6

    # full notation
    F = [ T[1]      T[6]/SR2  T[5]/SR2 ;
          T[6]/SR2  T[2]      T[4]/SR2 ;
          T[5]/SR2  T[4]/SR2  T[3]     ]
    L, _ = eigen!(F, permute=false, scale=false)

    # put biggest eigenvalue first
    return sort!(L, rev=true)
end

"""
This function is not precise enouugh...
"""
function eigvals2(T::Vec6)::Vect
    @assert length(T) == 6

    t11, t22, t33, t12, t23, t13 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2

    i1 = t11 + t22 + t33
    i2 = t11*t22 + t22*t33 + t11*t33 - t12*t12 - t23*t23 - t13*t13

    i1==0.0 && i2==0.0 && return zeros(3)

    i3  = t11*(t22*t33 - t23*t23) - t12*(t12*t33 - t23*t13) + t13*(t12*t23 - t22*t13)
    val = (2*i1^3 - 9*i1*i2 + 27*i3 )/( 2*(i1^2 - 3*i2)^(3/2) )
    val = clamp(val, -1.0, 1.0) # to avoid 1.000000000000001

    θ = 1/3*acos( val )

    r = 2/3*√(i1^2-3*i2)

    s1 = i1/3 + r*cos(θ)
    s2 = i1/3 + r*cos(θ - 2*π/3)
    s3 = i1/3 + r*cos(θ - 4*π/3)

    # sorting
    if s1<s2; s1,s2 = s2,s1 end
    if s2<s3; s2,s3 = s3,s2 end
    if s1<s2; s1,s2 = s2,s1 end

    return Vec3(s1, s2, s3)
end


"""
Computes eigenvalues and eigenvectors of a second order tensor written in Mandel notation.
The first eigenvalues corresponds to the highest.
The eigenvectors are returned columnwise and disposed in a clockwise coordinate system.
"""
function LinearAlgebra.eigen(T::Vec6)
    @assert length(T) == 6

    # full notation
    T4oR2 = T[4]/SR2
    T5oR2 = T[5]/SR2
    T6oR2 = T[6]/SR2
    F = [ T[1]   T6oR2  T5oR2
          T6oR2  T[2]   T4oR2
          T5oR2  T4oR2  T[3]     ]
    L, V = eigen!(F, permute=false, scale=false)

    # put biggest eigenvalue first
    p = sortperm(L, rev=true)
    L = L[p]
    V = V[:,p]
    V[:,3] .= normalize(cross(V[:,1], V[:,2]))

    return L, V
end


function tfull(T::Vec6)
    t1, t2, t3, t4, t5, t6 = T.*M2V
    return [
        t1 t6 t5
        t6 t2 t4
        t5 t4 t3 ]
end


function matrix2Mandel(M::Array{Float64,2})
    t11 = M[1,1]
    t22 = M[2,2]
    t33 = M[3,3]
    t23 = M[2,3]*SR2
    t13 = M[1,3]*SR2
    t12 = M[1,2]*SR2
    return [t11, t22, t33, t23, t13, t12]
end


LinearAlgebra.norm(T::Vec6) = √dot(T,T)


# function dyad(T1::Vec6, T2::Vec6)
#     return T1 * T2'
# end


# ⊗ = dyad


# function inner(T1::Tensor4, T2::Tensor4)
#     return sum(T1 .* T2)
# end


# function inner(T1::Tensor4, T2::Vec6)
#     return T1 * T2
# end


# function inner(T1::Vec6, T2::Tensor4)
#     return T2*T1
# end


# function inner(T1::Vec6, T2::Tensor4, T3::Vec6)
#     return dot(T2*T1, T3)
# end

# ∷ = inner


# function set_tensor_rot!(V::Array{Float64,2}, R::Tensor4)
#     # V : second order tensor with direction cosines (new system axes in old system coordinates)
#     # R : fourth order tensor

#     lx, ly, lz = V[1,:]
#     mx, my, mz = V[2,:]
#     nx, ny, nz = V[3,:]

#     R[1,1] =     lx*lx;  R[1,2] =     ly*ly;  R[1,3] =     lz*lz;   R[1,4] =   SR2*ly*lz;  R[1,5] =   SR2*lz*lx;  R[1,6] =   SR2*lx*ly;
#     R[2,1] =     mx*mx;  R[2,2] =     my*my;  R[2,3] =     mz*mz;   R[2,4] =   SR2*my*mz;  R[2,5] =   SR2*mz*mx;  R[2,6] =   SR2*mx*my;
#     R[3,1] =     nx*nx;  R[3,2] =     ny*ny;  R[3,3] =     nz*nz;   R[3,4] =   SR2*ny*nz;  R[3,5] =   SR2*nz*nx;  R[3,6] =   SR2*nx*ny;
#     R[4,1] = SR2*mx*nx;  R[4,2] = SR2*my*ny;  R[4,3] = SR2*mz*nz;   R[4,4] = my*nz+ny*mz;  R[4,5] = mx*nz+nx*mz;  R[4,6] = mx*ny+nx*my;
#     R[5,1] = SR2*nx*lx;  R[5,2] = SR2*ny*ly;  R[5,3] = SR2*nz*lz;   R[5,4] = ny*lz+ly*nz;  R[5,5] = nx*lz+lx*nz;  R[5,6] = nx*ly+lx*ny;
#     R[6,1] = SR2*lx*mx;  R[6,2] = SR2*ly*my;  R[6,3] = SR2*lz*mz;   R[6,4] = ly*mz+my*lz;  R[6,5] = lx*mz+mx*lz;  R[6,6] = lx*my+mx*ly;
#     return R
# end


function rotation_tensor!(V::Array{Float64,2})
    # V : second order tensor with direction cosines (new system axes in old system coordinates)
    # R : fourth order tensor

    lx, ly, lz = V[1,:]
    mx, my, mz = V[2,:]
    nx, ny, nz = V[3,:]

    return @Mat6x6 [
            lx*lx      ly*ly      lz*lz    SR2*ly*lz    SR2*lz*lx    SR2*lx*ly
            mx*mx      my*my      mz*mz    SR2*my*mz    SR2*mz*mx    SR2*mx*my
            nx*nx      ny*ny      nz*nz    SR2*ny*nz    SR2*nz*nx    SR2*nx*ny
        SR2*mx*nx  SR2*my*ny  SR2*mz*nz  my*nz+ny*mz  mx*nz+nx*mz  mx*ny+nx*my
        SR2*nx*lx  SR2*ny*ly  SR2*nz*lz  ny*lz+ly*nz  nx*lz+lx*nz  nx*ly+lx*ny
        SR2*lx*mx  SR2*ly*my  SR2*lz*mz  ly*mz+my*lz  lx*mz+mx*lz  lx*my+mx*ly ]
end



# function tensor_rot(V::Array{Float64,2})
#     T = zeros(6,6)
#     tensor_rot!(V, T)
#     return T
# end


"""
Return a dictionary with conventional stress and stress values
from stress and strain tensors defined in Mandel notation.
"""
@inline function stress_strain_dict(σ::Vec6, ε::Vec6, stressmodel::String)
    
    if stressmodel in ("plane-stress","plane-strain")
        s1, _, s3 = eigvals(σ)
        return OrderedDict{Symbol,Float64}(
            :sxx => σ[1],
            :syy => σ[2],
            :szz => σ[3],
            :sxy => σ[6]/SR2,
            :s1  => s1,
            :s3  => s3,
            :exx => ε[1],
            :eyy => ε[2],
            :ezz => ε[3],
            :exy => ε[6]/SR2,
        )
    elseif stressmodel=="axisymmetric"
        return OrderedDict{Symbol,Float64}(
            :srr => σ[1],
            :syy => σ[2],
            :stt => σ[3],
            :sry => σ[6]/SR2,
            :err => ε[1],
            :eyy => ε[2],
            :ett => ε[3],
            :ery => ε[6]/SR2,
        )
    else
        s1, s2, s3 = eigvals(σ)
        return OrderedDict{Symbol,Float64}(
            :sxx => σ[1],
            :syy => σ[2],
            :szz => σ[3],
            :syz => σ[4]/SR2,
            :sxz => σ[5]/SR2,
            :sxy => σ[6]/SR2,
            :s1  => s1,
            :s2  => s2,
            :s3  => s3,
            :exx => ε[1],
            :eyy => ε[2],
            :ezz => ε[3],
            :eyz => ε[4]/SR2,
            :exz => ε[5]/SR2,
            :exy => ε[6]/SR2,
        )
    end

end