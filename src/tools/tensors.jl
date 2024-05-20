# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Tensor definitions using Mandel notation

import DataStructures.OrderedDict

const I2 = SVector(1., 1., 1., 0., 0., 0.)
const I4 = SMatrix{6,6}(I)

const Psd = @SArray [
    2/3. -1/3. -1/3. 0. 0. 0.
   -1/3.  2/3. -1/3. 0. 0. 0.
   -1/3. -1/3.  2/3. 0. 0. 0.
      0.    0.    0. 1. 0. 0.
      0.    0.    0. 0. 1. 0.
      0.    0.    0. 0. 0. 1. ]

dev(T::Vec6) = Psd*T # deviatoric tensor

# Tensor invariants
LinearAlgebra.tr(σ::Vec6) = σ[1]+σ[2]+σ[3]

# Deviatoric tensor invariants (Mandel notation)
function J2(σ::Vec6)
    t11, t22, t33, t23, t13, t12 = σ[1], σ[2], σ[3], σ[4]/SR2, σ[5]/SR2, σ[6]/SR2
    return 1/6*( (t11-t22)^2 + (t22-t33)^2 + (t33-t11)^2 ) + t23^2 + t13^2 + t12^2
end


# Third invariant of the deviatoric tensor
function J3(σ::Vec6)
    t11, t22, t33 = σ[1], σ[2], σ[3]
    s23, s13, s12 = σ[4]/SR2, σ[5]/SR2, σ[6]/SR2
    p = t11+t22+t33
    s11 = t11 - 1/3*p
    s22 = t22 - 1/3*p
    s33 = t33 - 1/3*p
    return s11*s22*s33 + 2*s12*s23*s13 - s11*s23^2 - s22*s13^2 - s33*s12^2
end


"""
This function is not precise enouugh...
"""
function eigvals2(T::Vec6)
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
Computes the eigenvalues of a second order tensor written in Mandel notation.
The eigenvalues are sorted from highest to lowest
"""
function eigvals(T::Vec6)
    t11, t22, t33, t23, t13, t12 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2
    
    # full notation
    F = @SArray[ t11  t12  t13
                 t12  t22  t23
                 t13  t23  t33 ]

    L, _ = eigen(F, permute=false, scale=false)

    # put biggest eigenvalue first
    return sort(L, rev=true)
end


"""
Computes eigenvalues and eigenvectors of a second order tensor written in Mandel notation.
The first eigenvalues corresponds to the highest.
The eigenvectors are returned columnwise and disposed in a clockwise coordinate system.
"""
function LinearAlgebra.eigen(T::Vec6)
    t11, t22, t33, t23, t13, t12 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2
    
    # full notation
    F = @SArray[ t11  t12  t13
                 t12  t22  t23
                 t13  t23  t33 ]

    L, V = eigen(F, permute=false, scale=false)

    # put biggest eigenvalue first
    p = sortperm(L, rev=true)
    L = L[p]
    V = V[:, p]
    V = [ V[:,1] V[:,2] normalize(cross(V[:,1], V[:,2])) ]
    # V[:,3] = normalize(cross(V[:,1], V[:,2]))

    return L, V
end


function LinearAlgebra.inv(T::Vec6)
    t11, t22, t33, t23, t13, t12 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2
    
    # full notation
    F = @SArray[ t11  t12  t13
                 t12  t22  t23
                 t13  t23  t33 ]

    G = inv(F)
    return Vec6( G[1,1], G[2,2], G[3,3], SR2*G[2,3], SR2*G[1,3], SR2*G[1,2] )
end


# TODO: do we need this?
LinearAlgebra.norm(T::Vec6) = √dot(T,T)


function rotation_tensor!(V::Array{Float64,2})
    # V : second order tensor with direction cosines (new system axes in old system coordinates)
    # R : fourth order tensor

    lx, ly, lz = V[1,:]
    mx, my, mz = V[2,:]
    nx, ny, nz = V[3,:]

    return @SArray [
            lx*lx      ly*ly      lz*lz    SR2*ly*lz    SR2*lz*lx    SR2*lx*ly
            mx*mx      my*my      mz*mz    SR2*my*mz    SR2*mz*mx    SR2*mx*my
            nx*nx      ny*ny      nz*nz    SR2*ny*nz    SR2*nz*nx    SR2*nx*ny
        SR2*mx*nx  SR2*my*ny  SR2*mz*nz  my*nz+ny*mz  mx*nz+nx*mz  mx*ny+nx*my
        SR2*nx*lx  SR2*ny*ly  SR2*nz*lz  ny*lz+ly*nz  nx*lz+lx*nz  nx*ly+lx*ny
        SR2*lx*mx  SR2*ly*my  SR2*lz*mz  ly*mz+my*lz  lx*mz+mx*lz  lx*my+mx*ly ]
end



"""
Return a dictionary with conventional stress and stress values
from stress and strain tensors defined in Mandel notation.
"""
@inline function stress_strain_dict(σ::Vec6, ε::Vec6, stressmodel::Symbol)
    svm = √(3*J2(σ))

    if stressmodel in (:planestress,:planestrain)
        s1, _, s3 = eigvals(σ)
        return OrderedDict{Symbol,Float64}(
            :sxx => σ[1],
            :syy => σ[2],
            :szz => σ[3],
            :syz => σ[4]/SR2,
            :sxz => σ[5]/SR2,
            :sxy => σ[6]/SR2,
            :svm => svm,
            :s1  => s1,
            :s3  => s3,
            :exx => ε[1],
            :eyy => ε[2],
            :ezz => ε[3],
            :exy => ε[6]/SR2,
        )
    elseif stressmodel==:axisymmetric
        return OrderedDict{Symbol,Float64}(
            :srr => σ[1],
            :syy => σ[2],
            :stt => σ[3],
            :sry => σ[6]/SR2,
            :svm => svm,
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
            :svm => svm,
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