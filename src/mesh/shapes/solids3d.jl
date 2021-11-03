#    Amaru - Finite Element Library


# TET4 shape
# ==========


# natural coordinates
const coords_TET4 =
[  0.0  0.0  0.0
   1.0  0.0  0.0
   0.0  1.0  0.0
   0.0  0.0  1.0 ]

const facet_idxs_TET4 =
    [ [1, 4, 3],     [1, 2, 4],      [1, 3, 2],     [2, 3, 4]  ]

const edge_idxs_TET4 =
    [ [1, 2],    [2, 3],    [3, 1],     [1, 4],      [2, 4],      [3, 4] ]

function shape_func_TET4(R::AbstractArray{<:Float64,1})
    r, s, t = R

    N = Array{Float64}(undef,4)
    N[1] = 1.0-r-s-t
    N[2] = r
    N[3] = s
    N[4] = t
    return N
end

function shape_deriv_TET4(R::AbstractArray{<:Float64,1})
    r, s, t = R

    D = Array{Float64}(undef,4,3)
    D[1,1] = -1.0;   D[1,2]= -1.0;   D[1,3]= -1.0
    D[2,1] =  1.0;   D[2,2]=  0.0;   D[2,3]=  0.0
    D[3,1] =  0.0;   D[3,2]=  1.0;   D[3,3]=  0.0
    D[4,1] =  0.0;   D[4,2]=  0.0;   D[4,3]=  1.0
    return D
end

# constructor
function MakeTET4()
    shape             = CellShape()
    shape.name        = "TET4"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 4
    shape.basic_shape = shape
    shape.vtk_type    = VTK_TETRA
    shape.facet_idxs  = facet_idxs_TET4
    shape.edge_idxs   = edge_idxs_TET4
    shape.facet_shape = TRI3
    shape.nat_coords  = coords_TET4
    shape.quadrature  = Dict( 0 => TET_IP4,  1 => TET_IP1,  4 => TET_IP4,   5 => TET_IP5,  11 => TET_IP11 )
    shape.func        = shape_func_TET4
    shape.deriv       = shape_deriv_TET4
    return shape
end


# Registration
const  TET4 = MakeTET4()
export TET4



# TET10 shape
# ===========

#                       t
#                       |
#                       |
#                       | 4
#                       @,
#                      /|`
#                      ||  `,
#                     / |    ',
#                     | |      \
#                    /  |       `.
#                    |  |         `,  10
#                   /   @ 8         `@
#                   |   |             \
#                  /    |              `.
#                  |    |                ',
#               9 @     |                  \
#                 |     @.,,_       7       `.
#                |     / 1   ``'-.,,@_        `.
#                |    /              ``''-.,,_  ', 3
#               |    /                        ``'@.,,,
#               |   '                       ,.-``     ``''- s
#              |  ,@ 5                 _,-'`
#              ' /                 ,.'`
#             | /             _.@``
#             '/          ,-'`   6
#            |/      ,.-``
#            /  _,-``
#          .@ '`
#         / 2
#        /
#       /
#      r
#

# natural coordinates
const coords_TET10 =
[  0.0   0.0   0.0
   1.0   0.0   0.0
   0.0   1.0   0.0
   0.0   0.0   1.0
   0.5   0.0   0.0
   0.5   0.5   0.0
   0.0   0.5   0.0
   0.0   0.0   0.5
   0.5   0.0   0.5
   0.0   0.5   0.5 ]


const facet_idxs_TET10 = [ [1, 4, 3, 8, 10, 7], [1, 2, 4, 5, 9, 8], [1, 3, 2, 7, 6, 5], [2, 3, 4, 6, 10, 9] ]
const edge_idxs_TET10 = [ [1, 2, 5], [2, 3, 6], [3, 1, 7], [1, 4, 8], [2, 4, 9], [3, 4, 10] ]

function shape_func_TET10(R::AbstractArray{<:Float64,1})
    r, s, t = R

    N = Array{Float64}(undef,10)

    u = 1.0 - r - s - t

    # corners
    N[1] = u*(2.0*u - 1.0)
    N[2] = r*(2.0*r - 1.0)
    N[3] = s*(2.0*s - 1.0)
    N[4] = t*(2.0*t - 1.0)

    # midedge
    N[5] = 4.0 * u * r
    N[6] = 4.0 * r * s
    N[7] = 4.0 * s * u
    N[8] = 4.0 * u * t
    N[9] = 4.0 * r * t
    N[10] = 4.0 * s * t

    return N
end

function shape_deriv_TET10(R::AbstractArray{<:Float64,1})
    r, s, t = R

    D = Array{Float64}(undef,10,3)

    # r-derivatives: dN0/dr to dN9/dr
    D[1,1]  =  4.0*(r + s + t) - 3.0
    D[2,1]  =  4.0*r - 1.0
    D[3,1]  =  0.0
    D[4,1]  =  0.0
    D[5,1]  =  4.0 - 8.0*r - 4.0*s - 4.0*t
    D[6,1]  =  4.0*s
    D[7,1]  = -4.0*s
    D[8,1]  = -4.0*t
    D[9,1]  =  4.0*t
    D[10,1] =  0.0

    # s-derivatives: dN0/ds to dN9/ds
    D[1,2]  =  4.0*(r + s + t) - 3.0
    D[2,2]  =  0.0
    D[3,2]  =  4.0*s - 1.0
    D[4,2]  =  0.0
    D[5,2]  = -4.0*r
    D[6,2]  =  4.0*r
    D[7,2]  =  4.0 - 4.0*r - 8.0*s - 4.0*t
    D[8,2]  = -4.0*t
    D[9,2]  =  0.0
    D[10,2] =  4.0*t

    # t-derivatives: dN0/dt to dN9/dt
    D[1,3]  =  4.0*(r + s + t) - 3.0
    D[2,3]  =  0.0
    D[3,3]  =  0.0
    D[4,3]  =  4.0*t - 1.0
    D[5,3]  = -4.0*r
    D[6,3]  =  0.0
    D[7,3]  = -4.0*s
    D[8,3]  =  4.0 - 4.0*r - 4.0*s - 8.0*t
    D[9,3]  =  4.0*r
    D[10,3] =  4.0*s

    return D
end

# constructor
function MakeTET10()
    shape             = CellShape()
    shape.name        = "TET10"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 10
    shape.basic_shape = TET4
    shape.vtk_type    = VTK_QUADRATIC_TETRA
    shape.facet_idxs  = facet_idxs_TET10
    shape.edge_idxs   = edge_idxs_TET10
    shape.facet_shape = TRI6
    shape.nat_coords  = coords_TET10
    shape.quadrature  = Dict( 0 => TET_IP4, 1 => TET_IP1, 4 => TET_IP4, 5 => TET_IP5, 11 => TET_IP11 )
    shape.func        = shape_func_TET10
    shape.deriv       = shape_deriv_TET10
    return shape
end


# Registration
const  TET10 = MakeTET10()
export TET10

# PYR5 shape
# ==========


# natural coordinates
const coords_PYR5 =
[ -1.0 -1.0  0.0
   1.0 -1.0  0.0
   1.0  1.0  0.0
  -1.0  1.0  0.0
   0.0  0.0  1.0  ]

const facet_idxs_PYR5 =
    [ [1, 4, 3, 2],     [1, 2, 5],     [1, 5, 4],      [3, 4, 5],     [2, 3, 5] ]

const edge_idxs_PYR5 =
    [ [1, 2],    [2, 3],    [3, 4],     [4, 1],      [1, 5],      [2, 5],      [3, 5],      [4, 5] ]

function shape_func_PYR5(R::AbstractArray{<:Float64,1})
    r, s, t = R
    w = t==1.0 ? 0.0 : 1/(1-t)

    N = Array{Float64}(undef,5)
    N[1] = 0.25*(1-r-s-t+r*s*w)
    N[2] = 0.25*(1+r-s-t-r*s*w)
    N[3] = 0.25*(1+r+s-t+r*s*w)
    N[4] = 0.25*(1-r+s-t-r*s*w)
    N[5] = t
    return N
end

function shape_deriv_PYR5(R::AbstractArray{<:Float64,1})
    r, s, t = R
    w = t==1.0 ? 0.0 : 1/(1-t)


    D = Array{Float64}(undef,5,3)
    D[1,1] = 0.25*(-1+s*w)  ;   D[1,2]= 0.25*(-1+r*w)  ;   D[1,3]= 0.25*(-1+r*s*w^2)
    D[2,1] = 0.25*(1-s*w)   ;   D[2,2]= 0.25*(-1-r*w)  ;   D[2,3]= 0.25*(-1-r*s*w^2)
    D[3,1] = 0.25*(1+s*w)   ;   D[3,2]= 0.25*(1+r*w)   ;   D[3,3]= 0.25*(-1+r*s*w^2)
    D[4,1] = 0.25*(-1-s*w)  ;   D[4,2]= 0.25*(1-r*w)   ;   D[4,3]= 0.25*(-1-r*s*w^2)
    D[5,1] =  0.0           ;   D[5,2]=  0.0           ;   D[5,3]=  1.0
    return D
end

# constructor
function MakePYR5()
    shape             = CellShape()
    shape.name        = "PYR5"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 5
    shape.basic_shape = shape
    shape.vtk_type    = VTK_PYRAMID
    shape.facet_idxs  = facet_idxs_PYR5
    shape.edge_idxs   = edge_idxs_PYR5
    shape.facet_shape = (QUAD4, TRI3, TRI3, TRI3, TRI3)
    shape.nat_coords  = coords_PYR5
    shape.quadrature  = Dict( 0 => PYR_IP5,  5 => PYR_IP5,  8 => PYR_IP8 )
    shape.func        = shape_func_PYR5
    shape.deriv       = shape_deriv_PYR5
    return shape
end

# Registration
const  PYR5 = MakePYR5()
export PYR5


# PYR13 shape
# ===========

# G. Bedrosian. Shape functions and integration formulas for
# three-dimensional finite element analysis.
# Int. J. Numerical Methods Engineering, vol 35, p. 95-108, 1992.

# natural coordinates
const coords_PYR13 =
[ -1.0 -1.0  0.0
   1.0 -1.0  0.0
   1.0  1.0  0.0
  -1.0  1.0  0.0
   0.0  0.0  1.0
   0.0 -1.0  0.0
   1.0  0.0  0.0
   0.0  1.0  0.0
  -1.0  0.0  0.0
  -0.5 -0.5  0.5
   0.5 -0.5  0.5
   0.5  0.5  0.5
  -0.5  0.5  0.5 ]

const facet_idxs_PYR13 =
    [ [1, 4, 3, 2, 9, 8, 7, 6],     [1, 2, 5, 6, 11, 10],     [1, 5, 4, 10, 13, 9],      [3, 4, 5, 8, 13, 12],     [2, 3, 5, 7, 12, 11] ]

const edge_idxs_PYR13 =
    [ [1, 2, 6],    [2, 3, 7],    [3, 4, 8],     [4, 1, 9],      [1, 5, 10],      [2, 5, 11],      [3, 5, 12],      [4, 5, 13] ]

function shape_func_PYR13(R::AbstractArray{<:Float64,1})
    r, s, t = R
    if r==s==0.0 && t==1.0
        return [ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    end
    t==1 && (t=0.9999999999999999) # undetermined at t==1; tested for Float64
    w = 1/(1-t)

    N = Array{Float64}(undef,13)
    N[1]  = 1/4*(-r-s-1)*((1-r)*(1-s)-t+r*s*t*w)
    N[2]  = 1/4*(-s+r-1)*((1+r)*(1-s)-t-r*s*t*w)
    N[3]  = 1/4*(r+s-1)*((1+r)*(1+s)-t+r*s*t*w)
    N[4]  = 1/4*(s-r-1)*((1-r)*(1+s)-t-r*s*t*w)
    N[5]  = t*(2*t-1)
    N[6]  = 1/2*(1+r-t)*(1-r-t)*(1-s-t)*w
    N[7]  = 1/2*(1+s-t)*(1-s-t)*(1+r-t)*w
    N[8]  = 1/2*(1+r-t)*(1-r-t)*(1+s-t)*w
    N[9]  = 1/2*(1+s-t)*(1-s-t)*(1-r-t)*w
    N[10] = t*(1-r-t)*(1-s-t)*w
    N[11] = t*(1+r-t)*(1-s-t)*w
    N[12] = (1+s-t)*(1+r-t)*t*w
    N[13] = t*(1-r-t)*(1+s-t)*w
    return N
end

function shape_deriv_PYR13(R::AbstractArray{<:Float64,1})
    r, s, t = R
    if r==s==0.0 && t==1.0
        return [
             0.25   0.25   0.25
            -0.25   0.25   0.25
            -0.25  -0.25   0.25
             0.25  -0.25   0.25
             0.0    0.0    3.0
             0.0    0.0    0.0
             0.0    0.0    0.0
             0.0    0.0    0.0
             0.0    0.0    0.0
            -1.0   -1.0   -1.0
             1.0   -1.0   -1.0
             1.0    1.0   -1.0
            -1.0    1.0   -1.0
        ]       
    end

    t==1 && (t=0.99999) # undetermined at t==1; tested for Float64
    w = 1/(1-t)

    D = Array{Float64}(undef,13,3)
    D[1,1]  = -1/4*(-t-s+2*s*t-2*r+2*t*r+2*s*r+t^2+s^2)*w
    D[2,1]  =  1/4*(-t-s+2*s*t+2*r-2*t*r-2*s*r+t^2+s^2)*w
    D[3,1]  =  1/4*(-t+s-2*s*t+2*r-2*t*r+2*s*r+t^2+s^2)*w
    D[4,1]  = -1/4*(-t+s-2*s*t-2*r+2*t*r-2*s*r+t^2+s^2)*w
    D[5,1]  =  0.0
    D[6,1]  =  (-1+s+t)*r*w
    D[7,1]  = -1/2*(-1+s+t)*(1+s-t)*w
    D[8,1]  = -(1+s-t)*r*w
    D[9,1]  =  1/2*(-1+s+t)*(1+s-t)*w
    D[10,1] =  (-1+s+t)*t*w
    D[11,1] = -(-1+s+t)*t*w
    D[12,1] =  (1+s-t)*t*w
    D[13,1] = -(1+s-t)*t*w

    D[1,2]  = -1/4*(-t-2*s+2*s*t-r+2*t*r+2*s*r+t^2+r^2)*w
    D[2,2]  =  1/4*(t+2*s-2*s*t-r+2*t*r+2*s*r-t^2-r^2)*w
    D[3,2]  =  1/4*(-t+2*s-2*s*t+r-2*t*r+2*s*r+t^2+r^2)*w
    D[4,2]  = -1/4*(t-2*s+2*s*t+r-2*t*r+2*s*r-t^2-r^2)*w
    D[5,2]  =  0.0
    D[6,2]  =  1/2*(-1+r+t)*(1+r-t)*w
    D[7,2]  = -(1+r-t)*s*w
    D[8,2]  = -1/2*(-1+r+t)*(1+r-t)*w
    D[9,2]  =  (-1+r+t)*s*w
    D[10,2] =  (-1+r+t)*t*w
    D[11,2] = -(1+r-t)*t*w
    D[12,2] =  (1+r-t)*t*w
    D[13,2] = -(-1+r+t)*t*w

    D[1,3]  = -1/4*(r+s+1)*(-1+2*t-t^2+s*r)*w^2
    D[2,3]  =  1/4*(s-r+1)*(1-2*t+t^2+s*r)*w^2
    D[3,3]  =  1/4*(r+s-1)*(-1+2*t-t^2+s*r)*w^2
    D[4,3]  = -1/4*(s-r-1)*(1-2*t+t^2+s*r)*w^2
    D[5,3]  =  4*t-1
    D[6,3]  =  1/2*(-2+s+6*t+s*r^2+s*t^2-6*t^2+2*t^3-2*s*t)*w^2
    D[7,3]  = -1/2*(2-6*t+r+r*t^2+s^2*r+6*t^2-2*t^3-2*t*r)*w^2
    D[8,3]  = -1/2*(2+s-6*t+s*r^2+s*t^2+6*t^2-2*t^3-2*s*t)*w^2
    D[9,3]  =  1/2*(-2+6*t+r+r*t^2+s^2*r-6*t^2+2*t^3-2*t*r)*w^2
    D[10,3] =  (1-s-4*t-r-r*t^2-s*t^2+s*r+5*t^2-2*t^3+2*s*t+2*t*r)*w^2
    D[11,3] = -(-1+s+4*t-r-r*t^2+s*t^2+s*r-5*t^2+2*t^3-2*s*t+2*t*r)*w^2
    D[12,3] =  (1+s-4*t+r+r*t^2+s*t^2+s*r+5*t^2-2*t^3-2*s*t-2*t*r)*w^2
    D[13,3] = -(-1-s+4*t+r+r*t^2-s*t^2+s*r-5*t^2+2*t^3+2*s*t-2*t*r)*w^2

    return D
end

# constructor
function MakePYR13()
    shape             = CellShape()
    shape.name        = "PYR13"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 13
    shape.basic_shape = PYR5
    shape.vtk_type    = VTK_QUADRATIC_PYRAMID
    shape.facet_idxs  = facet_idxs_PYR13
    shape.edge_idxs   = edge_idxs_PYR13
    shape.facet_shape = (QUAD8, TRI6, TRI6, TRI6, TRI6)
    shape.nat_coords  = coords_PYR13
    shape.quadrature  = Dict( 0 => PYR_IP5,  5 => PYR_IP5,  8 => PYR_IP8 )
    shape.func        = shape_func_PYR13
    shape.deriv       = shape_deriv_PYR13
    return shape
end


# Registration
const  PYR13 = MakePYR13()
export PYR13

# HEX8 shape
# ==========

# Local IDs
#                  Nodes                                   Faces
#     z
#     |           5                  8
#    ,+--y         @________________@                    +________________+
#  x'            ,'|              ,'|                  ,'|              ,'|
#              ,'  |            ,'  |                ,'  |  ___       ,'  |
#            ,'    |          ,'    |              ,'    |,'5,'  [0],'    |
#      6   ,'      |      7 ,'      |            ,'      |~~~     ,'      |
#        @'===============@'        |          +'===============+'  ,'|   |
#        |         |      |         |          |   ,'|   |      |   |3|   |
#        |         |      |         |          |   |2|   |      |   |,'   |
#        |       1 @______|_________@          |   |,'   +______|_________+
#        |       ,'       |       ,' 4         |       ,'       |       ,'
#        |     ,'         |     ,'             |     ,' [1]  ___|     ,'
#        |   ,'           |   ,'               |   ,'      ,'4,'|   ,'
#        | ,'             | ,'                 | ,'        ~~~  | ,'
#        @________________@'                   +________________+'
#      2                   3

# natural coordinates
const coords_HEX8 =
[ -1.0 -1.0 -1.0
   1.0 -1.0 -1.0
   1.0  1.0 -1.0
  -1.0  1.0 -1.0
  -1.0 -1.0  1.0
   1.0 -1.0  1.0
   1.0  1.0  1.0
  -1.0  1.0  1.0 ]

const facet_idxs_HEX8 = [ [1, 5, 8, 4], [2, 3, 7, 6], [1, 2, 6, 5], [3, 4, 8, 7], [1, 4, 3, 2], [5, 6, 7, 8] ]
const edge_idxs_HEX8 = [ [1, 2], [2, 3], [3, 4], [4, 1], [5, 6], [6, 7], [7, 8], [8, 5], [1, 5], [2, 6], [3, 7], [4, 8] ]

function shape_func_HEX8(R::AbstractArray{<:Float64,1})
    r, s, t = R[1:3]
    N = Array{Float64}(undef,8)
    N[1] = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t)
    N[2] = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t)
    N[3] = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t)
    N[4] = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t)
    N[5] = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t)
    N[6] = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t)
    N[7] = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t)
    N[8] = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t)
    return N
end

function shape_deriv_HEX8(R::AbstractArray{<:Float64,1})
    r, s, t = R
    st = s*t
    rt = r*t
    rs = r*s
    D = Array{Float64}(undef,8,3)
    D[1,1] = -1.0+s+t-st;   D[1,2]=-1.0+r+t-rt;   D[1,3]=-1.0+r+s-rs
    D[2,1] = +1.0-s-t+st;   D[2,2]=-1.0-r+t+rt;   D[2,3]=-1.0-r+s+rs
    D[3,1] = +1.0+s-t-st;   D[3,2]=+1.0+r-t-rt;   D[3,3]=-1.0-r-s-rs
    D[4,1] = -1.0-s+t+st;   D[4,2]=+1.0-r-t+rt;   D[4,3]=-1.0+r-s+rs
    D[5,1] = -1.0+s-t+st;   D[5,2]=-1.0+r-t+rt;   D[5,3]=+1.0-r-s+rs
    D[6,1] = +1.0-s+t-st;   D[6,2]=-1.0-r-t-rt;   D[6,3]=+1.0+r-s-rs
    D[7,1] = +1.0+s+t+st;   D[7,2]=+1.0+r+t+rt;   D[7,3]=+1.0+r+s+rs
    D[8,1] = -1.0-s-t-st;   D[8,2]=+1.0-r+t-rt;   D[8,3]=+1.0-r+s-rs
    D = 0.125*D

    return D
end

# constructor
function MakeHEX8()
    shape             = CellShape()
    shape.name        = "HEX8"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 8
    shape.basic_shape = shape
    shape.vtk_type    = VTK_HEXAHEDRON
    shape.facet_idxs  = facet_idxs_HEX8
    shape.edge_idxs   = edge_idxs_HEX8
    shape.facet_shape = QUAD4
    shape.nat_coords  = coords_HEX8
    shape.quadrature  = Dict( 0 => HEX_IP2, 8 => HEX_IP2, 27 => HEX_IP3 )
    shape.func        = shape_func_HEX8
    shape.deriv       = shape_deriv_HEX8
    return shape
end


# Registration
const  HEX8 = MakeHEX8()
export HEX8



# HEX20 shape
# ===========

# Local IDs
#                   Vertices                               Faces
#     t
#     |           5        16        8
#    ,+--s         @-------@--------@                   +----------------+
#  r'            ,'|              ,'|                 ,'|              ,'|
#           13 @'  |         15 ,'  |               ,'  |  ___       ,'  |
#            ,'    |17        ,@    |20           ,'    |,'6,'  [1],'    |
#      6   ,'      @      7 ,'      @           ,'      |~~~     ,'      |
#        @'=======@=======@'        |         +'===============+'  ,'|   |
#        |      14 |      |         |         |   ,'|   |      |   |4|   |
#        |         |      |  12     |         |   |3|   |      |   |,'   |
#     18 |       1 @- - - | @- - - -@         |   |,'   +- - - | +- - - -+
#        @       ,'       @       ,' 4        |       ,'       |       ,'
#        |   9 @'      19 |     ,'            |     ,' [2]  ___|     ,'
#        |   ,'           |   ,@ 11           |   ,'      ,'5,'|   ,'
#        | ,'             | ,'                | ,'        ~~~  | ,'
#        @-------@--------@'                  +----------------+'
#      2        10         3


# natural coordinates
const coords_HEX20 =
[ -1.0 -1.0 -1.0
   1.0 -1.0 -1.0
   1.0  1.0 -1.0
  -1.0  1.0 -1.0
  -1.0 -1.0  1.0
   1.0 -1.0  1.0
   1.0  1.0  1.0
  -1.0  1.0  1.0

   0.0 -1.0 -1.0
   1.0  0.0 -1.0
   0.0  1.0 -1.0
  -1.0  0.0 -1.0
   0.0 -1.0  1.0
   1.0  0.0  1.0
   0.0  1.0  1.0
  -1.0  0.0  1.0

  -1.0 -1.0  0.0
   1.0 -1.0  0.0
   1.0  1.0  0.0
  -1.0  1.0  0.0 ]

const facet_idxs_HEX20 = [ [1, 5, 8, 4,17,16,20,12], [2, 3, 7, 6, 10,19,14,18], [1, 2, 6, 5, 9,18,13,17], [3, 4, 8, 7,11,20,15,19], [1, 4, 3, 2,12,11, 10, 9], [5, 6, 7, 8,13,14,15,16] ]
const edge_idxs_HEX20 = [ [1, 2, 9], [2, 3, 10], [3, 4, 11], [4, 1, 12], [5, 6, 13], [6, 7, 14], [7, 8, 15], [8, 5, 16], [1, 5, 17], [2, 6, 18], [3, 7, 19], [4, 8, 20] ]

function shape_func_HEX20(R::AbstractArray{<:Float64,1})
    r, s, t = R

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    N = Array{Float64}(undef,20)
    N[ 1] = 0.125*rm1*sm1*tm1*(-r-s-t-2.0)
    N[ 2] = 0.125*rp1*sm1*tm1*( r-s-t-2.0)
    N[ 3] = 0.125*rp1*sp1*tm1*( r+s-t-2.0)
    N[ 4] = 0.125*rm1*sp1*tm1*(-r+s-t-2.0)
    N[ 5] = 0.125*rm1*sm1*tp1*(-r-s+t-2.0)
    N[ 6] = 0.125*rp1*sm1*tp1*( r-s+t-2.0)
    N[ 7] = 0.125*rp1*sp1*tp1*( r+s+t-2.0)
    N[ 8] = 0.125*rm1*sp1*tp1*(-r+s+t-2.0)
    N[ 9] = 0.25*(1.0-r*r)*sm1*tm1
    N[10] = 0.25*rp1*(1.0-s*s)*tm1
    N[11] = 0.25*(1.0-r*r)*sp1*tm1
    N[12] = 0.25*rm1*(1.0-s*s)*tm1
    N[13] = 0.25*(1.0-r*r)*sm1*tp1
    N[14] = 0.25*rp1*(1.0-s*s)*tp1
    N[15] = 0.25*(1.0-r*r)*sp1*tp1
    N[16] = 0.25*rm1*(1.0-s*s)*tp1
    N[17] = 0.25*rm1*sm1*(1.0-t*t)
    N[18] = 0.25*rp1*sm1*(1.0-t*t)
    N[19] = 0.25*rp1*sp1*(1.0-t*t)
    N[20] = 0.25*rm1*sp1*(1.0-t*t)
    return N
end

function shape_deriv_HEX20(R::AbstractArray{<:Float64,1})
    r, s, t = R

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    D = Array{Float64}(undef,20,3)
    # Derivatives with respect to r
    D[ 1,1] = -0.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[ 2,1] =  0.125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1
    D[ 3,1] =  0.125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[ 4,1] = -0.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[ 5,1] = -0.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[ 6,1] =  0.125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[ 7,1] =  0.125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[ 8,1] = -0.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1
    D[ 9,1] = -0.5*r*sm1*tm1
    D[10,1] =  0.25*(1-s*s)*tm1
    D[11,1] = -0.5*r*sp1*tm1
    D[12,1] = -0.25*(1-s*s)*tm1
    D[13,1] = -0.5*r*sm1*tp1
    D[14,1] =  0.25*(1-s*s)*tp1
    D[15,1] = -0.5*r*sp1  *tp1
    D[16,1] = -0.25*(1-s*s)*tp1
    D[17,1] = -0.25*sm1*(1-t*t)
    D[18,1] =  0.25*sm1*(1-t*t)
    D[19,1] =  0.25*sp1*(1-t*t)
    D[20,1] = -0.25*sp1*(1-t*t)

    # Derivatives with respect to s
    D[ 1,2] = -0.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[ 2,2] = -0.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[ 3,2] =  0.125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[ 4,2] =  0.125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1
    D[ 5,2] = -0.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[ 6,2] = -0.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1
    D[ 7,2] =  0.125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[ 8,2] =  0.125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[ 9,2] = -0.25*(1-r*r)*tm1
    D[10,2] = -0.5*s*rp1*tm1
    D[11,2] =  0.25*(1-r*r)*tm1
    D[12,2] = -0.5*s*rm1*tm1
    D[13,2] = -0.25*(1-r*r)*tp1
    D[14,2] = -0.5*s*rp1*tp1
    D[15,2] =  0.25*(1-r*r)*tp1
    D[16,2] = -0.5*s*rm1*tp1
    D[17,2] = -0.25*rm1*(1-t*t)
    D[18,2] = -0.25*rp1*(1-t*t)
    D[19,2] =  0.25*rp1*(1-t*t)
    D[20,2] =  0.25*rm1*(1-t*t)

    # Derivatives with respect to t
    D[ 1,3] = -0.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[ 2,3] = -0.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[ 3,3] = -0.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1
    D[ 4,3] = -0.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[ 5,3] =  0.125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1
    D[ 6,3] =  0.125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[ 7,3] =  0.125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[ 8,3] =  0.125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[ 9,3] = -0.25*(1-r*r)*sm1
    D[10,3] = -0.25*rp1*(1-s*s)
    D[11,3] = -0.25*(1-r*r)*sp1
    D[12,3] = -0.25*rm1*(1-s*s)
    D[13,3] =  0.25*(1-r*r)*sm1
    D[14,3] =  0.25*rp1*(1-s*s)
    D[15,3] =  0.25*(1-r*r)*sp1
    D[16,3] =  0.25*rm1*(1-s*s)
    D[17,3] = -0.5*t*rm1*sm1
    D[18,3] = -0.5*t*rp1*sm1
    D[19,3] = -0.5*t*rp1*sp1
    D[20,3] = -0.5*t*rm1*sp1

    return D
end

# constructor
function MakeHEX20()
    shape             = CellShape()
    shape.name        = "HEX20"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 20
    shape.basic_shape = HEX8
    shape.vtk_type    = VTK_QUADRATIC_HEXAHEDRON
    shape.facet_idxs  = facet_idxs_HEX20
    shape.edge_idxs   = edge_idxs_HEX20
    shape.facet_shape = QUAD8
    shape.nat_coords  = coords_HEX20
    shape.quadrature  = Dict( 0 => HEX_IP2, 8 => HEX_IP2, 27 => HEX_IP3 ) # Note: Reduced integration may produce spurius modes
    shape.func        = shape_func_HEX20
    shape.deriv       = shape_deriv_HEX20
    return shape
end


# Registration
const  HEX20 = MakeHEX20()
export HEX20




# HEX27 shape
# ===========

# Local IDs
#                   Vertices                               Faces
#     t
#     |           5        16        8
#    ,+--s         @-------@--------@                   +----------------+
#  r'            ,'|    26        ,'|                 ,'|              ,'|
#           13 @'  |     #   15 ,'  |               ,'  |  ___       ,'  |
#            ,'    |17   21   ,@    |             ,'    |,'6,'  [1],'    |
#      6   ,'      @      # ,'      @ 20        ,'      |~~~     ,'      |
#        @'=======@=======@'   24   |         +'===============+'  ,'|   |
#        |      14 |  27  |7   #    |         |   ,'|   |      |   |4|   |
#        |    #    |  @   |  12     |         |   |3|   |      |   |,'   |
#        |   23  1 @- - - | @- - - -@         |   |,'   +- - - | +- - - -+
#     18 @       ,'#      @       ,' 4        |       ,'       |       ,'
#        |   9 @'  22  #  |19   ,'            |     ,' [2]  ___|     ,'
#        |   ,'       25  |   ,@ 11           |   ,'      ,'5,'|   ,'
#        | ,'             | ,'                | ,'        ~~~  | ,'
#        @-------@--------@'                  +----------------+'
#      2        10         3


# natural coordinates
const coords_HEX27 =
[
  # corner points
  -1.0  -1.0  -1.0
   1.0  -1.0  -1.0
   1.0   1.0  -1.0
  -1.0   1.0  -1.0
  -1.0  -1.0   1.0
   1.0  -1.0   1.0
   1.0   1.0   1.0
  -1.0   1.0   1.0

  # mid-edge points
   0.0  -1.0  -1.0
   1.0   0.0  -1.0
   0.0   1.0  -1.0
  -1.0   0.0  -1.0
   0.0  -1.0   1.0
   1.0   0.0   1.0
   0.0   1.0   1.0
  -1.0   0.0   1.0

  -1.0  -1.0   0.0
   1.0  -1.0   0.0
   1.0   1.0   0.0
  -1.0   1.0   0.0

  # face center points
  -1.0   0.0   0.0
   1.0   0.0   0.0
   0.0  -1.0   0.0
   0.0   1.0   0.0
   0.0   0.0  -1.0
   0.0   0.0   1.0

  # element center point
   0.0   0.0   0.0 ]


const facet_idxs_HEX27 = [[1, 5, 8, 4, 17, 16, 20, 12, 21], [2, 3, 7, 6, 10, 19, 14, 18, 22], [1, 2, 6, 5, 9, 18, 13, 17, 23], [4, 8, 7, 3, 20, 15, 19, 11, 24], [1, 4, 3, 2, 12, 11, 10, 9, 25], [5, 6, 7, 8, 13, 14, 15, 16, 26]]
const edge_idxs_HEX27 = [[1, 2, 9], [2, 3, 10], [4, 3, 11], [1, 4, 12], [5, 6, 13], [6, 7, 14], [8, 7, 15], [5, 8, 16], [1, 5, 17], [2, 6, 18], [4, 8, 20], [3, 7, 19]]


function shape_func_HEX27(R::AbstractArray{<:Float64,1})
    r, s, t = R
  
    g1r = -0.5 * r * (1 - r)
    g1s = -0.5 * s * (1 - s)
    g1t = -0.5 * t * (1 - t)
  
    g2r = (1 + r)*(1 - r)
    g2s = (1 + s)*(1 - s)
    g2t = (1 + t)*(1 - t)
  
    g3r = 0.5 * r * (1 + r)
    g3s = 0.5 * s * (1 + s)
    g3t = 0.5 * t * (1 + t)
  
    N = Array{Float64}(undef,27)

    N[1] = g1r * g1s * g1t
    N[2] = g3r * g1s * g1t
    N[3] = g3r * g3s * g1t
    N[4] = g1r * g3s * g1t
    N[5] = g1r * g1s * g3t
    N[6] = g3r * g1s * g3t
    N[7] = g3r * g3s * g3t
    N[8] = g1r * g3s * g3t
   
    N[9] =  g2r * g1s * g1t
    N[10] =  g3r * g2s * g1t
    N[11] = g2r * g3s * g1t
    N[12] = g1r * g2s * g1t
    N[13] = g2r * g1s * g3t
    N[14] = g3r * g2s * g3t
    N[15] = g2r * g3s * g3t
    N[16] = g1r * g2s * g3t
    N[17] = g1r * g1s * g2t
    N[18] = g3r * g1s * g2t
    N[19] = g3r * g3s * g2t
    N[20] = g1r * g3s * g2t
   
    N[21] = g1r * g2s * g2t
    N[22] = g3r * g2s * g2t
    N[23] = g2r * g1s * g2t
    N[24] = g2r * g3s * g2t
    N[25] = g2r * g2s * g1t
    N[26] = g2r * g2s * g3t
   
    N[27] = g2r * g2s * g2t

    return N
end

function shape_deriv_HEX27(R::AbstractArray{<:Float64,1})
    r, s, t = R

    g1r = -0.5 * r * (1 - r);
    g1s = -0.5 * s * (1 - s);
    g1t = -0.5 * t * (1 - t);
   
    g2r = (1 + r)*(1 - r);
    g2s = (1 + s)*(1 - s);
    g2t = (1 + t)*(1 - t);
   
    g3r = 0.5 * r * (1 + r);
    g3s = 0.5 * s * (1 + s);
    g3t = 0.5 * t * (1 + t);
   
    g1r_r = r - 0.5;
    g1s_s = s - 0.5;
    g1t_t = t - 0.5;
   
    g2r_r = -2*r;
    g2s_s = -2*s;
    g2t_t = -2*t;
   
    g3r_r = r + 0.5;
    g3s_s = s + 0.5;
    g3t_t = t + 0.5;
   
    D = Array{Float64}(undef,27,3)
   
    # r-derivatives
    D[ 1,1] = g1r_r * g1s * g1t
    D[ 2,1] = g3r_r * g1s * g1t
    D[ 3,1] = g3r_r * g3s * g1t
    D[ 4,1] = g1r_r * g3s * g1t
    D[ 5,1] = g1r_r * g1s * g3t
    D[ 6,1] = g3r_r * g1s * g3t
    D[ 7,1] = g3r_r * g3s * g3t
    D[ 8,1] = g1r_r * g3s * g3t
    D[ 9,1] = g2r_r * g1s * g1t
    D[10,1] = g3r_r * g2s * g1t
    D[11,1] = g2r_r * g3s * g1t
    D[12,1] = g1r_r * g2s * g1t
    D[13,1] = g2r_r * g1s * g3t
    D[14,1] = g3r_r * g2s * g3t
    D[15,1] = g2r_r * g3s * g3t
    D[16,1] = g1r_r * g2s * g3t
    D[17,1] = g1r_r * g1s * g2t
    D[18,1] = g3r_r * g1s * g2t
    D[19,1] = g3r_r * g3s * g2t
    D[20,1] = g1r_r * g3s * g2t
    D[21,1] = g1r_r * g2s * g2t
    D[22,1] = g3r_r * g2s * g2t
    D[23,1] = g2r_r * g1s * g2t
    D[24,1] = g2r_r * g3s * g2t
    D[25,1] = g2r_r * g2s * g1t
    D[26,1] = g2r_r * g2s * g3t
    D[27,1] = g2r_r * g2s * g2t
   
    # s-derivatives
    D[ 1,2] =  g1r * g1s_s * g1t
    D[ 2,2] =  g3r * g1s_s * g1t
    D[ 3,2] =  g3r * g3s_s * g1t
    D[ 4,2] =  g1r * g3s_s * g1t
    D[ 5,2] =  g1r * g1s_s * g3t
    D[ 6,2] =  g3r * g1s_s * g3t
    D[ 7,2] =  g3r * g3s_s * g3t
    D[ 8,2] =  g1r * g3s_s * g3t
    D[ 9,2] =  g2r * g1s_s * g1t
    D[10,2] =  g3r * g2s_s * g1t
    D[11,2] =  g2r * g3s_s * g1t
    D[12,2] =  g1r * g2s_s * g1t
    D[13,2] =  g2r * g1s_s * g3t
    D[14,2] =  g3r * g2s_s * g3t
    D[15,2] =  g2r * g3s_s * g3t
    D[16,2] =  g1r * g2s_s * g3t
    D[17,2] =  g1r * g1s_s * g2t
    D[18,2] =  g3r * g1s_s * g2t
    D[19,2] =  g3r * g3s_s * g2t
    D[20,2] =  g1r * g3s_s * g2t
    D[21,2] =  g1r * g2s_s * g2t
    D[22,2] =  g3r * g2s_s * g2t
    D[23,2] =  g2r * g1s_s * g2t
    D[24,2] =  g2r * g3s_s * g2t
    D[25,2] =  g2r * g2s_s * g1t
    D[26,2] =  g2r * g2s_s * g3t
    D[27,2] =  g2r * g2s_s * g2t
   
    # t-derivatives
    D[ 1,3] =  g1r * g1s * g1t_t
    D[ 2,3] =  g3r * g1s * g1t_t
    D[ 3,3] =  g3r * g3s * g1t_t
    D[ 4,3] =  g1r * g3s * g1t_t
    D[ 5,3] =  g1r * g1s * g3t_t
    D[ 6,3] =  g3r * g1s * g3t_t
    D[ 7,3] =  g3r * g3s * g3t_t
    D[ 8,3] =  g1r * g3s * g3t_t
    D[ 9,3] =  g2r * g1s * g1t_t
    D[10,3] =  g3r * g2s * g1t_t
    D[11,3] =  g2r * g3s * g1t_t
    D[12,3] =  g1r * g2s * g1t_t
    D[13,3] =  g2r * g1s * g3t_t
    D[14,3] =  g3r * g2s * g3t_t
    D[15,3] =  g2r * g3s * g3t_t
    D[16,3] =  g1r * g2s * g3t_t
    D[17,3] =  g1r * g1s * g2t_t
    D[18,3] =  g3r * g1s * g2t_t
    D[19,3] =  g3r * g3s * g2t_t
    D[20,3] =  g1r * g3s * g2t_t
    D[21,3] =  g1r * g2s * g2t_t
    D[22,3] =  g3r * g2s * g2t_t
    D[23,3] =  g2r * g1s * g2t_t
    D[24,3] =  g2r * g3s * g2t_t
    D[25,3] =  g2r * g2s * g1t_t
    D[26,3] =  g2r * g2s * g3t_t
    D[27,3] =  g2r * g2s * g2t_t

    return D
end

# constructor
function MakeHEX27()
    shape             = CellShape()
    shape.name        = "HEX27"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 27
    shape.basic_shape = HEX8
    shape.vtk_type    = VTK_TRIQUADRATIC_HEXAHEDRON
    shape.facet_idxs  = facet_idxs_HEX27
    shape.edge_idxs   = edge_idxs_HEX27
    shape.facet_shape = QUAD9
    shape.nat_coords  = coords_HEX27
    shape.quadrature  = Dict( 0 => HEX_IP3, 8 => HEX_IP2, 27 => HEX_IP3 ) # Note: Reduced integration may produce spurius modes
    shape.func        = shape_func_HEX27
    shape.deriv       = shape_deriv_HEX27
    return shape
end


# Registration
const  HEX27 = MakeHEX27()
export HEX27



# WED6 shape
# ==========

# natural coordinates
const coords_WED6 =
[ 0.0  0.0 -1.0
  1.0  0.0 -1.0
  0.0  1.0 -1.0
  0.0  0.0  1.0
  1.0  0.0  1.0
  0.0  1.0  1.0 ]

const facet_idxs_WED6 = [ [1, 4, 6, 3], [1, 2, 5, 4], [2, 3, 6, 5], [1, 3, 2], [4, 5, 6]]
const edge_idxs_WED6 = [ [1, 2], [2, 3], [3, 1], [4, 5], [5, 6], [6, 4], [1, 4], [2, 5], [3, 6] ]

function shape_func_WED6(R::AbstractArray{<:Float64,1})
    r, s, t = R[1:3]
    N = Array{Float64}(undef,6)
    N[1] = 0.5*(1.0-r-s-t+r*t+s*t)
    N[2] = 0.5*(r-r*t)
    N[3] = 0.5*(s-s*t)
    N[4] = 0.5*(1.0-r-s+t-r*t-s*t)
    N[5] = 0.5*(r+r*t)
    N[6] = 0.5*(s+s*t)
    return N
end

function shape_deriv_WED6(R::AbstractArray{<:Float64,1})
    r, s, t = R
    D = Array{Float64}(undef,6,3)
    D[1,1] = 0.5*(-1.0+t);  D[1,2] = 0.5*(-1.0+t);  D[1,3] = 0.5*(-1.0+r+s)
    D[2,1] = 0.5*(1.0-t) ;  D[2,2] = 0.0         ;  D[2,3] = 0.5*(-r)
    D[3,1] = 0.0         ;  D[3,2] = 0.5*(1.0-t) ;  D[3,3] = 0.5*(-s)
    D[4,1] = 0.5*(-1.0-t);  D[4,2] = 0.5*(-1.0-t);  D[4,3] = 0.5*(1.0-r-s)
    D[5,1] = 0.5*(1.0+t) ;  D[5,2] = 0.0         ;  D[5,3] = 0.5*(r)
    D[6,1] = 0.0         ;  D[6,2] = 0.5*(1.0+t) ;  D[6,3] = 0.5*(s)

    return D
end

# constructor
function MakeWED6()
    shape             = CellShape()
    shape.name        = "WED6"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 6
    shape.basic_shape = shape
    shape.vtk_type    = VTK_WEDGE
    shape.facet_idxs  = facet_idxs_WED6
    shape.edge_idxs   = edge_idxs_WED6
    shape.facet_shape = (QUAD4, QUAD4, QUAD4, TRI3, TRI3)
    shape.nat_coords  = coords_WED6
    shape.quadrature  = Dict( 0 => WED_IP9, 2 => WED_IP2, 9 => WED_IP9, 18 => WED_IP18 )
    shape.func        = shape_func_WED6
    shape.deriv       = shape_deriv_WED6
    return shape
end


# Registration
const  WED6 = MakeWED6()
export WED6



# WED15 shape
# ==========

# natural coordinates
const coords_WED15 =
[ 0.0  0.0 -1.0
  1.0  0.0 -1.0
  0.0  1.0 -1.0
  0.0  0.0  1.0
  1.0  0.0  1.0
  0.0  1.0  1.0

  0.5  0.0 -1.0
  0.5  0.5 -1.0
  0.0  0.5 -1.0
  0.5  0.0  1.0
  0.5  0.5  1.0
  0.0  0.5  1.0
  0.0  0.0  0.0
  1.0  0.0  0.0
  0.0  1.0  0.0 ]


const facet_idxs_WED15 = [ [1, 4, 6, 3, 13, 12, 15, 9], [1, 2, 5, 4, 7, 14, 10, 13], [2, 3, 6, 5, 8, 15, 11, 14], [1, 3, 2, 9, 8, 7], [4, 5, 6, 10, 11, 12]]
const edge_idxs_WED15 = [ [1, 2, 7], [2, 3, 8], [3, 1, 9], [4, 5, 10], [5, 6, 11], [6, 4, 12], [1, 4, 13], [2, 5, 14], [3, 6, 15] ]

function shape_func_WED15(R::AbstractArray{<:Float64,1})
    r, s, t = R[1:3]
    N = Array{Float64}(undef,15)
    N[1]  = 0.5*(-2.0*(r^2.0)*t-4.0*r*s*t-r*(t^2+0)-2.0*(s^2.0)*t-s*(t^2.0)+2.0*(r^2.0)+4.0*r*s+3.0*r*t+2.0*(s^2.0)+3.0*s*t+(t^2.0)-2.0*r-2.0*s-t)
    N[2]  = 0.5*(-2.0*(r^2.0)*t+r*(t^2.0)+2.0*(r^2.0)+r*t-2*r)
    N[3]  = 0.5*(-2.0*(s^2.0)*t+s*(t^2.0)+2.0*(s^2.0)+s*t-2.0*s)
    N[4]  = 0.5*(2.0*(r^2.0)*t+4.0*r*s*t-r*(t^2.0)+2.0*(s^2.0)*t-s*(t^2.0)+2.0*(r^2.0)+4.0*r*s-3.0*r*t+2.0*(s^2.0)-3.0*s*t+(t^2.0)-2.0*r-2.0*s+t)
    N[5]  = 0.5*(2.0*(r^2.0)*t+r*(t^2.0)+2.0*(r^2.0)-r*t-2.0*r)
    N[6]  = 0.5*(2.0*(s^2.0)*t+s*(t^2.0)+2.0*(s^2.0)-s*t-2.0*s)
    N[7]  = (2.0*(r^2.0)*t+2*r*s*t-2.0*(r^2.0)-2.0*r*s-2.0*r*t+2*r)
    N[8]  = (-2.0*r*s*t+2.0*r*s)
    N[9]  = (2.0*r*s*t+2.0*(s^2.0)*t-2.0*r*s-2.0*(s^2.0)-2.0*s*t+2.0*s)
    N[10] = (-2.0*(r^2.0)*t-2.0*r*s*t-2.0*(r^2.0)-2.0*r*s+2.0*r*t+2.0*r)
    N[11] = (2.0*r*s*t+2.0*r*s)
    N[12] = (-2.0*r*s*t-2.0*(s^2.0)*t-2.0*r*s-2*(s^2.0)+2.0*s*t+2.0*s)
    N[13] = (r*(t^2.0)+s*(t^2.0)-(t^2.0)-r-s+1.0)
    N[14] = (-r*(t^2.0)+r)
    N[15] = (-s*(t^2.0)+s)
    return N
end

function shape_deriv_WED15(R::AbstractArray{<:Float64,1})
    r, s, t = R
    D = Array{Float64}(undef,15,3)
    # Derivatives with respect to r
    D[ 1,1] = 0.5*(-4.0*r*t-4.0*s*t-(t^2.0)+4.0*r+4.0*s+3.0*t-2.0)
    D[ 2,1] = 0.5*(-4.0*r*t+(t^2.0)+4.0*r+t-2)
    D[ 3,1] = 0.0
    D[ 4,1] = 0.5*(4.0*r*t+4.0*s*t-(t^2.0)+4.0*r+4.0*s-3.0*t-2.0)
    D[ 5,1] = 0.5*(4.0*r*t+(t^2.0)+4.0*r-t-2.0)
    D[ 6,1] = 0.0
    D[ 7,1] = (4.0*r*t+2.0*s*t-4.0*r-2.0*s-2.0*t+2.0)
    D[ 8,1] = (-2.0*s*t+2.0*s)
    D[ 9,1] = (2.0*s*t-2.0*s)
    D[10,1] = (-4.0*r*t-2.0*s*t-4.0*r-2.0*s+2.0*t+2.0)
    D[11,1] = (2.0*s*t+2.0*s)
    D[12,1] = (-2.0*s*t-2.0*s)
    D[13,1] = ((t^2.0)-1.0)
    D[14,1] = ((-t^2.0)+1.0)
    D[15,1] = 0.0

    # Derivatives with respect to s
    D[ 1,2] = 0.5*(-4.0*r*t-4.0*s*t-(t^2.0)+4.0*r+4.0*s+3.0*t-2)
    D[ 2,2] = 0.0
    D[ 3,2] = 0.5*(-4.0*s*t+(t^2.0)+4*s+t-2.0)
    D[ 4,2] = 0.5*(4.0*r*t+4.0*s*t-(t^2.0)+4.0*r+4.0*s-3.0*t-2.0)
    D[ 5,2] = 0.0
    D[ 6,2] = 0.5*(4.0*s*t+(t^2.0)+4.0*s-t-2.0)
    D[ 7,2] = (2.0*r*t-2.0*r)
    D[ 8,2] = (-2.0*r*t+2.0*r)
    D[ 9,2] = (2.0*r*t+4.0*s*t-2.0*r-4.0*s-2.0*t+2.0)
    D[10,2] = (-2.0*r*t-2.0*r)
    D[11,2] = (2.0*r*t+2.0*r)
    D[12,2] = (-2.0*r*t-4.0*s*t-2.0*r-4.0*s+2.0*t+2.0)
    D[13,2] = ((t^2.0)-1.0)
    D[14,2] = 0.0
    D[15,2] = ((-t^2.0)+1.0)


    # Derivatives with respect to t
    D[ 1,3] = 0.5*(-2.0*(r^2.0)-4.0*r*s-2.0*r*t-2.0*(s^2.0)-2.0*s*t+3.0*r+3.0*s+2.0*t-1.0)
    D[ 2,3] = 0.5*(-2.0*(r^2.0)+2.0*r*t+r)
    D[ 3,3] = 0.5*(-2.0*(s^2.0)+2*s*t+s)
    D[ 4,3] = 0.5*(2.0*(r^2.0)+4.0*r*s-2.0*r*t+2.0*(s^2.0)-2.0*s*t-3.0*r-3.0*s+2.0*t+1)
    D[ 5,3] = 0.5*(2.0*(r^2)+2.0*r*t-r)
    D[ 6,3] = 0.5*(2*(s^2.0)+2*s*t-s)
    D[ 7,3] = (2.0*(r^2.0)+2.0*r*s-2*r)
    D[ 8,3] = (-2.0*r*s)
    D[ 9,3] = (2.0*r*s+2.0*(s^2.0)-2.0*s)
    D[10,3] = (2.0*r-2.0*(r^2.0)-2.0*r*s)
    D[11,3] = (2.0*r*s)
    D[12,3] = (2.0*s-2.0*s*r-2.0*(s^2.0))
    D[13,3] = (-2.0*t+2.0*t*r+2.0*t*s)
    D[14,3] = (-2.0*r*t)
    D[15,3] = (-2.0*s*t)

    return D
end

# constructor
function MakeWED15()
    shape             = CellShape()
    shape.name        = "WED15"
    shape.family      = SOLID_CELL
    shape.ndim        = 3
    shape.npoints     = 15
    shape.basic_shape = WED6
    shape.vtk_type    = VTK_QUADRATIC_WEDGE
    shape.facet_idxs  = facet_idxs_WED15
    shape.edge_idxs   = edge_idxs_WED15
    shape.facet_shape = (QUAD8, QUAD8, QUAD8, TRI6, TRI6, TRI6)
    shape.nat_coords  = coords_WED15
    shape.quadrature  = Dict( 0 => WED_IP9, 2 => WED_IP2, 9 => WED_IP9 )
    shape.func        = shape_func_WED15
    shape.deriv       = shape_deriv_WED15
    return shape
end


# Registration
const  WED15 = MakeWED15()
export WED15
