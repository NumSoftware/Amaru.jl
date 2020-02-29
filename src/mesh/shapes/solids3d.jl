##############################################################################
#    FemLab - Finite Element Library                                         #
#    Copyright (C) Raul Durand <raul.durand at gmail.com>                    #
#    All rights reserved.                                                    #
#                                                                            #
#    This software is distributed WITHOUT ANY WARRANTY; without even         #
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR     #
#    PURPOSE. See the GNU General Public License for more details.           #
#                                                                            #
##############################################################################




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

function shape_func_TET4(R::Array{Float64,1})
    r, s, t = R

    N = Array{Float64}(undef,4)
    N[1] = 1.0-r-s-t
    N[2] = r
    N[3] = s
    N[4] = t
    return N
end

function shape_deriv_TET4(R::Array{Float64,1})
    r, s, t = R

    D = Array{Float64}(undef,3, 4)
    D[1,1] = -1.0;   D[2,1]= -1.0;   D[3,1]= -1.0
    D[1,2] =  1.0;   D[2,2]=  0.0;   D[3,2]=  0.0
    D[1,3] =  0.0;   D[2,3]=  1.0;   D[3,3]=  0.0
    D[1,4] =  0.0;   D[2,4]=  0.0;   D[3,4]=  1.0
    return D
end

# constructor
function MakeTET4()
    shape             = ShapeType()
    shape.name        = "TET4"
    shape.family      = SOLID_SHAPE
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
[  0.0      0.0   0.0
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

function shape_func_TET10(R::Array{Float64,1})
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

function shape_deriv_TET10(R::Array{Float64,1})
    r, s, t = R

    D = Array{Float64}(undef,3, 10)

    # r-derivatives: dN0/dr to dN9/dr
    D[1,1]  =  4.0*(r + s + t) - 3.0
    D[1,2]  =  4.0*r - 1.0
    D[1,3]  =  0.0
    D[1,4]  =  0.0
    D[1,5]  =  4.0 - 8.0*r - 4.0*s - 4.0*t
    D[1,6]  =  4.0*s
    D[1,7]  = -4.0*s
    D[1,8]  = -4.0*t
    D[1,9]  =  4.0*t
    D[1,10] =  0.0

    # s-derivatives: dN0/ds to dN9/ds
    D[2,1]  =  4.0*(r + s + t) - 3.0
    D[2,2]  =  0.0
    D[2,3]  =  4.0*s - 1.0
    D[2,4]  =  0.0
    D[2,5]  = -4.0*r
    D[2,6]  =  4.0*r
    D[2,7]  =  4.0 - 4.0*r - 8.0*s - 4.0*t
    D[2,8]  = -4.0*t
    D[2,9]  =  0.0
    D[2,10] =  4.0*t

    # t-derivatives: dN0/dt to dN9/dt
    D[3,1]  =  4.0*(r + s + t) - 3.0
    D[3,2]  =  0.0
    D[3,3]  =  0.0
    D[3,4]  =  4.0*t - 1.0
    D[3,5]  = -4.0*r
    D[3,6]  =  0.0
    D[3,7]  = -4.0*s
    D[3,8]  =  4.0 - 4.0*r - 4.0*s - 8.0*t
    D[3,9]  =  4.0*r
    D[3,10] =  4.0*s

    return D
end

# constructor
function MakeTET10()
    shape             = ShapeType()
    shape.name        = "TET10"
    shape.family      = SOLID_SHAPE
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

function shape_func_PYR5(R::Array{Float64,1})
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

function shape_deriv_PYR5(R::Array{Float64,1})
    r, s, t = R
    w = t==1.0 ? 0.0 : 1/(1-t)


    D = Array{Float64}(undef,3, 5)
    D[1,1] = 0.25*(-1+s*w)  ;   D[2,1]= 0.25*(-1+r*w)  ;   D[3,1]= 0.25*(-1+r*s*w^2)
    D[1,2] = 0.25*(1-s*w)   ;   D[2,2]= 0.25*(-1-r*w)  ;   D[3,2]= 0.25*(-1-r*s*w^2)
    D[1,3] = 0.25*(1+s*w)   ;   D[2,3]= 0.25*(1+r*w)   ;   D[3,3]= 0.25*(-1+r*s*w^2)
    D[1,4] = 0.25*(-1-s*w)  ;   D[2,4]= 0.25*(1-r*w)   ;   D[3,4]= 0.25*(-1-r*s*w^2)
    D[1,5] =  0.0           ;   D[2,5]=  0.0           ;   D[3,5]=  1.0
    return D
end

# constructor
function MakePYR5()
    shape             = ShapeType()
    shape.name        = "PYR5"
    shape.family      = SOLID_SHAPE
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

function shape_func_HEX8(R::Array{Float64,1})
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

function shape_deriv_HEX8(R::Array{Float64,1})
    r, s, t = R
    st = s*t
    rt = r*t
    rs = r*s
    D = Array{Float64}(undef,3, 8)
    D[1,1] = -1.0+s+t-st;   D[2,1]=-1.0+r+t-rt;   D[3,1]=-1.0+r+s-rs
    D[1,2] = +1.0-s-t+st;   D[2,2]=-1.0-r+t+rt;   D[3,2]=-1.0-r+s+rs
    D[1,3] = +1.0+s-t-st;   D[2,3]=+1.0+r-t-rt;   D[3,3]=-1.0-r-s-rs
    D[1,4] = -1.0-s+t+st;   D[2,4]=+1.0-r-t+rt;   D[3,4]=-1.0+r-s+rs
    D[1,5] = -1.0+s-t+st;   D[2,5]=-1.0+r-t+rt;   D[3,5]=+1.0-r-s+rs
    D[1,6] = +1.0-s+t-st;   D[2,6]=-1.0-r-t-rt;   D[3,6]=+1.0+r-s-rs
    D[1,7] = +1.0+s+t+st;   D[2,7]=+1.0+r+t+rt;   D[3,7]=+1.0+r+s+rs
    D[1,8] = -1.0-s-t-st;   D[2,8]=+1.0-r+t-rt;   D[3,8]=+1.0-r+s-rs
    D = 0.125*D

    return D
end

# constructor
function MakeHEX8()
    shape             = ShapeType()
    shape.name        = "HEX8"
    shape.family      = SOLID_SHAPE
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

function shape_func_HEX20(R::Array{Float64,1})
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

function shape_deriv_HEX20(R::Array{Float64,1})
    r, s, t = R

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    D = Array{Float64}(undef,3, 20)
    # Derivatives with respect to r
    D[1, 1] = -0.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[1, 2] =  0.125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1
    D[1, 3] =  0.125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[1, 4] = -0.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[1, 5] = -0.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[1, 6] =  0.125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[1, 7] =  0.125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[1, 8] = -0.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1
    D[1, 9] = -0.5*r*sm1*tm1
    D[1,10] =  0.25*(1-s*s)*tm1
    D[1,11] = -0.5*r*sp1*tm1
    D[1,12] = -0.25*(1-s*s)*tm1
    D[1,13] = -0.5*r*sm1*tp1
    D[1,14] =  0.25*(1-s*s)*tp1
    D[1,15] = -0.5*r*sp1  *tp1
    D[1,16] = -0.25*(1-s*s)*tp1
    D[1,17] = -0.25*sm1*(1-t*t)
    D[1,18] =  0.25*sm1*(1-t*t)
    D[1,19] =  0.25*sp1*(1-t*t)
    D[1,20] = -0.25*sp1*(1-t*t)

    # Derivatives with respect to s
    D[2, 1] = -0.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[2, 2] = -0.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[2, 3] =  0.125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[2, 4] =  0.125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1
    D[2, 5] = -0.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[2, 6] = -0.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1
    D[2, 7] =  0.125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[2, 8] =  0.125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[2, 9] = -0.25*(1-r*r)*tm1
    D[2,10] = -0.5*s*rp1*tm1
    D[2,11] =  0.25*(1-r*r)*tm1
    D[2,12] = -0.5*s*rm1*tm1
    D[2,13] = -0.25*(1-r*r)*tp1
    D[2,14] = -0.5*s*rp1*tp1
    D[2,15] =  0.25*(1-r*r)*tp1
    D[2,16] = -0.5*s*rm1*tp1
    D[2,17] = -0.25*rm1*(1-t*t)
    D[2,18] = -0.25*rp1*(1-t*t)
    D[2,19] =  0.25*rp1*(1-t*t)
    D[2,20] =  0.25*rm1*(1-t*t)

    # Derivatives with respect to t
    D[3, 1] = -0.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[3, 2] = -0.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[3, 3] = -0.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1
    D[3, 4] = -0.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[3, 5] =  0.125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1
    D[3, 6] =  0.125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[3, 7] =  0.125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[3, 8] =  0.125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[3, 9] = -0.25*(1-r*r)*sm1
    D[3,10] = -0.25*rp1*(1-s*s)
    D[3,11] = -0.25*(1-r*r)*sp1
    D[3,12] = -0.25*rm1*(1-s*s)
    D[3,13] =  0.25*(1-r*r)*sm1
    D[3,14] =  0.25*rp1*(1-s*s)
    D[3,15] =  0.25*(1-r*r)*sp1
    D[3,16] =  0.25*rm1*(1-s*s)
    D[3,17] = -0.5*t*rm1*sm1
    D[3,18] = -0.5*t*rp1*sm1
    D[3,19] = -0.5*t*rp1*sp1
    D[3,20] = -0.5*t*rm1*sp1

    return D
end

# constructor
function MakeHEX20()
    shape             = ShapeType()
    shape.name        = "HEX20"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 3
    shape.npoints     = 20
    shape.basic_shape = HEX8
    shape.vtk_type    = VTK_QUADRATIC_HEXAHEDRON
    shape.facet_idxs  = facet_idxs_HEX20
    shape.edge_idxs   = edge_idxs_HEX20
    shape.facet_shape = QUAD8
    shape.nat_coords  = coords_HEX20
    shape.quadrature  = Dict( 0 => HEX_IP3, 8 => HEX_IP2, 27 => HEX_IP3 )
    shape.func        = shape_func_HEX20
    shape.deriv       = shape_deriv_HEX20
    return shape
end


# Registration
const  HEX20 = MakeHEX20()
export HEX20



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

function shape_func_WED6(R::Array{Float64,1})
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

function shape_deriv_WED6(R::Array{Float64,1})
    r, s, t = R
    D = Array{Float64}(undef,3, 6)
    D[1,1] = 0.5*(-1.0+t);  D[2,1] = 0.5*(-1.0+t);  D[3,1] = 0.5*(-1.0+r+s)
    D[1,2] = 0.5*(1.0-t) ;  D[2,2] = 0.0         ;  D[3,2] = 0.5*(-r)
    D[1,3] = 0.0         ;  D[2,3] = 0.5*(1.0-t) ;  D[3,3] = 0.5*(-s)
    D[1,4] = 0.5*(-1.0-t);  D[2,4] = 0.5*(-1.0-t);  D[3,4] = 0.5*(1.0-r-s)
    D[1,5] = 0.5*(1.0+t) ;  D[2,5] = 0.0         ;  D[3,5] = 0.5*(r)
    D[1,6] = 0.0         ;  D[2,6] = 0.5*(1.0+t) ;  D[3,6] = 0.5*(s)

    return D
end

# constructor
function MakeWED6()
    shape             = ShapeType()
    shape.name        = "WED6"
    shape.family      = SOLID_SHAPE
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

function shape_func_WED15(R::Array{Float64,1})
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

function shape_deriv_WED15(R::Array{Float64,1})
    r, s, t = R
    D = Array{Float64}(undef,3, 15)
    # Derivatives with respect to r
    D[1, 1] = 0.5*(-4.0*r*t-4.0*s*t-(t^2.0)+4.0*r+4.0*s+3.0*t-2.0)
    D[1, 2] = 0.5*(-4.0*r*t+(t^2.0)+4.0*r+t-2)
    D[1, 3] = 0.0
    D[1, 4] = 0.5*(4.0*r*t+4.0*s*t-(t^2.0)+4.0*r+4.0*s-3.0*t-2.0)
    D[1, 5] = 0.5*(4.0*r*t+(t^2.0)+4.0*r-t-2.0)
    D[1, 6] = 0.0
    D[1, 7] = (4.0*r*t+2.0*s*t-4.0*r-2.0*s-2.0*t+2.0)
    D[1, 8] = (-2.0*s*t+2.0*s)
    D[1, 9] = (2.0*s*t-2.0*s)
    D[1,10] = (-4.0*r*t-2.0*s*t-4.0*r-2.0*s+2.0*t+2.0)
    D[1,11] = (2.0*s*t+2.0*s)
    D[1,12] = (-2.0*s*t-2.0*s)
    D[1,13] = ((t^2.0)-1.0)
    D[1,14] = ((-t^2.0)+1.0)
    D[1,15] = 0.0

    # Derivatives with respect to s
    D[2, 1] = 0.5*(-4.0*r*t-4.0*s*t-(t^2.0)+4.0*r+4.0*s+3.0*t-2)
    D[2, 2] = 0.0
    D[2, 3] = 0.5*(-4.0*s*t+(t^2.0)+4*s+t-2.0)
    D[2, 4] = 0.5*(4.0*r*t+4.0*s*t-(t^2.0)+4.0*r+4.0*s-3.0*t-2.0)
    D[2, 5] = 0.0
    D[2, 6] = 0.5*(4.0*s*t+(t^2.0)+4.0*s-t-2.0)
    D[2, 7] = (2.0*r*t-2.0*r)
    D[2, 8] = (-2.0*r*t+2.0*r)
    D[2, 9] = (2.0*r*t+4.0*s*t-2.0*r-4.0*s-2.0*t+2.0)
    D[2,10] = (-2.0*r*t-2.0*r)
    D[2,11] = (2.0*r*t+2.0*r)
    D[2,12] = (-2.0*r*t-4.0*s*t-2.0*r-4.0*s+2.0*t+2.0)
    D[2,13] = ((t^2.0)-1.0)
    D[2,14] = 0.0
    D[2,15] = ((-t^2.0)+1.0)


    # Derivatives with respect to t
    D[3, 1] = 0.5*(-2.0*(r^2.0)-4.0*r*s-2.0*r*t-2.0*(s^2.0)-2.0*s*t+3.0*r+3.0*s+2.0*t-1.0)
    D[3, 2] = 0.5*(-2.0*(r^2.0)+2.0*r*t+r)
    D[3, 3] = 0.5*(-2.0*(s^2.0)+2*s*t+s)
    D[3, 4] = 0.5*(2.0*(r^2.0)+4.0*r*s-2.0*r*t+2.0*(s^2.0)-2.0*s*t-3.0*r-3.0*s+2.0*t+1)
    D[3, 5] = 0.5*(2.0*(r^2)+2.0*r*t-r)
    D[3, 6] = 0.5*(2*(s^2.0)+2*s*t-s)
    D[3, 7] = (2.0*(r^2.0)+2.0*r*s-2*r)
    D[3, 8] = (-2.0*r*s)
    D[3, 9] = (2.0*r*s+2.0*(s^2.0)-2.0*s)
    D[3,10] = (2.0*r-2.0*(r^2.0)-2.0*r*s)
    D[3,11] = (2.0*r*s)
    D[3,12] = (2.0*s-2.0*s*r-2.0*(s^2.0))
    D[3,13] = (-2.0*t+2.0*t*r+2.0*t*s)
    D[3,14] = (-2.0*r*t)
    D[3,15] = (-2.0*s*t)

    return D
end

# constructor
function MakeWED15()
    shape             = ShapeType()
    shape.name        = "WED15"
    shape.family      = SOLID_SHAPE
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
