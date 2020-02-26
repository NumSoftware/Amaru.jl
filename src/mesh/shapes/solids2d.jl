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


# TRI3 shape
# ==========

#    s
#    ^
#    |
#  3
#    @,(0,1)
#    | ',
#    |   ',
#    |     ',
#    |       ',
#    |         ',
#    |           ',
#    |             ',
#    |               ',
#    |(0,0)            ', (1,0)
#    @-------------------@  --> r
#  1                      2
#

# natural coordinates
const coords_TRI3 =
[ 0.0  0.0
  1.0  0.0
  0.0  1.0 ]

const facet_idxs_TRI3 = [ [1, 2], [2, 3], [3, 1] ]

function shape_func_TRI3(R::Array{Float64,1})
    r, s = R[1:2]
    N = Array{Float64}(undef,3)
    N[1] = 1.0-r-s
    N[2] = r
    N[3] = s
    return N
end

const deriv_TRI3_mat =
[ -1.0 1.0 0.0
  -1.0 0.0 1.0 ]

function shape_deriv_TRI3(R::Array{Float64,1})
    return deriv_TRI3_mat
end

# constructor
function MakeTRI3()
    shape             = ShapeType()
    shape.name        = "TRI3"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 3
    shape.basic_shape = shape
    shape.vtk_type    = VTK_TRIANGLE
    shape.facet_idxs  = facet_idxs_TRI3
    shape.edge_idxs   = facet_idxs_TRI3
    shape.facet_shape = LIN2
    shape.nat_coords  = coords_TRI3
    shape.quadrature  = Dict( 0 => TRI_IP3,  1 => TRI_IP1,  3 => TRI_IP3,   6 => TRI_IP6 )
    shape.func        = shape_func_TRI3
    shape.deriv       = shape_deriv_TRI3
    return shape
end


# Registration
const  TRI3 = MakeTRI3()
export TRI3



# TRI6 shape
# ==========

#    s
#    ^
#    |
#  2
#    @,(0,1)
#    | ',
#    |   ',
#    |     ',
#    |       ',  4
#  5 @         '@
#    |           ',
#    |             ',
#    |               ',
#    |(0,0)            ', (1,0)
#    @---------@---------@  --> r
#  0           3          1
#

# natural coordinates
const coords_TRI6 =
[ 0.0  0.0
  1.0  0.0
  0.0  1.0
  0.5  0.0
  0.5  0.5
  0.0  0.5 ]

const facet_idxs_TRI6 = [ [1, 2, 4], [2, 3, 5], [3, 1, 6] ]

function shape_func_TRI6(R::Array{Float64,1})
    r, s = R

    N = Array{Float64}(undef,6)
    N[1] = 1.0-(r+s)*(3.0-2.0*(r+s))
    N[2] = r*(2.0*r-1.0)
    N[3] = s*(2.0*s-1.0)
    N[4] = 4.0*r*(1.0-(r+s))
    N[5] = 4.0*r*s
    N[6] = 4.0*s*(1.0-(r+s))

    return N
end

function shape_deriv_TRI6(R::Array{Float64,1})
    r, s = R

    D = Array{Float64}(undef,2, 6)
    D[1,1] = -3.0 + 4.0 * (r + s);       D[2,1] = -3.0 + 4.0*(r + s)
    D[1,2] =  4.0 * r - 1.0;              D[2,2] =  0.0
    D[1,3] =  0.0;                       D[2,3] =  4.0 * s - 1.0
    D[1,4] =  4.0 - 8.0 * r - 4.0 * s;   D[2,4] = -4.0 * r
    D[1,5] =  4.0 * s;                   D[2,5] =  4.0 * r
    D[1,6] = -4.0 * s;                   D[2,6] =  4.0 - 4.0 * r - 8.0*s

    return D
end

# constructor
function MakeTRI6()
    shape             = ShapeType()
    shape.name        = "TRI6"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 6
    shape.basic_shape = TRI3
    shape.vtk_type    = VTK_QUADRATIC_TRIANGLE
    shape.facet_idxs  = facet_idxs_TRI6
    shape.edge_idxs   = []
    shape.facet_shape = LIN3
    shape.nat_coords  = coords_TRI6
    shape.quadrature  = Dict( 0 => TRI_IP3,  3 => TRI_IP3,  6 => TRI_IP6 )
    shape.func        = shape_func_TRI6
    shape.deriv       = shape_deriv_TRI6
    return shape
end


# Registration
const  TRI6 = MakeTRI6()
export TRI6




# TRI9 shape
# ==========


# natural coordinates
const coords_TRI9 = []

const facet_idxs_TRI9 = [ [1, 2, 4, 7], [2, 3, 5, 8], [3, 1, 6, 9] ]

function shape_func_TRI9(R::Array{Float64,1})
    error("TRI9 shape not fully implemented")
    r, s = R

    N = Array{Float64}(undef,9)
    return N
end

function shape_deriv_TRI9(R::Array{Float64,1})
    error("TRI9 shape not fully implemented")
    r, s = R

    D = Array{Float64}(undef,2, 9)
    return D
end

# constructor
function MakeTRI9()
    shape             = ShapeType()
    shape.name        = "TRI9"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 9
    shape.basic_shape = TRI3
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = facet_idxs_TRI9
    shape.edge_idxs   = facet_idxs_TRI9
    shape.facet_shape = LIN4
    shape.nat_coords  = coords_TRI9
    shape.quadrature  = Dict( 0 => TRI_IP6, 3 => TRI_IP3, 6 => TRI_IP6 )
    shape.func        = shape_func_TRI9
    shape.deriv       = shape_deriv_TRI9
    return shape
end


# Registration
const  TRI9 = MakeTRI9()
export TRI9




# TRI10 shape
# ===========


# natural coordinates
const coords_TRI10 = []

const facet_idxs_TRI10 = []

function shape_func_TRI10(R::Array{Float64,1})
    error("TRI10 shape not fully implemented")
    r, s = R

    N = Array{Float64}(undef,10)
    return N
end

function shape_deriv_TRI10(R::Array{Float64,1})
    error("TRI10 shape not fully implemented")
    r, s = R

    D = Array{Float64}(undef,2, 10)
    return D
end

# constructor
function MakeTRI10()
    shape             = ShapeType()
    shape.name        = "TRI10"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 10
    shape.basic_shape = TRI3
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = facet_idxs_TRI10
    shape.edge_idxs   = facet_idxs_TRI10
    shape.facet_shape = LIN4
    shape.nat_coords  = coords_TRI10
    shape.quadrature  = Dict( 0 => TRI_IP6, 3 => TRI_IP3, 6 => TRI_IP6 )
    shape.func        = shape_func_TRI10
    shape.deriv       = shape_deriv_TRI10
    return shape
end


# Registration
const  TRI10 = MakeTRI10()
export TRI10


# QUAD4 shape
# ===========

#     4                        3
#       @--------------------@
#       |               (1,1)|
#       |       s ^          |
#       |         |          |
#       |         |          |
#       |         +----> r   |
#       |       (0,0)        |
#       |                    |
#       |                    |
#       |(-1,-1)             |
#       @--------------------@
#     1                        2
#

# natural coordinates
const coords_QUAD4 =
[ -1.0 -1.0
   1.0 -1.0
   1.0  1.0
  -1.0  1.0  ]

const facet_idxs_QUAD4 = [ [1, 2], [2, 3], [3, 4], [4, 1] ]

function shape_func_QUAD4(R::Array{Float64,1})
    r, s = R[1:2]
    N = Array{Float64}(undef,4)
    N[1] = 0.25*(1.0-r-s+r*s)
    N[2] = 0.25*(1.0+r-s-r*s)
    N[3] = 0.25*(1.0+r+s+r*s)
    N[4] = 0.25*(1.0-r+s-r*s)
    return N
end

function shape_deriv_QUAD4(R::Array{Float64,1})
    r, s = R[1:2]
    D = Array{Float64}(undef,2, 4)
    D[1,1] = 0.25*(-1.0+s);   D[2,1] = 0.25*(-1.0+r)
    D[1,2] = 0.25*(+1.0-s);   D[2,2] = 0.25*(-1.0-r)
    D[1,3] = 0.25*(+1.0+s);   D[2,3] = 0.25*(+1.0+r)
    D[1,4] = 0.25*(-1.0-s);   D[2,4] = 0.25*(+1.0-r)
    return D
end

# constructor
function MakeQUAD4()
    shape             = ShapeType()
    shape.name        = "QUAD4"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 4
    shape.basic_shape = shape
    shape.vtk_type    = VTK_QUAD
    shape.facet_idxs  = facet_idxs_QUAD4
    shape.edge_idxs   = []
    shape.facet_shape = LIN2
    shape.nat_coords  = coords_QUAD4
    shape.quadrature  = Dict( 0 => QUAD_IP2, 1 => QUAD_IP1, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4 )
    shape.func        = shape_func_QUAD4
    shape.deriv       = shape_deriv_QUAD4
    return shape
end


# Registration
const QUAD4 = MakeQUAD4()
export QUAD4



# QUAD8 shape
# ===========

#     4           7            3
#       @---------@----------@
#       |               (1,1)|
#       |       s ^          |
#       |         |          |
#       |         |          |
#     8 @         +----> r   @ 6
#       |       (0,0)        |
#       |                    |
#       |                    |
#       |(-1,-1)             |
#       @---------@----------@
#     1           5            2
#

# natural coordinates
const coords_QUAD8 =
[ -1.0 -1.0
   1.0 -1.0
   1.0  1.0
  -1.0  1.0
   0.0 -1.0
   1.0  0.0
   0.0  1.0
  -1.0  0.0 ]

const facet_idxs_QUAD8 = [ [1, 2, 5], [2, 3, 6], [3, 4, 7], [4, 1, 8] ]

function shape_func_QUAD8(R::Array{Float64,1})
    r, s = R[1:2]
    N = Array{Float64}(undef,8)
    rp1=1.0+r; rm1=1.0-r;
    sp1=1.0+s; sm1=1.0-s;
    N[1] = 0.25*rm1*sm1*(rm1+sm1-3.0)
    N[2] = 0.25*rp1*sm1*(rp1+sm1-3.0)
    N[3] = 0.25*rp1*sp1*(rp1+sp1-3.0)
    N[4] = 0.25*rm1*sp1*(rm1+sp1-3.0)
    N[5] = 0.50*sm1*(1.0-r*r)
    N[6] = 0.50*rp1*(1.0-s*s)
    N[7] = 0.50*sp1*(1.0-r*r)
    N[8] = 0.50*rm1*(1.0-s*s)
    return N
end

function shape_deriv_QUAD8(R::Array{Float64,1})
    r, s = R[1:2]
    D = Array{Float64}(undef,2, 8)
    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s

    D[1,1] = -0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0)
    D[1,2] =  0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0)
    D[1,3] =  0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0)
    D[1,4] = -0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0)
    D[1,5] = -r * sm1
    D[1,6] =  0.50 * (1.0 - s * s)
    D[1,7] = -r * sp1
    D[1,8] = -0.5 * (1.0 - s * s)

    D[2,1] = -0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0)
    D[2,2] = -0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0)
    D[2,3] =  0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0)
    D[2,4] =  0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0)
    D[2,5] = -0.50 * (1.0 - r * r)
    D[2,6] = -s * rp1
    D[2,7] =  0.50 * (1.0 - r * r)
    D[2,8] = -s * rm1
    return D
end

# constructor
function MakeQUAD8()
    shape             = ShapeType()
    shape.name        = "QUAD8"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 8
    shape.basic_shape = QUAD4
    shape.vtk_type    = VTK_QUADRATIC_QUAD
    shape.facet_idxs  = facet_idxs_QUAD8
    shape.edge_idxs   = []
    shape.facet_shape = LIN3
    shape.nat_coords  = coords_QUAD8
    shape.quadrature  = Dict( 0 => QUAD_IP2, 1 => QUAD_IP1, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4 )
    shape.func        = shape_func_QUAD8
    shape.deriv       = shape_deriv_QUAD8
    return shape
end


# Registration
const QUAD8 = MakeQUAD8()
export QUAD8


# QUAD9 shape
# ===========

#     4           7            3
#       @---------@----------@
#       |               (1,1)|
#       |       s ^          |
#       |         |          |
#       |         | 9        |
#     8 @         +----> r   @ 6
#       |       (0,0)        |
#       |                    |
#       |                    |
#       |(-1,-1)             |
#       @---------@----------@
#     1           5            2
#

# natural coordinates
const coords_QUAD9 =
[ -1.0 -1.0
   1.0 -1.0
   1.0  1.0
  -1.0  1.0
   0.0 -1.0
   1.0  0.0
   0.0  1.0
  -1.0  0.0
   0.0  0.0 ]

const facet_idxs_QUAD9 = [ [1, 2, 5], [2, 3, 6], [3, 4, 7], [4, 1, 8] ]

function shape_func_QUAD9(R::Array{Float64,1})
    r, s = R[1:2]
    N = Array{Float64}(undef,9)
    rp1=r+1.0; rm1=r-1.0
    sp1=s+1.0; sm1=s-1.0

    N[1] =  0.25*r*s*sm1*rm1
    N[2] =  0.25*r*s*sm1*rp1
    N[3] =  0.25*r*s*sp1*rp1
    N[4] =  0.25*r*s*sp1*rm1
    N[5] = -0.50*s*sm1*rp1*rm1
    N[6] = -0.50*r*sp1*sm1*rp1
    N[7] = -0.50*s*sp1*rp1*rm1
    N[8] = -0.50*r*sp1*sm1*rm1
    N[9] = sp1*sm1*rp1*rm1
    return N
end

function shape_deriv_QUAD9(R::Array{Float64,1})
    r, s = R[1:2]
    D = Array{Float64}(undef,2,9)
    rp1=r+1.0; rm1=r-1.0
    sp1=s+1.0; sm1=s-1.0

    D[1,1] = (r + rm1)*s*sm1/4.0
    D[1,2] = (r + rp1)*s*sm1/4.0
    D[1,3] = (r + rp1)*s*sp1/4.0
    D[1,4] = (r + rm1)*s*sp1/4.0
    D[1,5] = -(r + r)*s*(sm1)/2.0
    D[1,6] = -(r + rp1)*(s*s - 1.0)/2.0
    D[1,7] = -(r + r)*s*(sp1)/2.0
    D[1,8] = -(r + rm1)*(s*s - 1.0)/2.0
    D[1,9] = 2.0*r*(s*s - 1.0)

    D[2,1] = r*rm1*(s+sm1)/4.0
    D[2,2] = r*rp1*(s+sm1)/4.0
    D[2,3] = r*rp1*(s+sp1)/4.0
    D[2,4] = r*rm1*(s+sp1)/4.0
    D[2,5] = -(r*r - 1.0)*(s + sm1)/2.0
    D[2,6] = -r*rp1*(s + s)/2.0
    D[2,7] = -(r*r - 1.0)*(s + sp1)/2.0
    D[2,8] = -r*rm1*(s + s)/2.0
    D[2,9] = 2.0*s*(r*r - 1.0)

    return D
end

# constructor
function MakeQUAD9()
    shape             = ShapeType()
    shape.name        = "QUAD9"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 9
    shape.basic_shape = QUAD4
    shape.vtk_type    = VTK_BIQUADRATIC_QUAD
    shape.facet_idxs  = facet_idxs_QUAD9
    shape.edge_idxs   = []
    shape.facet_shape = LIN3
    shape.nat_coords  = coords_QUAD9
    shape.quadrature  = Dict( 0 => QUAD_IP3, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4 )
    shape.func        = shape_func_QUAD9
    shape.deriv       = shape_deriv_QUAD9
    return shape
end

# Registration
const QUAD9 = MakeQUAD9()
export QUAD9



# QUAD12 shape
# ============

#
#    4      11       7        3
#      @-----@-------@------@
#      |               (1,1)|
#      |       s ^          |
#    8 @         |          @ 10
#      |         |          |
#      |         +----> r   |
#      |       (0,0)        |
#   12 @                    @ 6
#      |                    |
#      |(-1,-1)             |
#      @-----@-------@------@
#    1       5       9        2
#

# natural coordinates
const coords_QUAD12 =
[   -1.0       -1.0
     1.0       -1.0
     1.0        1.0
    -1.0        1.0
    -1.0/3.0   -1.0
     1.0       -1.0/3.0
     1.0/3.0    1.0
    -1.0        1.0/3.0
     1.0/3.0   -1.0
     1.0        1.0/3.0
    -1.0/3.0    1.0
    -1         -1.0/3.0 ]

const facet_idxs_QUAD12 = [ [1, 2, 5, 9], [2, 3, 6, 10], [3, 4, 7, 11], [4, 1, 8, 12] ]

function shape_func_QUAD12(R::Array{Float64,1})
    r, s = R[1:2]
    N = Array{Float64}(undef,12)

    RM = 1.0 - r
    RP = 1.0 + r
    SM = 1.0 - s
    SP = 1.0 + s
    N[1]  = RM*SM*( 9.0*(r*r + s*s) - 10.0)/32.0
    N[2]  = RP*SM*( 9.0*(r*r + s*s) - 10.0)/32.0
    N[3]  = RP*SP*( 9.0*(r*r + s*s) - 10.0)/32.0
    N[4]  = RM*SP*( 9.0*(r*r + s*s) - 10.0)/32.0
    N[5]  = 9.0*(1.0 - r*r)*(1.0 - 3.0*r)*SM/32.0
    N[6]  = 9.0*(1.0 - s*s)*(1.0 - 3.0*s)*RP/32.0
    N[7]  = 9.0*(1.0 - r*r)*(1.0 + 3.0*r)*SP/32.0
    N[8]  = 9.0*(1.0 - s*s)*(1.0 + 3.0*s)*RM/32.0
    N[9]  = 9.0*(1.0 - r*r)*(1.0 + 3.0*r)*SM/32.0
    N[10] = 9.0*(1.0 - s*s)*(1.0 + 3.0*s)*RP/32.0
    N[11] = 9.0*(1.0 - r*r)*(1.0 - 3.0*r)*SP/32.0
    N[12] = 9.0*(1.0 - s*s)*(1.0 - 3.0*s)*RM/32.0

    return N
end

function shape_deriv_QUAD12(R::Array{Float64,1})
    r, s = R[1:2]
    D = Array{Float64}(undef,2, 12)

    RP = 1.0 + r
    RM = 1.0 - r
    SP = 1.0 + s
    SM = 1.0 - s

    D[1,1]  =  SM*(9.0*(2.0*r - 3.0*r*r - s*s) + 10.0)/32.0
    D[1,2]  =  SM*(9.0*(2.0*r + 3.0*r*r + s*s) - 10.0)/32.0
    D[1,3]  =  SP*(9.0*(2.0*r + 3.0*r*r + s*s) - 10.0)/32.0
    D[1,4]  =  SP*(9.0*(2.0*r - 3.0*r*r - s*s) + 10.0)/32.0
    D[1,5]  =  9.0*SM*(9.0*r*r - 2.0*r - 3.0)/32.0
    D[1,6]  =  9.0*(1.0 - s*s)*(1.0 - 3.0*s)/32.0
    D[1,7]  =  9.0*SP*(-9.0*r*r - 2.0*r + 3.0)/32.0
    D[1,8]  = -9.0*(1.0 - s*s)*(1.0 + 3.0*s)/32.0
    D[1,9]  =  9.0*SM*(-9.0*r*r - 2.0*r + 3.0)/32.0
    D[1,10] =  9.0*(1.0 - s*s)*(1.0 + 3.0*s)/32.0
    D[1,11] =  9.0*SP*(9.0*r*r - 2.0*r - 3.0)/32.0
    D[1,12] = -9.0*(1.0 - s*s)*(1.0 - 3.0*s)/32.0
    D[2,1]  =  RM*(9.0*(2.0*s - 3.0*s*s - r*r) + 10.0)/32.0
    D[2,2]  =  RP*(9.0*(2.0*s - 3.0*s*s - r*r) + 10.0)/32.0
    D[2,3]  =  RP*(9.0*(2.0*s + 3.0*s*s + r*r) - 10.0)/32.0
    D[2,4]  =  RM*(9.0*(2.0*s + 3.0*s*s + r*r) - 10.0)/32.0
    D[2,5]  = -9.0*(1.0 - r*r)*(1.0 - 3.0*r)/32.0
    D[2,6]  =  9.0*RP*(9.0*s*s - 2.0*s - 3.0)/32.0
    D[2,7]  =  9.0*(1.0 - r*r)*(1.0 + 3.0*r)/32.0
    D[2,8]  =  9.0*RM*(-9.0*s*s - 2.0*s + 3.0)/32.0
    D[2,9]  = -9.0*(1.0 - r*r)*(1.0 + 3.0*r)/32.0
    D[2,10] =  9.0*RP*(-9.0*s*s - 2.0*s + 3.0)/32.0
    D[2,11] =  9.0*(1.0 - r*r)*(1.0 - 3.0*r)/32.0
    D[2,12] =  9.0*RM*(9.0*s*s - 2.0*s - 3.0)/32.0

    return D
end

# constructor
function MakeQUAD12()
    shape             = ShapeType()
    shape.name        = "QUAD12"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 12
    shape.basic_shape = QUAD4
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = facet_idxs_QUAD12
    shape.edge_idxs   = []
    shape.facet_shape = LIN4
    shape.nat_coords  = coords_QUAD12
    shape.quadrature  = Dict( 0 => QUAD_IP3, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4 )
    shape.func        = shape_func_QUAD12
    shape.deriv       = shape_deriv_QUAD12
    return shape
end

# Registration
const QUAD12 = MakeQUAD12()
export QUAD12



# QUAD16 shape
# ============


# natural coordinates
const coords_QUAD16 = []

const facet_idxs_QUAD16 = [ [1, 2, 5, 9], [2, 3, 6, 10], [3, 4, 7, 11], [4, 1, 8, 12] ]

function shape_func_QUAD16(R::Array{Float64,1})
    error("TRI9 shape not fully implemented")
    r, s = R[1:2]
    N = Array{Float64}(undef,16)

    return N
end

function shape_deriv_QUAD16(R::Array{Float64,1})
    error("TRI9 shape not fully implemented")
    r, s = R[1:2]
    D = Array{Float64}(undef,2, 16)

    return D
end

# constructor
function MakeQUAD16()
    shape             = ShapeType()
    shape.name        = "QUAD16"
    shape.family      = SOLID_SHAPE
    shape.ndim        = 2
    shape.npoints     = 16
    shape.basic_shape = QUAD4
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = facet_idxs_QUAD16
    shape.edge_idxs   = []
    shape.facet_shape = LIN4
    shape.nat_coords  = coords_QUAD16
    shape.quadrature  = Dict( 0 => QUAD_IP4, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4 )
    shape.func        = shape_func_QUAD16
    shape.deriv       = shape_deriv_QUAD16
    return shape
end

# Registration
const QUAD16 = MakeQUAD16()
export QUAD16


