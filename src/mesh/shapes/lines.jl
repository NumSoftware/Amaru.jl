##############################################################################
#    FemLab - Finite Element Library                                         #
#    Copyright (C) 2014 Raul Durand <raul.durand at gmail.com>               #
#                                                                            #
#    This file is part of FemLab.                                            #
#                                                                            #
#    FemLab is free software: you can redistribute it and/or modify          #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    any later version.                                                      #
#                                                                            #
#    FemLab is distributed in the hope that it will be useful,               #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with FemLab.  If not, see <http://www.gnu.org/licenses/>.         #
##############################################################################


# LIN2 shape
# ==========

# natural coordinates
const coords_LIN2 = [ -1.0, 1.0 ]

# shape functions
function shape_func_LIN2(R::Array{Float64,1})
    r = R[1]
    N = Array{Float64}(undef,2)
    N[1] = 0.5*(1-r)
    N[2] = 0.5*(1+r)
    return N
end

# shape derivatives
const deriv_LIN2_mat = [ -0.5  0.5 ]

function shape_deriv_LIN2(R::Array{Float64,1})
    return deriv_LIN2_mat
end


# constructor
function MakeLIN2()
    shape             = ShapeType()
    shape.name        = "LIN2"
    shape.family      = LINE_SHAPE
    shape.ndim        = 1
    shape.npoints     = 2
    shape.basic_shape = shape
    shape.vtk_type    = VTK_LINE
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = ()
    shape.nat_coords  = coords_LIN2
    shape.quadrature  = Dict{Int, Array}(0=>LIN_IP2, 1=>ALL_IP1, 2=>LIN_IP2, 3=>LIN_IP3, 4=>LIN_IP4)
    shape.func        = shape_func_LIN2
    shape.deriv       = shape_deriv_LIN2
    return shape
end


# Registration
const LIN2 = MakeLIN2()
export LIN2



# LIN3 shape
# ==========

# natural coordinates
const coords_LIN3 = [ -1.0, 1.0,  0.0]

# shape functions
function shape_func_LIN3(R::Array{Float64,1})
    r = R[1]
    N = Array{Float64}(undef,3)
    N[1] = 0.5*(r*r - r)
    N[2] = 0.5*(r*r + r)
    N[3] = 1.0 - r*r
    return N
end

# shape derivatives
function shape_deriv_LIN3(R::Array{Float64,1})
    r = R[1]
    D = Array{Float64}(undef,1, 3)
    D[1, 1] = r - 0.5
    D[1, 2] = r + 0.5
    D[1, 3] = -2.0*r
    return D
end

# constructor
function MakeLIN3()
    shape             = ShapeType()
    shape.name        = "LIN3"
    shape.family      = LINE_SHAPE
    shape.ndim        = 1
    shape.npoints     = 3
    shape.basic_shape = LIN2
    shape.vtk_type    = VTK_QUADRATIC_EDGE
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = ()
    shape.nat_coords  = coords_LIN3
    shape.quadrature  = Dict( 0=>LIN_IP2, 1=>ALL_IP1, 2=>LIN_IP2, 3=>LIN_IP3, 4=>LIN_IP4 )
    shape.func        = shape_func_LIN3
    shape.deriv       = shape_deriv_LIN3
    return shape
end


# Registration
const LIN3 = MakeLIN3()
export LIN3



# LIN4 shape
# ==========

# natural coordinates
const coords_LIN4 = [ -1.0,  1.0,  -1.0/3.0,  1.0/3.0 ]

# shape functions
function shape_func_LIN4(R::Array{Float64,1})
    #   (-1)            '   (+1)
    #    @------@-----@------@  --> r
    #    1      3     4      2

    r = R[1]
    N = Array{Float64}(undef,4)
    N[1] = 1.0/16.0*( -9.0*r^3 + 9.0*r*r +     r - 1.0)
    N[2] = 1.0/16.0*(  9.0*r^3 + 9.0*r*r -     r - 1.0)
    N[3] = 1.0/16.0*( 27.0*r^3 - 9.0*r*r - 27.0*r + 9.0)
    N[4] = 1.0/16.0*(-27.0*r^3 - 9.0*r*r + 27.0*r + 9.0)
    return N
end

# shape derivatives
function shape_deriv_LIN4(R::Array{Float64,1})
    r = R[1]
    D = Array{Float64}(undef,1, 4)
    D[1,1] = 1.0/16.0*( -27.0*r*r + 18.0*r + 1.0 )
    D[1,2] = 1.0/16.0*(  27.0*r*r + 18.0*r - 1.0 )
    D[1,3] = 1.0/16.0*(  81.0*r*r - 18.0*r - 27.0)
    D[1,4] = 1.0/16.0*( -81.0*r*r - 18.0*r + 27.0)
    return D
end

# constructor
function MakeLIN4()
    shape             = ShapeType()
    shape.name        = "LIN4"
    shape.family      = LINE_SHAPE
    shape.ndim        = 1
    shape.npoints     = 4
    shape.basic_shape = LIN2
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = ()
    shape.nat_coords  = coords_LIN4
    shape.quadrature  = Dict( 0=>LIN_IP2, 1=>ALL_IP1, 2=>LIN_IP2, 3=>LIN_IP3, 4=>LIN_IP4 )
    shape.func        = shape_func_LIN4
    shape.deriv       = shape_deriv_LIN4
    return shape
end


# Registration
const LIN4 = MakeLIN4()
export LIN4
