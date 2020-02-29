# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

@enum(ShapeFamily,
LINE_SHAPE    = 1,
SOLID_SHAPE   = 2,
JOINT_SHAPE   = 3,
JOINT1D_SHAPE = 4,
VERTEX_SHAPE  = 5,
EMBEDDED      = 6
)

# Export
for s in instances(ShapeFamily)
    @eval export $(Symbol(s))
end

mutable struct ShapeType
    name       ::String
    family     ::ShapeFamily
    ndim       ::Int
    npoints    ::Int
    basic_shape::ShapeType
    vtk_type   ::VTKCellType
    facet_idxs ::Array
    edge_idxs  ::Array
    facet_shape::Union{ShapeType, Tuple}
    nat_coords ::Array
    quadrature ::Dict{Int, Array}
    func       ::Function
    deriv      ::Function
    function ShapeType()
        return new()
    end
end

include("shapes/lines.jl")
include("shapes/solids2d.jl")
include("shapes/solids3d.jl")
include("shapes/joints.jl")

# Shape for unknown polyvertex
function MakePOLIV()
    shape             = ShapeType()
    shape.name        = "POLYV"
    shape.family      = VERTEX_SHAPE
    shape.ndim        = 0
    shape.npoints     = 0
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.facet_shape = ()
    shape.quadrature  = Dict( 0 => [] )
    return shape
end
const POLYV = MakePOLIV()
export POLYV


function get_ip_coords(shape::ShapeType, nips=0)
    !haskey(shape.quadrature, nips) && error("Cannot set $nips integ. points for shape $shape. Available numbers are $(collect(keys(shape.quadrature))).")
    return shape.quadrature[nips]
end


# Available shape types
const VTK2SHAPE = Dict{VTKCellType,ShapeType}(
    VTK_POLY_VERTEX          => POLYV,
    VTK_LINE                 => LIN2,
    VTK_TRIANGLE             => TRI3,
    VTK_QUAD                 => QUAD4,
    VTK_PYRAMID              => PYR5,
    VTK_TETRA                => TET4,
    VTK_HEXAHEDRON           => HEX8,
    VTK_WEDGE                => WED6,
    VTK_QUADRATIC_EDGE       => LIN3,
    VTK_QUADRATIC_TRIANGLE   => TRI6,
    VTK_QUADRATIC_QUAD       => QUAD8,
    VTK_QUADRATIC_TETRA      => TET10,
    VTK_QUADRATIC_HEXAHEDRON => HEX20,
    VTK_QUADRATIC_WEDGE      => WED15,
    VTK_BIQUADRATIC_QUAD     => QUAD9
)

# All bulk cell shapes
const ALL_SHAPES = [
    LIN2
    LIN3
    LIN4
    TRI3
    TRI6
    QUAD4
    QUAD8
    QUAD9
    QUAD12
    PYR5
    TET4
    TET10
    WED6
    WED15
    HEX8
    HEX20
]


function get_shape_from_vtk(vtk_type::VTKCellType, npoints::Int64, ndim::Int64, layers::Int64=0)::ShapeType
    # vtk_type: VTK cell code
    # npoints : total number of cell points
    # ndim    : analysis dimension
    # layers  : number of layers for joint elements

    vtk_type!=VTK_POLY_VERTEX && return VTK2SHAPE[vtk_type]

    # Check if it is a joint cell with layers
    if layers in (2,3)
        # dictionary (ndim,nfpoints) => basic_shape
        shapedict = Dict( (2,2)=>:LIN2, (2,3)=>:LIN3, (2,4)=>:LIN4, (3,3)=>:TRI3, (3,4)=>:QUAD4, (3,6)=>:TRI6, (3,8)=>:QUAD8 )
        nfpoints = div(npoints,layers)
        if haskey(shapedict, (ndim, nfpoints))
            numstr = layers==3 ? "3" : ""
            shape = Symbol("J$(numstr)$(shapedict[(ndim, nfpoints)])")
            return eval(shape)
        end
    end

    # Check for other cells
    if     npoints==2   return JLINK2
    elseif npoints==3   return JLINK3
    elseif npoints==9   return TRI9
    elseif npoints==10  return TRI10
    elseif npoints==12
        #if ndim==2 return QUAD12 end
    elseif npoints==16
        if ndim==2 return QUAD16 end
    end

    error("get_shape_from_vtk: Unknown shape for vtk_type $vtk_type and npoints $npoints with ndim $ndim")
end


function get_shape_from_ndim_npoints(npoints::Int64, ndim::Int64)::ShapeType
    # npoints : total number of cell points
    # ndim    : analysis dimension

    if ndim==3
        npoints==2  &&  return LIN2
        npoints==4  &&  return TET4
        npoints==10 &&  return TET10
        npoints==8  &&  return HEX8
        npoints==20 &&  return HEX20
    else
        npoints==2  && return LIN2
        npoints==3  && return TRI3
        npoints==6  && return TRI6
        npoints==4  && return QUAD4
        npoints==8  && return QUAD8
        npoints==12 && return QUAD12
    end

    npoints==1 && return POLYV

    error("get_shape_from_ndim_npoints: Cannot get the cell shape from npoints=$npoints and ndim=$ndim")
end


# For compatibility with gofem input file
function get_shape_from_geo(geo, npoints=0)
    types = [ LIN2, LIN3, -1, TRI3, TRI6, -1, QUAD4, QUAD8, QUAD9, PYR5, TET4, TET10, HEX8, HEX20, -2, LIN4, TRI10, QUAD12, QUAD16]
    shapetype = types[geo+1]
    if shapetype==-2 #joint
        shapetype = npoints==2 ?  JLINK2 : JLINK3
    end

    return shapetype
end


function bdistance(shape::ShapeType, R::Array{Float64,1})
    # Returns a real value which is a pseudo distance from a point to the border of an element
    # Arguments:
    #     R - a vector containing the point coordinates
    # Returns:
    #     a real value: if possitive then the point is inside the element and negative otherwise

    r, s, t = R
    bshape = shape.basic_shape
    if bshape == TRI3  return min(r, s, 1.0-r-s) end
    if bshape == QUAD4 return min(1.0 - abs(r), 1.0 - abs(s)) end
    if bshape == TET4  return min(r, s, t, 1.0-r-s-t) end
    if bshape == HEX8  return min(1.0 - abs(r), 1.0 - abs(s), 1.0 - abs(t)) end
    if bshape == WED6
        return min(r, s, 1.0-r-s, 1.0-abs(t))
    end
    error("No boundary distance for shape ($shape)")
end


function inverse_map(shape::ShapeType, coords::Array{Float64,2}, X0::Array{Float64,1}, Tol=1.0e-7)
    MAXIT = 20
    ndim  = shape.ndim
    R = zeros(ndim)
    C = coords
    local ΔX::Array{Float64,1}

    X = X0
    if size(coords,2)==2
        X = X0[1:ndim]
    end

    k = 0
    for k=1:MAXIT
        # calculate Jacobian
        D = shape.deriv(R)
        J = D*C

        # calculate trial of real coordinates
        N  = shape.func(R)
        Xt = C'*N # interpolating

        # calculate the error
        ΔX = Xt - X
        ΔR = pinv(J)'*ΔX

        # updating local coords R
        R -= ΔR
        if norm(ΔX) < Tol; break end
    end

    # TODO: Improve accuracy of inverse_map function in elements with non regular shape
    #k==MAXIT && println("Warning: max iterations (MAXIT=$MAXIT) reached in inverse mapping. norm(ΔX)=$(norm(ΔX))")

    if ndim==2
        R = vcat( R, 0.0 )
    end
    return R
end


function is_inside(shape::ShapeType, C::Array{Float64,2}, X::Array{Float64,1}, Tol = 1.e-7)
    if shape.family!=SOLID_SHAPE return false end

    # Testing with bounding box
    ndim = size(C,1)
    Cmin = vec(minimum(C,dims=1))
    Cmax = vec(maximum(C,dims=1))
    maxl = maximum(Cmax-Cmin)
    ttol = 0.1*maxl # 10% is important for curved elements

    if any(X .< Cmin .- ttol) || any(X .> Cmax .+ ttol)
        return false
    end

    # Testing with inverse mapping
    R = inverse_map(shape, C, X, Tol)
    if bdistance(shape, R) > -Tol
        return true
    else
        return false
    end
end



function extrapolator(shape::ShapeType, nips::Int)
    #  Returns a numpy matrix E that extrapolates ip values to nodal values as:
    #
    #                 NodalValues = E * IpValues;
    # where:
    #                             +                +              +
    #                        E = N * (I - EPS1*EPS1 ) + EPS * EPS1
    #
    # and            N = [shape functions matrix]
    #                         1        2        ...        nNodes
    #                  1    [[N_11 N_12
    #                  2     [N_21
    #                  :     [
    #                 nIP    [N_ ...                    ]]


    npoints = shape.npoints
    IP      = get_ip_coords(shape, nips)
    ndim    = shape.ndim # not related with the analysis ndim

    #filling N matrix with shape functions of all ips
    N = Array{Float64}(undef, nips, npoints)
    for i=1:nips
        N[i,:] = shape.func(vec(IP[i,:]))
    end

    #calculate extrapolator matrix
    if nips==npoints
        return inv(N)
    elseif nips>=npoints
        return pinv(N)
    elseif nips==1
        return pinv(N) # Correction procedure is not applicable for nips==1
    end

    # Case for nips<npoints
    #I = eye(nips)

    # εip matrix: Local ip coordinates of integration points
    εip = [ IP[:,1:ndim] ones(nips) ]

    # ε matrix: Local coordinates of nodal points
    ε = [ shape.nat_coords ones(npoints) ] # increase a column of ones

    E = pinv(N)*(I - εip*pinv(εip)) + ε*pinv(εip)

    return E
end
