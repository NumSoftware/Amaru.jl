# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

@enum(CellFamily,
    VERTEXCELL    = 0,
    LINECELL      = 1,
    BULKCELL      = 2,
    EMBEDDEDCELL  = 3,
    JOINTCELL     = 4,
    LINEJOINTCELL = 5,
    TIPJOINTCELL  = 6,
)
const SOLIDCELL = BULKCELL

# Export
for s in instances(CellFamily)
    @eval export $(Symbol(s))
end
export SOLIDCELL

mutable struct CellShape
    name       ::String
    family     ::CellFamily
    ndim       ::Int
    npoints    ::Int
    basic_shape::CellShape
    vtk_type   ::VTKCellType
    facet_idxs ::Array
    edge_idxs  ::Array
    facet_shape::Union{CellShape, Tuple}
    nat_coords ::Array
    quadrature ::Dict{Int, Array}
    func       ::Function
    deriv      ::Function
    deriv2     ::Function
    function CellShape()
        return new()
    end
end

include("shapes/lines.jl")
include("shapes/solids2d.jl")
include("shapes/solids3d.jl")
include("shapes/joints.jl")

# Shape for unknown polyvertex
function MakePOLYVERTEX()
    shape             = CellShape()
    shape.name        = "POLYVERTEX"
    shape.family      = VERTEXCELL
    shape.ndim        = 0
    shape.npoints     = 0
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.facet_shape = ()
    shape.quadrature  = Dict( 0 => [] )
    return shape
end
const POLYVERTEX = MakePOLYVERTEX()
export POLYVERTEX


function get_ip_coords(shape::CellShape, nips=0)
    !haskey(shape.quadrature, nips) && error("Cannot set $nips integ. points for shape $shape. Available numbers are $(collect(keys(shape.quadrature))).")
    return shape.quadrature[nips]
end


# Available VTK shapes
const VTK2SHAPE = Dict{VTKCellType,CellShape}(
    VTK_POLY_VERTEX             => POLYVERTEX,
    VTK_LINE                    => LIN2,
    VTK_TRIANGLE                => TRI3,
    VTK_QUAD                    => QUAD4,
    VTK_PYRAMID                 => PYR5,
    VTK_TETRA                   => TET4,
    VTK_HEXAHEDRON              => HEX8,
    VTK_WEDGE                   => WED6,
    VTK_QUADRATIC_EDGE          => LIN3,
    VTK_QUADRATIC_TRIANGLE      => TRI6,
    VTK_QUADRATIC_QUAD          => QUAD8,
    VTK_QUADRATIC_TETRA         => TET10,
    VTK_QUADRATIC_HEXAHEDRON    => HEX20,
    VTK_TRIQUADRATIC_HEXAHEDRON => HEX27,
    VTK_QUADRATIC_WEDGE         => WED15,
    VTK_BIQUADRATIC_QUAD        => QUAD9
)

# All isoparametric cell shapes
const ALL_ISO_SHAPES = [
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
    PYR13
    TET4
    TET10
    WED6
    WED15
    HEX8
    HEX20
    HEX27
]

# dictionary (ndim,nfpoints,nlayers) => joint_shape
_joint_ndim_nfpoints_nlayers_dict = Dict( 
    (2,2,2)=>:JLIN2 , (2,3,2)=>:JLIN3 , (2,4,2)=>:JLIN4 , (3,3,2)=>:JTRI3 , (3,4,2)=>:JQUAD4 , (3,6,2)=>:JTRI6 , (3,8,2)=>:JQUAD8,
    (2,2,3)=>:J3LIN2, (2,3,2)=>:J3LIN3, (2,4,2)=>:J3LIN4, (3,3,2)=>:J3TRI3, (3,4,2)=>:J3QUAD4, (3,6,2)=>:J3TRI6, (3,8,2)=>:J3QUAD8,
)

function get_shape_from_vtk(vtk_type::VTKCellType, npoints::Int64, ndim::Int64, nlayers::Int64=0)::CellShape
    # vtk_type: VTK cell code
    # npoints : total number of cell points
    # ndim    : analysis dimension
    # nlayers : number of layers for joint elements

    if vtk_type==VTK_POLYGON
        if npoints==12
            return QUAD12    
        end
    end

    vtk_type==VTK_POLY_VERTEX || return VTK2SHAPE[vtk_type]

    # Check if it is a joint cell with layers
    if nlayers in (2,3)
        # dictionary (ndim,nfpoints,nlayers) => joint_shape
        shapedict = _joint_ndim_nfpoints_nlayers_dict
        nfpoints = div(npoints,nlayers)
        if haskey(shapedict, (ndim, nfpoints, nlayers))
            return shapedict[(ndim, nfpoints, nlayers)]
        end
    end

    # Check for other cells
    if     npoints==2   return JLINK2
    elseif npoints==3   return JLINK3
    elseif npoints==9   return TRI9
    elseif npoints==10  return TRI10
    end

    error("get_shape_from_vtk: Unknown shape for vtk_type $vtk_type and npoints $npoints with ndim $ndim")
end


function get_shape_from_ndim_npoints(npoints::Int64, ndim::Int64)::CellShape
    # npoints : total number of cell points
    # ndim    : analysis dimension

    if ndim==3
        npoints==2  &&  return LIN2
        npoints==4  &&  return TET4
        npoints==10 &&  return TET10
        npoints==8  &&  return HEX8
        npoints==20 &&  return HEX20
        npoints==27 &&  return HEX27
    else
        npoints==2  && return LIN2
        npoints==3  && return TRI3
        npoints==6  && return TRI6
        npoints==4  && return QUAD4
        npoints==8  && return QUAD8
        npoints==12 && return QUAD12
    end

    npoints==1 && return POLYVERTEX

    error("get_shape_from_ndim_npoints: Cannot get the cell shape from npoints=$npoints and ndim=$ndim")
end


function bdistance(shape::CellShape, R::Array{Float64,1})
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


function inverse_map(shape::CellShape, coords::Array{Float64,2}, X0::AbstractArray{Float64,1}, tol=1.0e-7)
    maxits = 20
    ndim   = shape.ndim
    R      = zeros(ndim)
    C      = coords[:,1:ndim]
    X      = X0[1:ndim]

    local ΔX::Array{Float64,1}
    
    for k in 1:maxits
        # calculate Jacobian
        D = shape.deriv(R)
        J = C'*D

        # calculate trial of real coordinates
        N  = shape.func(R)
        Xt = C'*N # interpolating

        # calculate the error
        ΔX = Xt - X
        # ΔR = pinv(J)'*ΔX
        ΔR = inv(J)*ΔX # don't use pinv

        # updating local coords R
        R -= ΔR

        if norm(ΔX) < tol; break end
    end

    # TODO: Improve accuracy of inverse_map function in elements with non regular shape
    #k==maxits && println("Warning: max iterations (maxits=$maxits) reached in inverse mapping. norm(ΔX)=$(norm(ΔX))")

    if ndim==2
        R = vcat( R, 0.0 )
    end
    return R
end


function is_inside(shape::CellShape, C::Array{Float64,2}, X::Array{Float64,1}, tol = 1.e-7)
    if shape.family!=SOLIDCELL return false end

    # Testing with bounding box
    Cmin = vec(minimum(C, dims=1))
    Cmax = vec(maximum(C, dims=1))
    maxl = maximum(Cmax-Cmin)
    ttol = 0.1*maxl # 10% is important for curved elements

    if any(X .< Cmin .- ttol) || any(X .> Cmax .+ ttol)
        return false
    end

    # Testing with inverse mapping
    R = inverse_map(shape, C, X, tol)
    if bdistance(shape, R) > -tol
        return true
    else
        return false
    end
end



function extrapolator(shape::CellShape, nips::Int)
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
    for i in 1:nips
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
