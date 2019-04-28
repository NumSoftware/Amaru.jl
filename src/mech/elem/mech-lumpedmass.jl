# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechLumpedMass<:Mechanical
    id    ::Int
    shape ::ShapeType
    cell  ::Cell
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function MechLumpedMass()
        return new() 
    end
end

matching_shape_family(::Type{MechLumpedMass}) = VERTEX_SHAPE


function elem_mass(elem::MechLumpedMass)
    ndim = elem.env.ndim
    M = mat.m*Matrix{Float64}(I, ndim, ndim)

    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return M, map, map
end
            
                        

function elem_vals(elem::MechLumpedMass)
    vals = OrderedDict()
    return vals
end
