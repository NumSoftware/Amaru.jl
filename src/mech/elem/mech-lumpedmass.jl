# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechLumpedMassElem<:MechElem
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function MechLumpedMass()
        return new()
    end
end

matching_shape_family(::Type{MechLumpedMassElem}) = VERTEXCELL


function elem_stiffness(elem::MechLumpedMassElem)
    ndim = elem.env.ndim
    matparams  = elem.matparams
    K = zeros(ndim, ndim)

    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_mass(elem::MechLumpedMassElem)
    ndim = elem.env.ndim
    matparams  = elem.matparams

    M = matparams.m*Matrix{Float64}(I, ndim, ndim)

    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return M, map, map
end


function update_elem!(elem::MechLumpedMassElem, U::Array{Float64,1}, Δt::Float64)
    return Float64[], Int[], success()
end


function elem_vals(elem::MechLumpedMassElem)
    vals = OrderedDict()
    return vals
end
