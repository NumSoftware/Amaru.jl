# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechSpring<:Mechanical
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

    function MechSpring()
        return new() 
    end
end

matching_shape_family(::Type{MechSpring}) = LINE_SHAPE

function elem_stiffness(elem::MechSpring)
    ndim   = elem.env.ndim

    if ndim==1
        kx = mat.kx
        K = [ kx -kx
             -kx  kx ]
    elseif ndim==2
        kx, ky = mat.kx, mat.ky
        K = [ kx    0  -kx    0
               0   ky    0  -ky
             -kx    0   kx    0
               0  -ky    0   ky ]
    else 
        kx, ky, kz = mat.kx, mat.ky, mat.kz

        K = [ kx    0    0  -kx    0    0
               0   ky    0    0  -ky    0
               0    0   kz    0    0  -kz
             -kx    0    0   kx    0    0
               0  -ky    0    0   ky    0
               0    0  -kz    0    0   kz ]
    end

    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end
            

#function elem_mass(elem::MechSpring) 
                #
    #ndim   = elem.env.ndim
#
    #M = zeros(2*ndim, 2*ndim)
#
    #keys = [:ux, :uy, :uz][1:ndim]
    #map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    #return M, map, map
#end
                        
                        
function elem_update!(elem::MechSpring, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim = elem.env.ndim
    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    K  = elem_stiffness(elem)
    dU = U[map]
    dF = K*dU

    F[map] .+= dF
    return
end

function elem_vals(elem::MechSpring)
    vals = OrderedDict(:fx => 0.0 )
    return vals
end
