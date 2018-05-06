# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinSeep

mutable struct LinSeepIpState<:IpState
    shared_data::SharedAnalysisData
    uw::Float64
    function LinSeepIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.uw = 0.0
        return this
    end
end

@show_function LinSeepIpState

mutable struct LinSeep<:Material
    k ::Float64
    gw::Float64

    function LinSeep(prms::Dict{Symbol,Float64})
        return  LinSeep(;prms...)
    end

    function LinSeep(;k=NaN, gw=NaN)
        isnan(k)     && error("Missing value for k")
        isnan(gw)    && error("Missing value for gw")
        !(gw>0)      && error("Invalid value for gw: $gw")
        this    = new(k, gw)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::LinSeep) = SeepSolid

# Create a new instance of Ip data
new_ip_state(mat::LinSeep, shared_data::SharedAnalysisData) = LinSeepIpState(shared_data)

function set_state(ipd::LinSeepIpState; uw=0.0)
end

function calcK(mat::LinSeep, ipd::LinSeepIpState) # Hydraulic conductivity matrix
    if ipd.shared_data.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function ip_state_vals(mat::LinSeep, ipd::LinSeepIpState)
    D = Dict{Symbol, Float64}()
    #D[:uw] = ipd.uw  # TODO

    return D
end
