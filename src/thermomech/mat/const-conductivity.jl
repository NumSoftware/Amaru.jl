# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ConstConductivity


mutable struct ConstConductivityState<:IpState
    ctx::Context
    ut::Float64
    Q::Array{Float64,1}
    function ConstConductivityState(ctx::Context)
        this = new(ctx)
        this.ut = 0.0
        this.Q  = zeros(ctx.ndim)
        return this
    end
end


ConstConductivity_params = [
    FunInfo(:ConstConductivity, "Constant thermal conductivity material model"),
    KwArgInfo(:k, "Conductivity", cond=:(k>=0)),
    KwArgInfo(:cv, "Specific heat", 0.0, cond=:(cv>=0))
]
@doc docstring(ConstConductivity_params) ConstConductivity(; kwargs...)


mutable struct ConstConductivity<:Material
    k ::Float64 # thermal conductivity with/m/K
    cv::Float64

    function ConstConductivity(; kwargs...)
        args = checkargs(kwargs, ConstConductivity_params)
        return new(args.k, args.cv)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{ConstConductivity}, ::Type{ThermoSolid}, ctx::Context) = ConstConductivityState
compat_state_type(::Type{ConstConductivity}, ::Type{ThermoShell}, ctx::Context) = ConstConductivityState


function calc_cv(mat::ConstConductivity, ut::Float64) # Specific heat
    return mat.cv
end

function calcK(mat::ConstConductivity, state::ConstConductivityState) # Thermal conductivity matrix
    if state.ctx.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state!(mat::ConstConductivity, state::ConstConductivityState, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(mat, state)
    q = -K*G
    state.Q  += q*Δt
    state.ut += Δut
    return q
end


function ip_state_vals(mat::ConstConductivity, state::ConstConductivityState)
    D = OrderedDict{Symbol, Float64}()
    #D[:qx] = state.q[1]
    #D[:qy] = state.q[2]
    #if state.ctx.ndim==3
        #D[:qz] = state.q[3]
    #end
    return D
end
