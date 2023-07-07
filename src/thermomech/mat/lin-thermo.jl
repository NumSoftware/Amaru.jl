# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinThermo

mutable struct LinThermoState<:IpState
    env::ModelEnv
    ut::Float64
    QQ::Array{Float64,1}
    D::Array{Float64,1} #
    function LinThermoState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.ut = 0.0
        this.QQ  = zeros(env.ndim)
        this.D  = zeros(env.ndim)
        return this
    end
end


mutable struct LinThermo<:MatParams
    k ::Float64 # thermal conductivity with/m/K
    #α ::Float64 #  coefficient of thermal expansion 1/K or 1/°C

    function LinThermo(;k=NaN)
        @check k>=0.0
        this = new(k)
        return this
    end
end



# Type of corresponding state structure
ip_state_type(::ThermoSolidElem, ::LinThermo) = LinThermoState


function calcK(matparams::LinThermo, state::LinThermoState) # Thermal conductivity matrix
    if state.env.ndim==2
        return matparams.k*eye(2)
    else
        return matparams.k*eye(3)
    end
end


function update_state!(matparams::LinThermo, state::LinThermoState, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(matparams, state)
    state.QQ   = -K*G
    state.D  += state.QQ*Δt
    state.ut += Δut
    return state.QQ
end


function ip_state_vals(matparams::LinThermo, state::LinThermoState)
    D = OrderedDict{Symbol, Float64}()
    #D[:qx] = state.QQ[1]
    #D[:qy] = state.QQ[2]
    #if state.env.ndim==3
        #D[:qz] = state.QQ[3]
    #end
    return D
end
