# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ConstConductivityShell

mutable struct ConstConductivityShellState<:IpState
    env::ModelEnv
    ut::Float64
    QQ::Array{Float64,1}
    D::Array{Float64,1} #
    function ConstConductivityShellState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.ut = 0.0
        this.QQ  = zeros(env.ndim)
        this.D  = zeros(env.ndim)
        return this
    end
end

mutable struct ConstConductivityShell<:Material
    k ::Float64 # Thermal conductivity w/m/k

    function ConstConductivityShell(; params...)
        names = (k="Conductivity")
        required = (:k, )
        @checkmissing params required names

        params  = values(params)
        k       = params.k

        @check k>=0.0
        return new(k)
    end
end



# Type of corresponding state structure
ip_state_type(::Type{ConstConductivityShell}) = ConstConductivityShellState


function calcK(mat::ConstConductivityShell, state::ConstConductivityShellState) # Thermal conductivity matrix
    if state.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state!(mat::ConstConductivityShell, state::ConstConductivityShellState, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(mat, state)
    state.QQ   = -K*G
    state.D  += state.QQ*Δt
    state.ut += Δut
    return state.QQ
end


function ip_state_vals(mat::ConstConductivityShell, state::ConstConductivityShellState)
    D = OrderedDict{Symbol, Float64}()
    #D[:qx] = state.QQ[1]
    #D[:qy] = state.QQ[2]
    #if state.env.ndim==3
        #D[:qz] = state.QQ[3]
    #end
    return D
end
