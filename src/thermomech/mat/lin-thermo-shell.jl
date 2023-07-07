# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinThermoShell

mutable struct LinThermoShellState<:IpState
    env::ModelEnv
    ut::Float64
    QQ::Array{Float64,1}
    D::Array{Float64,1} #
    function LinThermoShellState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.ut = 0.0
        this.QQ  = zeros(env.ndim)
        this.D  = zeros(env.ndim)
        return this
    end
end

mutable struct LinThermoShell<:MatParams
    k ::Float64 # Thermal conductivity w/m/k
    ρ ::Float64 # density Ton/m3
    cv::Float64 # Specific heat J/Ton/k
    thickness::Float64 # thickness
    #E ::Float64 # Young's modulus kPa
    #nu::Float64 # Poisson ratio
    #α ::Float64 #  coefficient of thermal expansion 1/K or 1/°C

    function LinThermoShell(prms::Dict{Symbol,Float64})
        return  LinThermoShell(;prms...)
    end

    function LinThermoShell(;k=NaN, rho=NaN, cv=NaN, thickness =NaN)
        k  >= 0.0 || error("Invalid value for k")
        rho>= 0.0 || error("Invalid value for rho")
        cv >= 0.0 || error("Invalid value for cv")
        thickness>0.0       || error("Invalid value for thick: $thickness")
        this = new(k,rho,cv,thickness)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::LinThermoShell) = ThermoShellElem

# Type of corresponding state structure
ip_state_type(matparams::LinThermoShell) = LinThermoShellState


function calcK(matparams::LinThermoShell, state::LinThermoShellState) # Thermal conductivity matrix
    if state.env.ndim==2
        return matparams.k*eye(2)
    else
        return matparams.k*eye(3)
    end
end


function update_state!(matparams::LinThermoShell, state::LinThermoShellState, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(matparams, state)
    state.QQ   = -K*G
    state.D  += state.QQ*Δt
    state.ut += Δut
    return state.QQ
end


function ip_state_vals(matparams::LinThermoShell, state::LinThermoShellState)
    D = OrderedDict{Symbol, Float64}()
    #D[:qx] = state.QQ[1]
    #D[:qy] = state.QQ[2]
    #if state.env.ndim==3
        #D[:qz] = state.QQ[3]
    #end
    return D
end
