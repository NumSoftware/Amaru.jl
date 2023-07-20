# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellThermo

mutable struct ElasticShellThermoState<:IpState
    env::ModelEnv
    σ::Vec6 # stress
    ε::Vec6 # strain
    QQ::Array{Float64,1} # heat flux
    D::Array{Float64,1}
    ut::Float64
    function ElasticShellThermoState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        this.QQ = zeros(env.ndim)
        this.D = zeros(env.ndim)
        this.ut = 0.0
        return this
    end
end


mutable struct ElasticShellThermo<:Material
    E ::Float64 # Young's Modulus kN/m2
    ν::Float64 # Poisson coefficient
    k ::Float64 # thermal conductivity  w/m/k
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C

    function ElasticShellThermo(; params...)
        names = (E="Young modulus", nu="Poisson ratio", k="Conductivity", alpha="Thermal expansion coefficient")
        required = (:E, :k, :nu, :alpha)
        @checkmissing params required names

        params = (; params...)
        E      = params.E
        nu     = params.nu
        k      = params.k
        alpha  = params.alpha

        @check E>=0.0
        @check 0<=nu<0.5
        @check k>0
        @check 0<=alpha<=1
        return new(E, nu, k, alpha)
    end
end


# Type of corresponding state structure
ip_state_type(::Type{ElasticShellThermo}) = ElasticShellThermoState

# Element types that work with this material
matching_elem_types(::Type{ElasticShellThermo}) = (TMShell,)


function calcD(mat::ElasticShellThermo, state::ElasticShellThermoState, stressmodel="shell")
    E = mat.E
    ν = mat.ν
    c = E/(1.0-ν^2)
    g = E/(1+ν)
    return [
        c    c*ν   0.0  0.0    0.0    0.0
        c*ν  c     0.0  0.0    0.0    0.0
        0.0  0.0   0.0  0.0    0.0    0.0
        0.0  0.0   0.0  5/6*g  0.0    0.0
        0.0  0.0   0.0  0.0    5/6*g  0.0
        0.0  0.0   0.0  0.0    0.0    g ]
    # ezz = -ν/E*(sxx+syy)
end


function calcK(mat::ElasticShellThermo, state::ElasticShellThermoState) # Thermal conductivity matrix
    if state.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state!(mat::ElasticShellThermo, state::ElasticShellThermoState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64, stressmodel="shell")
    De = calcD(mat, state)
    Δσ = De*Δε
    state.ε  += Δε
    state.σ  += Δσ

    K = calcK(mat, state)
    state.QQ = -K*G
    state.D  += state.QQ*Δt
    state.ut += Δut
    return Δσ, state.QQ, success()
end


function ip_state_vals(mat::ElasticShellThermo, state::ElasticShellThermoState, stressmodel="shell")
    D = stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
    #=
    D[:qx] = state.QQ[1] # VERIFICAR NECESSIDADE
    D[:qy] = state.QQ[2] # VERIFICAR NECESSIDADE
        if state.env.ndim==3 # VERIFICAR NECESSIDADE
            D[:qz] = state.QQ[3] # VERIFICAR NECESSIDADE
        end # VERIFICAR NECESSIDADE
    =#
    return D
    #return stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
end
