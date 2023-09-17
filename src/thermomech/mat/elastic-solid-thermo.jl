# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElasticThermo

mutable struct LinearElasticThermoState<:IpState
    env::ModelEnv
    σ::Vec6 # stress
    ε::Vec6 # strain
    QQ::Array{Float64,1} # heat flux
    D::Array{Float64,1}
    ut::Float64
    function LinearElasticThermoState(env::ModelEnv)
        this    = new(env)
        this.σ  = zeros(Vec6)
        this.ε  = zeros(Vec6)
        this.QQ = zeros(env.ndim)
        this.D  = zeros(env.ndim)
        this.ut = 0.0
        return this
    end
end


mutable struct LinearElasticThermo<:Material
    E ::Float64 # Young's Modulus kN/m2
    ν::Float64 # Poisson coefficient
    k ::Float64 # thermal conductivity  w/m/k
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C

    function LinearElasticThermo(; params...)
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
compat_state_type(::Type{LinearElasticThermo}, ::Type{TMSolid}, env::ModelEnv) = LinearElasticThermoState


function calcD(mat::LinearElasticThermo, state::LinearElasticThermoState)
    return calcDe(mat.E, mat.ν, state.env.ana.stressmodel) # function calcDe defined at elastic-solid.jl
end


function calcK(mat::LinearElasticThermo, state::LinearElasticThermoState) # Thermal conductivity matrix
    if state.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state!(mat::LinearElasticThermo, state::LinearElasticThermoState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
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


function ip_state_vals(mat::LinearElasticThermo, state::LinearElasticThermoState)
    D = stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)

    #D[:qx] = state.QQ[1] # VERIFICAR NECESSIDADE
    #D[:qy] = state.QQ[2] # VERIFICAR NECESSIDADE
    #if state.env.ndim==3 # VERIFICAR NECESSIDADE
        #D[:qz] = state.QQ[3] # VERIFICAR NECESSIDADE
    #end # VERIFICAR NECESSIDADE

    return D
end
