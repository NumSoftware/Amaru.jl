# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElasticSeep

mutable struct LinearElasticSeepState<:IpState
    env::ModelEnv
    σ::Vec6 # stress
    ε::Vec6 # strain
    V::Array{Float64,1} # fluid velocity
    D::Array{Float64,1} # distance traveled by the fluid
    uw::Float64         # pore pressure
    function LinearElasticSeepState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        this.V = zeros(env.ndim)
        this.D = zeros(env.ndim)
        this.uw = 0.0
        return this
    end
end


mutable struct LinearElasticSeep<:Material
    E ::Float64 # Young's modulus
    nu::Float64 # Poisson ratio
    k ::Float64 # specific permeability
    α ::Float64 # Biot's coefficient
    S ::Float64 # Storativity coefficient

    function LinearElasticSeep(;E=NaN, nu=NaN, k=NaN, kappa=NaN, alpha=1.0, S=0.0, n=NaN, Ks=NaN, Kw=NaN, eta=NaN)

        if isnan(k)
            k = (kappa*gammaw)/eta # specific permeability = (intrinsic permeability * fluid specific weight)/viscosity
        end

        if isnan(S)
            S = (alpha - n)/Ks + n/Kw # S = (alpha - porosity)/(bulk module of the solid) + (porosity)/(bulk module of the fluid)
        end

        @check E>0.0
        @check 0<=nu<0.5
        @check k>0
        @check 0<alpha<=1.0
        @check S>=0.0

        this = new(E, nu, k, alpha, S)
        return this
    end
end


# Type of corresponding state structure
ip_state_type(::HydromechSolid, ::LinearElasticSeep) = LinearElasticSeepState


function calcD(mat::LinearElasticSeep, state::LinearElasticSeepState)
    return calcDe(mat.E, mat.nu, state.env.ana.stressmodel) # function calcDe defined at elastic-solid.jl
end


function calcK(mat::LinearElasticSeep, state::LinearElasticSeepState) # Hydraulic conductivity matrix
    if state.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state(mat::LinearElasticSeep, state::LinearElasticSeepState, Δε::Array{Float64,1}, Δuw::Float64, G::Array{Float64,1}, Δt::Float64)
    De = calcD(mat, state)
    Δσ = De*Δε
    state.ε  += Δε
    state.σ  += Δσ
    K = calcK(mat, state)
    state.V   = -K*G
    state.D  += state.V*Δt
    state.uw += Δuw
    return Δσ, state.V
end


function ip_state_vals(mat::LinearElasticSeep, state::LinearElasticSeepState)
    D = stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)

    D[:vx] = state.V[1]
    D[:vy] = state.V[2]
    if state.env.ndim==3
        D[:vz] = state.V[3]
    end

    return D
end
