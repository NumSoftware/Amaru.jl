# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElasticThermo

mutable struct LinearElasticThermoState<:IpState
    env::ModelEnv
    σ::Array{Float64,1} # stress
    ε::Array{Float64,1} # strain
    QQ::Array{Float64,1} # heat flux
    D::Array{Float64,1}
    ut::Float64
    function LinearElasticThermoState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.QQ = zeros(env.ndim)
        this.D = zeros(env.ndim) 
        this.ut = 0.0
        return this
    end
end


mutable struct LinearElasticThermo<:MatParams
    E ::Float64 # Young's Modulus kN/m2
    nu::Float64 # Poisson coefficient
    k ::Float64 # thermal conductivity  w/m/k
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C

    function LinearElasticThermo(prms::Dict{Symbol,Float64})
        return  LinearElasticThermo(;prms...)
    end

    function LinearElasticThermo(;E=NaN, nu=NaN, k=NaN, alpha=1.0)
        @check E>0.0
        @check 0<=nu<0.5
        @check k>0
        @check 0<=alpha<=1

        this = new(E, nu, k, alpha)

        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::LinearElasticThermo) = TMSolidElem

# Type of corresponding state structure
ip_state_type(::TMSolidElem, ::LinearElasticThermo) = LinearElasticThermoState


function calcD(matparams::LinearElasticThermo, state::LinearElasticThermoState)
    return calcDe(matparams.E, matparams.nu, state.env.anaprops.stressmodel) # function calcDe defined at elastic-solid.jl
end


function calcK(matparams::LinearElasticThermo, state::LinearElasticThermoState) # Thermal conductivity matrix
    if state.env.ndim==2
        return matparams.k*eye(2)
    else
        return matparams.k*eye(3)
    end
end


function update_state(matparams::LinearElasticThermo, state::LinearElasticThermoState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    De = calcD(matparams, state)
    Δσ = De*Δε
    state.ε  += Δε
    state.σ  += Δσ
    K = calcK(matparams, state)
    state.QQ = -K*G
    state.D  += state.QQ*Δt
    state.ut += Δut
    return Δσ, state.QQ
end


function ip_state_vals(matparams::LinearElasticThermo, state::LinearElasticThermoState)
    D = stress_strain_dict(state.σ, state.ε, state.env.anaprops.stressmodel)

    #D[:qx] = state.QQ[1] # VERIFICAR NECESSIDADE
    #D[:qy] = state.QQ[2] # VERIFICAR NECESSIDADE
    #if state.env.ndim==3 # VERIFICAR NECESSIDADE
        #D[:qz] = state.QQ[3] # VERIFICAR NECESSIDADE
    #end # VERIFICAR NECESSIDADE

    return D
end
