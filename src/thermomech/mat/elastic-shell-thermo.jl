# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellThermo

mutable struct ElasticShellThermoState<:IpState
    env::ModelEnv
    σ::Array{Float64,1} # stress
    ε::Array{Float64,1} # strain
    QQ::Array{Float64,1} # heat flux
    D::Array{Float64,1}
    ut::Float64
    function ElasticShellThermoState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.QQ = zeros(env.ndim)
        this.D = zeros(env.ndim)
        this.ut = 0.0
        return this
    end
end


mutable struct ElasticShellThermo<:Material
    E ::Float64 # Young's Modulus kN/m2
    nu::Float64 # Poisson coefficient
    k ::Float64 # thermal conductivity  w/m/k
    ρ ::Float64 # density Ton/m3
    cv::Float64 # Specific heat J/Ton/k
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C
    thickness::Float64 # thickness


    function ElasticShellThermo(prms::Dict{Symbol,Float64})
        return  ElasticShellThermo(;prms...)
    end

    function ElasticShellThermo(;E=NaN, nu=NaN, k=NaN, rho=NaN, cv=NaN, alpha=1.0, thickness =NaN)
        E>0.0       || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        isnan(k)     && error("Missing value for k")
        cv<=0.0       && error("Invalid value for cv: $cv")
        0<=alpha<=1  || error("Invalid value for alpha: $alpha")
        thickness>0.0       || error("Invalid value for thick: $thickness")

        this = new(E, nu, k, rho, cv, alpha, thickness)

        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticShellThermo) = TMShell

# Type of corresponding state structure
ip_state_type(mat::ElasticShellThermo) = ElasticShellThermoState


#function set_state(state::ElasticShellLinSeepState; sig=zeros(0), eps=zeros(0))
#    sq2 = √2.0
#    mdl = [1, 1, 1, sq2, sq2, sq2]
#    if length(sig)==6
#        state.σ .= sig.*mdl
#    else
#        if length(sig)!=0; error("ElasticShellLinSeep: Wrong size for stress array: $sig") end
#    end
#    if length(eps)==6
#        state.ε .= eps.*mdl
#    else
#        if length(eps)!=0; error("ElasticShellLinSeep: Wrong size for strain array: $eps") end
#    end
#end

function calcD(mat::ElasticShellThermo, state::ElasticShellThermoState)
    E = mat.E
    ν = mat.nu
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

#=
function stress_update(mat::ElasticShellThermo, state::ElasticShellThermoState, dε::Array{Float64,1})
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end
=#

function stress_update(mat::ElasticShellThermo, state::ElasticShellThermoState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    De = calcD(mat, state)
    Δσ = De*Δε
    state.ε  += Δε
    state.σ  += Δσ

    K = calcK(mat, state)
    state.QQ = -K*G
    state.D  += state.QQ*Δt
    state.ut += Δut
    return Δσ, state.QQ
end

function ip_state_vals(mat::ElasticShellThermo, state::ElasticShellThermoState)
    D = stress_strain_dict(state.σ, state.ε, state.env.modeltype)
    #=
    D[:qx] = state.QQ[1] # VERIFICAR NECESSIDADE
    D[:qy] = state.QQ[2] # VERIFICAR NECESSIDADE
        if state.env.ndim==3 # VERIFICAR NECESSIDADE
            D[:qz] = state.QQ[3] # VERIFICAR NECESSIDADE
        end # VERIFICAR NECESSIDADE
    =#
    return D
    #return stress_strain_dict(state.σ, state.ε, state.env.modeltype)
end
