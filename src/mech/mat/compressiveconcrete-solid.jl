# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CompressiveConcrete

mutable struct CompressiveConcrete<:MatParams
    E0   ::Float64  # initial Young modulus
    ν    ::Float64
    fc   ::Float64
    εpeak::Float64
    β    ::Float64
    ρ    ::Float64

    function CompressiveConcrete(; E=NaN, nu=0.2, fc=NaN, epspeak=NaN, rho=0.0, beta=NaN)
        @check E>0.0  
        @check nu>=0.0 
        @check fc<0.0  
        @check epspeak<0.0
        @check rho>=0.0
        @check epspeak/(fc/E)>1.0
        if isnan(beta) 
            beta = clamp(inv(1-fc/(epspeak*E)), 1.0, 10.0) # limit the value of β
        end
        @check beta>=1.0

        this = new(E, nu, fc, epspeak, beta, rho)
        return this
    end
end


mutable struct CompressiveConcreteState<:IpState
    env        ::ModelEnv
    σ          ::Array{Float64,1}  # current stress
    ε          ::Array{Float64,1}  # current strain
    ε̅c         ::Float64
    ε̅min       ::Float64

    function CompressiveConcreteState(env::ModelEnv=ModelEnv())
        this      = new(env)
        this.σ    = zeros(6)
        this.ε    = zeros(6)
        this.ε̅c   = 0.0
        this.ε̅min = 0.0

        return this
    end
end




# Type of corresponding state structure
ip_state_type(matparams::CompressiveConcrete) = CompressiveConcreteState


function uniaxial_σ(matparams::CompressiveConcrete, state::CompressiveConcreteState, εi::Float64)
    # σp = eigvals(state.σ)
    # σ1c, σ2c, σ3c = neg.(σp)

    # compression: Popovics 1973; Carreira and Chu 1985
    # αc = 0.3
    # γ = 1.0 - αc/matparams.fc*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c)) # suggested values for coef: 0.2
    γ = 1.0
    fc = matparams.fc*γ

    β  = 1/(1-fc/(matparams.εpeak*matparams.E0))
    β  = max(min(β,10),2) # limit the value of β
    εr = εi/matparams.εpeak
    return fc*(β*εr)/(β - 1.0 + εr^β)

end


function uniaxial_E(matparams::CompressiveConcrete, state::CompressiveConcreteState, εi::Float64)
    # σp = eigvals(state.σ)
    # σ1c, σ2c, σ3c = neg.(σp)

    # αc = 0.3
    # γ = 1.0 - αc/matparams.fc*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c)) # suggested values for coef: 0.2
    γ = 1.0
    fc = matparams.fc*γ

    εpeak = matparams.εpeak*γ
    εr  = εi/εpeak

    β = 1/(1-fc/(εpeak*matparams.E0))
    β = max(min(β,10),2) # limit the value of β
    Ec = β*(fc/εpeak)/(β-1+εr^β) - β^2*(fc/εpeak)*εr^β/(β-1+εr^β)^2
    return Ec
end


function calcD(matparams::CompressiveConcrete, state::CompressiveConcreteState)

    if state.ε̅c > state.ε̅min
        E = matparams.E0
    else
        E = uniaxial_E(matparams, state, state.ε̅c)
        Emin = matparams.E0*1e-4
        abs(E)<Emin && (E=Emin)
    end

    D  = calcDe(E, matparams.ν, state.env.anaprops.stressmodel)
    return D
end


function update_state(matparams::CompressiveConcrete, state::CompressiveConcreteState, Δε::Array{Float64,1})
    # special function
    neg(x) = (-abs(x)+x)/2.0

    state.ε .+= Δε
    εp = eigvals(state.ε)
    state.ε̅c = -norm(neg.(εp), 2)

    if state.ε̅c > state.ε̅min
        E = matparams.E0
    else
        state.ε̅min = state.ε̅c
        E = uniaxial_E(matparams, state, state.ε̅c)
    end

    D  = calcDe(E, matparams.ν, state.env.anaprops.stressmodel)
    Δσ = D*Δε
    state.σ .+= Δσ

    return Δσ, success()
end


function ip_state_vals(matparams::CompressiveConcrete, state::CompressiveConcreteState)
    dict = stress_strain_dict(state.σ, state.ε, state.env.anaprops.stressmodel)
    dict[:Ec] = uniaxial_E(matparams, state, state.ε̅c)
    return dict
end
