# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolid


mutable struct ElasticSolid<:Material
    E ::Float64
    nu::Float64
    ρ::Float64

    function ElasticSolid(prms::Dict{Symbol,Float64})
        return  ElasticSolid(;prms...)
    end

    function ElasticSolid(;E=1.0, nu=0.0, rho=0.0)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        rho<0.0      && error("Invalid value for rho: $rho")
        this = new(E, nu, rho)
        return this
    end
end


mutable struct ElasticSolidIpState<:IpState
    env::ModelEnv
    σ::Array{Float64,1}
    ε::Array{Float64,1}
    function ElasticSolidIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end



# Returns the element type that works with this material model
matching_elem_type(::ElasticSolid) = MechSolid

# Type of corresponding state structure
ip_state_type(mat::ElasticSolid) = ElasticSolidIpState


function calcDe(E::Number, ν::Number, modeltype::String)
    if modeltype=="plane-stress"
        c = E/(1.0-ν^2)
        return [
            c    c*ν   0.0  0.0  0.0  0.0
            c*ν  c     0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  c*(1.0-ν) ]
        # ezz = -ν/E*(sxx+syy)
    else
        c = E/((1.0+ν)*(1.0-2.0*ν))
        return [
            c*(1-ν) c*ν     c*ν     0.0         0.0         0.0
            c*ν     c*(1-ν) c*ν     0.0         0.0         0.0
            c*ν     c*ν     c*(1-ν) 0.0         0.0         0.0
            0.0     0.0     0.0     c*(1-2*ν)   0.0         0.0
            0.0     0.0     0.0     0.0         c*(1-2*ν)   0.0
            0.0     0.0     0.0     0.0         0.0         c*(1-2*ν) ]
    end
end

function calcD(mat::ElasticSolid, ipd::ElasticSolidIpState)
    return calcDe(mat.E, mat.nu, ipd.env.modeltype)
end

function stress_update(mat::ElasticSolid, ipd::ElasticSolidIpState, dε::Array{Float64,1})
    De = calcDe(mat.E, mat.nu, ipd.env.modeltype)
    dσ = De*dε
    ipd.ε += dε
    ipd.σ += dσ
    return dσ, CallStatus(true)
end

function ip_state_vals(mat::ElasticSolid, ipd::ElasticSolidIpState)
    return stress_strain_dict(ipd.σ, ipd.ε, ipd.env)
end
