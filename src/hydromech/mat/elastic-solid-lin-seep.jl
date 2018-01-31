# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolidLinSeep

mutable struct ElasticSolidLinSeepIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Array{Float64,1}
    ε::Array{Float64,1}
    uw::Float64
    function ElasticSolidLinSeepIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.uw = 0.0
        return this
    end
end

@show_function ElasticSolidLinSeepIpState

mutable struct ElasticSolidLinSeep<:Material
    E ::Float64
    nu::Float64
    k ::Float64
    gw::Float64

    function ElasticSolidLinSeep(prms::Dict{Symbol,Float64})
        return  ElasticSolidLinSeep(;prms...)
    end

    function ElasticSolidLinSeep(;E=1.0, nu=0.0, k=NaN, gw=NaN)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        isnan(k)     && error("Missing value for k")
        isnan(gw)    && error("Missing value for gw")
        !(gw>0)      && error("Invalid value for gw: $gw")
        this    = new(E, nu, k, gw)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticSolidLinSeep) = HMSolid

# Create a new instance of Ip data
new_ip_state(mat::ElasticSolidLinSeep, shared_data::SharedAnalysisData) = ElasticSolidLinSeepIpState(shared_data)

function set_state(ipd::ElasticSolidLinSeepIpState; sig=zeros(0), eps=zeros(0))
    sq2 = √2.0
    mdl = [1, 1, 1, sq2, sq2, sq2]
    if length(sig)==6
        ipd.σ[:] = sig.*mdl
    else
        if length(sig)!=0; error("ElasticSolidLinSeep: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*mdl
    else
        if length(eps)!=0; error("ElasticSolidLinSeep: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState)
    return calcDe(mat.E, mat.nu, ipd.shared_data.model_type) # function calcDe defined at elastic-solid.jl
end

function calcK(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState) # Hydraulic conductivity matrix
    if ipd.shared_data.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function stress_update(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState, Δε::Array{Float64,1}, Δuw::Float64, G::Array{Float64,1})
    De = calcDe(mat.E, mat.nu, ipd.shared_data.model_type)
    Δσ = De*Δε
    ipd.ε += Δε
    ipd.σ += Δσ
    ipd.uw += Δuw
    K = calcK(mat, ipd)
    V = -K*G
    return Δσ, V
end

function ip_state_vals(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState)
    ndim = ipd.shared_data.ndim
    σ  = ipd.σ
    ε  = ipd.ε
    sr2  = 2.0^0.5

    if ndim==2;
        if ipd.shared_data.model_type == :plane_stress
            return Dict(
                :sxx => σ[1],
                :syy => σ[2],
                :szz => σ[3],
                :sxy => σ[4]/sr2,
                :syz => σ[5]/sr2,
                :sxz => σ[6]/sr2,
                :exx => ε[1],
                :eyy => ε[2],
                :ezz => ε[3],
                :exy => ε[4]/sr2,
                :eyz => ε[5]/sr2,
                :exz => ε[6]/sr2,
                :s_m => sum(σ[1:3])/3.0,
            )
        else
            return Dict(
                :sxx => σ[1],
                :syy => σ[2],
                :szz => σ[3],
                :sxy => σ[4]/sr2,
                :exx => ε[1],
                :eyy => ε[2],
                :ezz => ε[3],
                :exy => ε[4]/sr2,
                :s_m => sum(σ[1:3])/3.0
            )
        end
    else
        return Dict(
            :sxx => σ[1],
            :syy => σ[2],
            :szz => σ[3],
            :sxy => σ[4]/sr2,
            :syz => σ[5]/sr2,
            :sxz => σ[6]/sr2,
            :exx => ε[1],
            :eyy => ε[2],
            :ezz => ε[3],
            :exy => ε[4]/sr2,
            :eyz => ε[5]/sr2,
            :exz => ε[6]/sr2,
            :s_m => sum(σ[1:3])/3.0
        )
    end
end
