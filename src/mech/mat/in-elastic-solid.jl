# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export InElasticSolid

mutable struct InElasticSolidIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Array{Float64,1}
    ε::Array{Float64,1}
    function InElasticSolidIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end

@show_function InElasticSolidIpState

mutable struct InElasticSolid<:Material
    E ::Float64
    nu::Float64
    ρ::Float64

    function InElasticSolid(prms::Dict{Symbol,Float64})
        return  InElasticSolid(;prms...)
    end

    function InElasticSolid(;E=1.0, nu=0.0, rho=NaN)
        if E<=0.0      ; error("Invalid value for E: $E") end
        if !(0<=nu<0.5); error("Invalid value for nu: $nu") end
        this    = new(E, nu, rho)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::InElasticSolid) = InMechSolid

# Create a new instance of Ip data
new_ip_state(mat::InElasticSolid, shared_data::SharedAnalysisData) = InElasticSolidIpState(shared_data)

function set_state(ipd::InElasticSolidIpState; sig=zeros(0), eps=zeros(0))
    sq2 = √2.0
    mdl = [1, 1, 1, sq2, sq2, sq2]
    if length(sig)==6
        ipd.σ[:] = sig.*mdl
    else
        if length(sig)!=0; error("InElasticSolid: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*mdl
    else
        if length(eps)!=0; error("InElasticSolid: Wrong size for strain array: $eps") end
    end
end

function calcDe(E::Number, nu::Number, atype::Symbol)
    if atype==:plane_stress
        c = E/(1.0-nu^2)
        return [
            c    c*nu 0. 0. 0. 0.
            c*nu c    0. 0. 0. 0.
            0.   0.   0. 0. 0. 0.
            0.   0.   0. c*(1.-nu) 0. 0.
            0.   0.   0. 0. 0. 0.
            0.   0.   0. 0. 0. 0. ]
    else
        c = E/((1.0+nu)*(1.0-2.0*nu))
        return [
            c*(1-nu) c*nu     c*nu     0.         0.         0.
            c*nu     c*(1-nu) c*nu     0.         0.         0.
            c*nu     c*nu     c*(1-nu) 0.         0.         0.
            0.       0.       0.       c*(1-2*nu) 0.         0.
            0.       0.       0.       0.         c*(1-2*nu) 0.
            0.       0.       0.       0.         0.         c*(1-2*nu) ]
    end
end

function calcD(mat::InElasticSolid, ipd::InElasticSolidIpState)
    return calcDe(mat.E, mat.nu, ipd.shared_data.model_type)
    #return mat.De
end

function stress_update(mat::InElasticSolid, ipd::InElasticSolidIpState, dε::Array{Float64,1})
    De = calcDe(mat.E, mat.nu, ipd.shared_data.model_type)
    dσ = De*dε
    ipd.ε += dε
    ipd.σ += dσ
    return dσ
end

function ip_state_vals(mat::InElasticSolid, ipd::InElasticSolidIpState)
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
                :s_m => sum(σ[1:3])/3.0 
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
