# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint

mutable struct JointIpState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function JointIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

mutable struct ElasticJoint<:Material
    E::Float64 # Young modulus from bulk material
    ν::Float64 # Poisson ration from bulk material
    kn::Float64 # Normal stiffness (used only if E and ν are NaN)
    ks::Float64 # Shear stiffness (used only if E and ν are NaN)
    ζ::Float64 # elastic displacement scale factor (formerly α)

    function ElasticJoint(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJoint(;E=NaN, nu=NaN, kn=NaN, ks=NaN, zeta=1.0)
        # kn and ks are used only if E and ν are NaN

        if isnan(kn*ks)
            @check E>0.0    
            @check 0<=nu<0.5
        else
            @check kn>0.0    
            @check ks>0.0    
        end
        @check zeta>0

        this = new(E, nu, kn, ks, zeta)
        return this
    end

end

# Returns the element type that works with this material model
@static if @isdefined MechJoint
    matching_elem_type(::ElasticJoint) = MechJoint
end

# Type of corresponding state structure
ip_state_type(mat::ElasticJoint) = JointIpState


function mountD(mat::ElasticJoint, ipd::JointIpState)
    ndim = ipd.env.ndim
    if isnan(mat.kn*mat.ks)
        G  = mat.E/(1.0+mat.ν)/2.0
        kn = mat.E*mat.ζ/ipd.h
        ks =     G*mat.ζ/ipd.h
    else
        kn = mat.kn
        ks = mat.ks
    end

    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end

function stress_update(mat::ElasticJoint, ipd::JointIpState, Δu)
    ndim = ipd.env.ndim
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.w[1:ndim] += Δu
    ipd.σ[1:ndim] += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::ElasticJoint, ipd::JointIpState)
    ndim = ipd.env.ndim
    if ndim == 3
       return Dict(
          :jw1  => ipd.w[1],
          :jw2  => ipd.w[2],
          :jw3  => ipd.w[3],
          :js1  => ipd.σ[1],
          :js2  => ipd.σ[2],
          :js3  => ipd.σ[3],
          )
    else
        return Dict(
          :jw1  => ipd.w[1],
          :jw2  => ipd.w[2],
          :js1  => ipd.σ[1],
          :js2  => ipd.σ[2],
          )
    end
end

function output_keys(mat::ElasticJoint)
    return Symbol[:jw1, :jw1]
end