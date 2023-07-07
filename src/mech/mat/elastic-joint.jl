# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint

mutable struct JointState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function JointState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

mutable struct ElasticJoint<:MatParams
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

@static if @isdefined MechJoint
    
end

# Type of corresponding state structure
ip_state_type(matparams::ElasticJoint) = JointState


function mountD(matparams::ElasticJoint, state::JointState)
    ndim = state.env.ndim
    if isnan(matparams.kn*matparams.ks)
        G  = matparams.E/(1.0+matparams.ν)/2.0
        kn = matparams.E*matparams.ζ/state.h
        ks =     G*matparams.ζ/state.h
    else
        kn = matparams.kn
        ks = matparams.ks
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


function update_state(matparams::ElasticJoint, state::JointState, Δu)
    ndim = state.env.ndim
    D  = mountD(matparams, state)
    Δσ = D*Δu

    state.w[1:ndim] += Δu
    state.σ[1:ndim] += Δσ
    return Δσ, success()
end


function ip_state_vals(matparams::ElasticJoint, state::JointState)
    ndim = state.env.ndim
    if ndim == 3
       return Dict(
          :jw1  => state.w[1],
          :jw2  => state.w[2],
          :jw3  => state.w[3],
          :js1  => state.σ[1],
          :js2  => state.σ[2],
          :js3  => state.σ[3],
          )
    else
        return Dict(
          :jw1  => state.w[1],
          :jw2  => state.w[2],
          :js1  => state.σ[1],
          :js2  => state.σ[2],
          )
    end
end


function output_keys(matparams::ElasticJoint)
    return Symbol[:jw1, :jw1]
end