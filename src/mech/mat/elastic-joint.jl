# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint

mutable struct JointState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function JointState(env::ModelEnv)
        this = new(env)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

mutable struct ElasticJoint<:Material
    E::Float64 # Young modulus from bulk material
    ν::Float64 # Poisson ratio from bulk material
    kn::Float64 # Normal stiffness (used only if E and ν are NaN)
    ks::Float64 # Shear stiffness (used only if E and ν are NaN)
    ζ::Float64  # elastic displacement scale factor (formerly α)

    function ElasticJoint(; params...)
        names = (E="Young modulus", nu="Poisson ratio", kn="Normal stiffness per area", ks="shear stiffness per area", 
        zeta="elastic displacement scale factor")
        
        required = (:zeta,)
        @checkmissing params required names

        default = (E=NaN, nu=NaN, kn=NaN, ks=NaN)
        params  = merge(default, params)

        E       = params.E
        nu      = params.nu
        kn      = params.kn
        ks      = params.ks
        zeta      = params.zeta

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

# mat_arguments(::Type{ElasticJoint})  = [ 
            # @arg E=NaN E>0 "UM"
            # @arg nu=NaN nu>0 "UM"
            # @arg kn=NaN kn>0 "UM"
            # @arg ks=NaN ks>0 "UM"
            # @argcond ks*ks>0
            # @argopt (E, nu) (kn, ks) 
#             Arg(:E, "Young modulus", :(E>0), true),
#             Arg(:nu, "Poisson ratio", :(nu>0), true),
#             Arg(:kn, "Normal stiffness per area", :(kn>0), true),
#             ArgOpt((:E,:nu), (:kn,:ks))
#             ArgOpt(:kn,:E)
#         ]


# Type of corresponding state structure
compat_state_type(::Type{ElasticJoint}, ::Type{MechJoint}, env::ModelEnv) = JointState


function calcD(mat::ElasticJoint, state::JointState)
    ndim = state.env.ndim
    if isnan(mat.kn*mat.ks)
        G  = mat.E/(1.0+mat.ν)/2.0
        kn = mat.E*mat.ζ/state.h
        ks =     G*mat.ζ/state.h
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


function update_state!(mat::ElasticJoint, state::JointState, Δu)
    ndim = state.env.ndim
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.w[1:ndim] += Δu
    state.σ[1:ndim] += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::ElasticJoint, state::JointState)
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


function output_keys(mat::ElasticJoint)
    return Symbol[:jw1, :jw1]
end