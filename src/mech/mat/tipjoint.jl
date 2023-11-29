# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TipJoint

mutable struct TipJointState<:IpState
    env::ModelEnv
    f ::Float64
    w ::Float64
    function TipJointState(env::ModelEnv)
        this = new(env)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end


mutable struct TipJoint<:Material
    k::Float64
    fixed::Bool

    function TipJoint(;args...)
        args = checkargs(args, func_params(TipJoint))
        this = new(args.k, args.fixed)
        return this
    end
end

func_params(::Type{TipJoint}) = [
    ArgInfo( :k, "Elastic stiffness", 1.0, condition=:(k>=0) ),
    ArgInfo( :fixed, "Flag to control if the tip is fixed", false, type=Bool),
]


# Element types that work with this material
# compat_elem_types(::Type{TipJoint}) = (MechTipJoint,)

# Type of corresponding state structure
compat_state_type(::Type{TipJoint}, ::Type{MechTipJoint}, evn::ModelEnv) = TipJointState


function calcD(mat::TipJoint, state::TipJointState)
    if state.w>0.0 || mat.fixed
        return mat.k
    else
        return 0.0
    end
end


function update_state!(mat::TipJoint, state::TipJointState, Δw)
    fini = state.f
    ftr  = fini + mat.k*Δw
    
    if ftr>0.0 || mat.fixed
        f = ftr
    else
        f = 0.0
    end
    Δf = f - fini
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


function ip_state_vals(mat::TipJoint, state::TipJointState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
