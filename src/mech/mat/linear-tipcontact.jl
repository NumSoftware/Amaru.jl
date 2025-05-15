# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TipContact

mutable struct TipJointState<:IpState
    ctx::Context
    f ::Float64
    w ::Float64
    function TipJointState(ctx::Context)
        this = new(ctx)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end

TipJoint_params = [
    FunInfo( :TipContact, "A model for a bar tip contact."),
    KwArgInfo( :k, "Elastic stiffness", 1.0, cond=:(k>=0) ),
    KwArgInfo( :fixed, "Flag to control if the tip is fixed", false, type=Bool),
]
@doc docstring(TipJoint_params) TipContact

mutable struct TipContact<:Material
    k::Float64
    fixed::Bool

    function TipContact(;args...)
        args = checkargs(args, TipJoint_params)
        this = new(args.k, args.fixed)
        return this
    end
end


# Element types that work with this material
# compat_elem_types(::Type{TipContact}) = (MechTipJoint,)

# Type of corresponding state structure
compat_state_type(::Type{TipContact}, ::Type{MechTipJoint}, evn::Context) = TipJointState


function calcD(mat::TipContact, state::TipJointState)
    if state.w>0.0 || mat.fixed
        return mat.k
    else
        return 0.0
    end
end


function update_state!(mat::TipContact, state::TipJointState, Δw)
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


function ip_state_vals(mat::TipContact, state::TipJointState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
