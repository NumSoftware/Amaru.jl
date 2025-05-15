# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearTipContact

mutable struct LinearTipContactState<:IpState
    ctx::Context
    f ::Float64
    w ::Float64
    function LinearTipContactState(ctx::Context)
        this = new(ctx)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end

LinearTipContact_params = [
    FunInfo( :LinearTipContact, "A model for a bar tip contact."),
    KwArgInfo( :k, "Elastic stiffness", 1.0, cond=:(k>=0) ),
    KwArgInfo( :fixed, "Flag to control if the tip is fixed", false, type=Bool),
]
@doc docstring(LinearTipContact_params) LinearTipContact

mutable struct LinearTipContact<:Material
    k::Float64
    fixed::Bool

    function LinearTipContact(;args...)
        args = checkargs(args, LinearTipContact_params)
        this = new(args.k, args.fixed)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearTipContact}, ::Type{MechSlipTip}, evn::Context) = LinearTipContactState


function calcD(mat::LinearTipContact, state::LinearTipContactState)
    if state.w>0.0 || mat.fixed
        return mat.k
    else
        return 0.0
    end
end


function update_state!(mat::LinearTipContact, state::LinearTipContactState, Δw)
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


function ip_state_vals(mat::LinearTipContact, state::LinearTipContactState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
