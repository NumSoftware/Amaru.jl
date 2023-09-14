# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearAcousticFluid

mutable struct LinearAcousticFluidState<:IpState
    env::ModelEnv
    up ::Float64          # pressure
    # V  ::Array{Float64,1} # fluid velocity?
    
    function LinearAcousticFluidState(env::ModelEnv)
        this = new(env)
        # this.V  = zeros(env.ndim)
        this.up = 0.0
        return this
    end
end


mutable struct LinearAcousticFluid<:Material
    ρ::Float64 # densidade
    c::Float64 # sound speed
    # μ::Float64 # viscocity

    function LinearAcousticFluid(;rho=NaN, c=NaN)
        @check rho>0
        @check c>0

        this = new(rho, c)
        return this
    end
end



# Type of corresponding state structure
compat_state_type(::Type{LinearAcousticFluid}, ::Type{AcousticFluid}, env::ModelEnv) = LinearAcousticFluidState

# Element types that work with this material
# compat_elem_types(::Type{LinearAcousticFluid}) = (AcousticFluid,)


function update_state!(mat::LinearAcousticFluid, state::LinearAcousticFluidState, Δup::Float64, G::Array{Float64,1}, Δt::Float64)
    state.up += Δup
    return nothing
end


function ip_state_vals(mat::LinearAcousticFluid, state::LinearAcousticFluidState)
    D = OrderedDict{Symbol, Float64}()
    
    # D[:up] = state.up

    return D
end
