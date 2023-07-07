# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearAcousticFluid

mutable struct LinearAcousticFluidState<:IpState
    env::ModelEnv
    up ::Float64          # pressure
    # V  ::Array{Float64,1} # fluid velocity?
    
    function LinearAcousticFluidState(env::ModelEnv=ModelEnv())
        this = new(env)
        # this.V  = zeros(env.ndim)
        this.up = 0.0
        return this
    end
end


mutable struct LinearAcousticFluid<:MatParams
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
ip_state_type(::LinearAcousticFluid) = LinearAcousticFluidState


function update_state!(matparams::LinearAcousticFluid, state::LinearAcousticFluidState, Δup::Float64, G::Array{Float64,1}, Δt::Float64)
    return nothing
end


function ip_state_vals(matparams::LinearAcousticFluid, state::LinearAcousticFluidState)
    D = OrderedDict{Symbol, Float64}()
    D[:up] = state.up
    # D[:vx] = state.V[1]
    # D[:vy] = state.V[2]
    # if state.env.ndim==3
    #     D[:vz] = state.V[3]
    # end

    return D
end
