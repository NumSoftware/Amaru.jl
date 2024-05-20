mutable struct Stage
    id          ::Int
    "boundary conditions"
    bcs         ::AbstractArray
    nincs       ::Int
    nouts       ::Int
    tspan       ::Float64
    toactivate  ::Array{<:Element,1}
    todeactivate::Array{<:Element,1}
    status      ::Symbol # :pending, :solving, :done, :error

    function Stage(bcs         ::AbstractArray;
                   nincs       ::Int     = 1,
                   nouts       ::Int     = 0,
                   tspan       ::Number  = 0.0,
                   toactivate  ::Array{<:Element,1}=Element[],
                   todeactivate::Array{<:Element,1}=Element[],
    )
        @check nincs>0
        @check nouts>=0
        return new(-1, bcs, nincs, nouts, tspan, toactivate, todeactivate, :pending)
    end
end
