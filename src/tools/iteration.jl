
# Flatten a multilevel nested array
function unfold(N)
    typeof(N) <: Union{Tuple,AbstractArray} || return [N]
    U = []
    for x in N
        if typeof(x) <: Union{Tuple,AbstractArray}
            append!(U, unfold(x))
        else
            push!(U, x)
        end
    end
    return U
end
