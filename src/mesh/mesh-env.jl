# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MeshEnv
    ndim::Int                  # Mesh dimension

    function MeshEnv(ndim=3)
        return new(ndim)
    end
end