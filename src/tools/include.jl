# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

include("returnstatus.jl")
include("iteration.jl")

include("color.jl")
include("show.jl")

include("constants.jl")
include("math.jl")
include("signal.jl")

include("linalg.jl")
# include("vec3.jl")
# include("vec6.jl")
# include("mat6x6.jl")
include("quaternion.jl")
include("tensors.jl")

include("expr.jl")
include("arguments.jl")
include("threads.jl")
include("table.jl")
include("book.jl")
include("utils.jl")
include("stopwatch.jl")
include("xml.jl")

include("array-pool.jl")

include("tex.jl")
include("error.jl")

Base.show(io::IO, obj::XmlElement) = _show(io, obj, 3, "")
Base.show(io::IO, obj::XmlDocument)  = _show(io, obj, 3, "")
