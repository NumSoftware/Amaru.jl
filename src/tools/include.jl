include("returnstatus.jl")
include("iteration.jl")

include("color.jl")
include("show.jl")

include("constants.jl")
include("math.jl")

include("linalg.jl")
include("vec3.jl")
include("quaternion.jl")
include("tensors.jl")

include("expr.jl")
include("threads.jl")
include("table.jl")
include("book.jl")
include("utils.jl")
include("stopwatch.jl")
include("xml.jl")

include("tex.jl")

Base.show(io::IO, obj::Xnode) = _show(io, obj, 3, "")
Base.show(io::IO, obj::Xdoc)  = _show(io, obj, 3, "")
