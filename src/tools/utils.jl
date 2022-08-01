
function wrap(str::String; level=1)
    width = displaysize(stdout)[2]-2
    str = replace(str, r" +" => " ")
    parts = String[]

    for s in split(str, "\n")
        if (level-1)*2+length(str)>width
            len = width - 2*(level-1) 
            rx  = Regex(".{1,$len}( |\$)")
            ss  = SubstitutionString("  "^(level-1)*"\\0\\n")
            s   = replace(s, rx => ss)[1:end-1]
        else 
            s = "  "^(level-1)*s
        end
        push!(parts, s)
    end
    return join(parts, "\n")
end


function hint(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", msg, "\e[K\n", color=:light_black)
end

function message(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", msg, "\e[K\n", color=:white)
end

function info(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", "$msg\e[K\n", color=:cyan)
end

function notify(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", msg, "\e[K\n", color=:yellow)
end

function alert(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\n", "$msg\e[K\n", color=:light_red)
end

function warn(xs...)
    msg = join(string.(xs))
    msg = wrap(msg, level=2)
    printstyled("\rWarning: \e[K\n", color=:light_yellow)
    printstyled(msg, "\e[K\n", color=:yellow)
end

function headline(xs...; color=:cyan)
    msg = join(string.(xs))
    printstyled("\r", msg, "\e[K\n", color=color, bold=true)
end


export sound_alert
function sound_alert()
    try
        cmd = `beep`
        if Sys.islinux()
            file = joinpath(@__DIR__, "alert.oga")
            cmd = `paplay $file`
        end
        if Sys.iswindows()
            cmd = `cmd /c 'echo \^G'`
        end

        out = Pipe()
        err = Pipe()
        run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
        close(out.in)
        close(err.in)
    catch err
    end
end


struct StopException{T} <: Exception
    S::T
end

function Base.showerror(io::IO, ex::StopException, bt; backtrace=true)
    Base.with_output_color(get(io, :color, false) ? :green : :nothing, io) do io
        showerror(io, ex.S)
    end
end

export stop
stop(text="Stop point.") = throw(StopException(text))

function align(arr::Array{<:AbstractString,1}, pat::AbstractString)
    isempty(arr) && return String[]
    idxs   = [ findnext(pat, str, 1).start for str in arr ]
    maxidx = maximum(idxs)
    return [ replace(str, pat => " "^(maxidx-idxs[i])*pat, count=1)  for (i,str) in enumerate(arr) ]
end