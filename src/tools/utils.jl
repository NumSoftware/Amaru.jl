
function wrap(str::String; level=1)
    width = displaysize(stdout)[2]-2
    (level-1)*2+length(str)<=width && return "  "^(level-1)*str

    str = replace(str, r"(\s|\n)+" => " ")
    len = width - 2*(level-1) 
    rx  = Regex(".{1,$len}( |\$)")
    ss  = SubstitutionString("  "^(level-1)*"\\0\\n")
    str = replace(str, rx => ss)
    return str[1:end-1]
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
    #printstyled("\rWarning: ", color=:yellow, bold=true)
    #printstyled(msg, "\e[K\n\n", color=:yellow)
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


export @stop
macro stop()
    return :(error("Stopped here."))
end


function align(arr::Array{<:AbstractString,1}, pat::AbstractString)
    isempty(arr) && return String[]
    idxs   = [ findnext(pat, str, 1).start for str in arr ]
    maxidx = maximum(idxs)
    return [ replace(str, pat => " "^(maxidx-idxs[i])*pat, count=1)  for (i,str) in enumerate(arr) ]
end