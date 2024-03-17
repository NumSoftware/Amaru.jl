# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function getfullpath(dir, filename)
    filename=="" && return ""
    isabspath(filename) && return filename

    fullpath = joinpath(dir, filename)
    dir = dirname(fullpath)
    dir!="" && !isdir(dir) && error("getfullpath: Directory $(dir) does not exist.")

    return fullpath
end


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


function align(arr::Array{<:AbstractString,1}, pat::AbstractString)
    isempty(arr) && return String[]
    idxs   = [ findnext(pat, str, 1).start for str in arr ]
    maxidx = maximum(idxs)
    return [ replace(str, pat => " "^(maxidx-idxs[i])*pat, count=1)  for (i,str) in enumerate(arr) ]
end

export latex
function latex(M::Array; digits=3)
    l = IOBuffer()

    m = size(M, 1)
    n = size(M, 2)
    

    # widths calculation
    etype = eltype(M)
    width = etype<:Integer ? 6 : 12-digits

    # printing header
    level = 1
    indent = "    "
    println(l, indent^level, raw"\begin{pmatrix}" )
    level = 2
    # printing body
    for i in 1:m
        print(l, indent^level)
        for j in 1:n
            item = M[i,j]
            if etype<:AbstractFloat
                if isnan(item)
                    item = "NaN"
                else
                    item = sprintf("%$width.$(digits)f", item)
                end
                print(l, lpad(string(item), width))
            elseif etype<:Integer
                item = @sprintf("%6d", item)
                print(l, lpad(item, width))
            else
                print(l, rpad(item, width))
            end
            j<n && print(l, " & ")
        end
        println(l, raw" \\\\")
    end
    # printing ending
    level = 1
    println(l, indent^level, raw"\end{pmatrix}")

    return String(take!(l))
end

