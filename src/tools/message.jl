
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
    printstyled("\r", msg, "\033[K\n", color=:light_black)
end

function message(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", msg, "\033[K\n", color=:white)
end

function info(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", "$msg\033[K\n", color=:cyan)
end

function notify(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", msg, "\033[K\n", color=:yellow)
end

function alert(xs...; level=2)
    msg = join(string.(xs))
    msg = wrap(msg, level=level)
    printstyled("\r", "$msg\033[K\n", color=:light_red)
end

function warn(xs...)
    msg = join(string.(xs))
    msg = wrap(msg, level=2)
    printstyled("\rWarning: \033[K\n", color=:light_yellow)
    printstyled(msg, "\033[K\n", color=:yellow)
    #printstyled("\rWarning: ", color=:yellow, bold=true)
    #printstyled(msg, "\033[K\n\n", color=:yellow)
end

function headline(xs...; color=:cyan)
    msg = join(string.(xs))
    printstyled("\r", msg, "\033[K\n", color=color, bold=true)
end
