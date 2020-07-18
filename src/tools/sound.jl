
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
