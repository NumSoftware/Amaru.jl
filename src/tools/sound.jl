
export sound_alert

function sound_alert()
    try
        file = joinpath(@__DIR__, "alert.oga")
        #@show file
        run(`paplay $file`);
        #read(`paplay $file`);
    catch
        nothing
    end
end
