mutable struct StopWatch
    lapse::Float64
    lasttime::Float64
    stopped::Bool
    function StopWatch()
        new(0.0, time(), false)
    end
end


function stop!(sw::StopWatch)
    if sw.stopped == false
        sw.lapse += time() - sw.lasttime
        sw.stopped = true
    end
    return nothing
end


function resume!(sw::StopWatch)
    if sw.stopped
        sw.lasttime = time()
        sw.stopped  = false
    end
    return nothing
end


function getlapse(sw::StopWatch)
    sw.stopped && return sw.lapse
    t = time()
    if t-sw.lasttime>30 # ignore if time interval is greater that 30 seconds
        sw.lasttime = t
    end
    sw.lapse = sw.lapse + t - sw.lasttime
    sw.lasttime = t
    return sw.lapse
end


function formatlapse(lapse::Float64; format::Symbol=:short)
    # format may be :hms, :ms, :s or :short
    if format==:s
        return @sprintf("%5.3fs", lapse)
    elseif format==:ms
        m, s = divrem(lapse, 60)
        return @sprintf("%1dm %5.2fs", m, s)
    else
        h, r = divrem(lapse, 3600)
        m, s = divrem(r, 60)
        #@show format
        if format==:hms
            return @sprintf("%1dh %1dm %5.2fs", h, m, s)
        else
            return @sprintf("%1d:%1d:%3.2f", h, m, s)
        end
    end
end


function see(sw::StopWatch; format::Symbol=:short)
    return formatlapse(getlapse(sw), format=format)
end