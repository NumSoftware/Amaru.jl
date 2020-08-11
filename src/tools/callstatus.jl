mutable struct CallStatus
    success::Bool
    message::String
    function CallStatus(success::Bool=true, message::String="")
        return new(success, message)
    end
end
