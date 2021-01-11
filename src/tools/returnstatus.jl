mutable struct ReturnStatus
    success::Bool
    message::String
    function ReturnStatus(success::Bool=true, message::String="")
        return new(success, message)
    end
end

failed(rs::ReturnStatus)    = !rs.success
succeeded(rs::ReturnStatus) =  rs.success

failure(msg::String="")      = ReturnStatus(false, msg)
Base.success(msg::String="") = ReturnStatus(true, msg)
