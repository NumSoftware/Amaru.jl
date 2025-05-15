# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function encode_string_to_uint64(s::AbstractString)
    codes = codeunits(s)
    @assert length(codes) ≤ 8 "Can only encode up to 8 ASCII characters"
    x = UInt64(0)
    for (i, b) in enumerate(codeunits(s))
        x |= UInt64(b) << (8 * (i - 1))
    end
    return x
end


function decode_uint64_to_string(x::UInt64)
    bytes = UInt8[]
    for i in 0:7
        b = UInt8((x >> (8 * i)) & 0xff)
        b == 0x00 && break
        push!(bytes, b)
    end
    return String(bytes)
end


function safe_string_cut(s::AbstractString, n_bytes::Int)
    if n_bytes > ncodeunits(s)
        return (s, "")
    end

    p = prevind(s, 8)
    if nextind(s,p)==n_bytes
        return (s[1:n_bytes], s[n_bytes+1:end])
    else
        return (s[1:p], s[p+1:end])
    end
end