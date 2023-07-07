export Color, mix, rgb, rgba, Rgb, Rgba

struct Color
    r::Float64
    g::Float64
    b::Float64
    a::Float64

    function Color(c::Color)
        return new(c.r, c.g, c.b)
    end

    function Color(r::Float64, g::Float64, b::Float64, a::Float64=1.0)
        return new(
            clamp(r, 0.0, 1.0), 
            clamp(g, 0.0, 1.0), 
            clamp(b, 0.0, 1.0), 
            clamp(a, 0.0, 1.0)
        )
    end

    function Color(r::Int, g::Int, b::Int, a::Float64=1.0)
        return new(
            clamp(r, 0, 255)/255, 
            clamp(g, 0, 255)/255, 
            clamp(b, 0, 255)/255, 
            clamp(a, 0.0, 1.0)
        )
    end

    function Color(t::Tuple)
        a::Float64 = length(t)==4 ? t[4] : 1.0
        return Color(t[1], t[2], t[3], a)
    end

    # Color((r, g, b)::NTuple{4,Float64}) = Color(r, g, b)
    # Color((r, g, b, a)::NTuple{4,Float64}) = Color(r, g, b, a)
    # Color((r, g, b)::Tuple{Int,Int,Int}) = Color(r, g, b, a)
    # Color((r, g, b, a)::Tuple{Int,Int,Int,Float64}) = Color(r, g, b, a)
end


function mix(c1::Color, a::Number, c2::Color=Color(1.0,1.0,1.0))
    return Color(
        c1.r*a + c2.r*(1-a),
        c1.g*a + c2.g*(1-a),
        c1.b*a + c2.b*(1-a),
        c1.a*a + c2.a*(1-a)
    )
end


function mix(c1, a::Number, c2=Color(1.0,1.0,1.0))
    return mix(Color(c1), a, Color(c2))
end

rgb(c::Color)  = (c.r, c.g, c.b)
rgba(c::Color) = (c.r, c.g, c.b, c.a)
Rgb(c::Color)  = round.(Int, (c.r*255, c.g*255, c.b*255))
Rgba(c::Color) = round.(Int, (c.r*255, c.g*255, c.b*255, c.a))


rgb(c1, a=1.0, c2=Color(1.0, 1.0, 1.0))  = rgb(mix(c1, a, c2))
rgba(c1, a=1.0, c2=Color(1.0, 1.0, 1.0)) = rgba(mix(c1, a, c2))
Rgb(c1, a=1.0, c2=Color(1.0, 1.0, 1.0))  = Rgb(mix(c1, a, c2))
Rgba(c1, a=1.0, c2=Color(1.0, 1.0, 1.0)) = Rgba(mix(c1, a, c2))

# rgb(blue, 0.5
# base.:*(a, c::rbg) = rgb(a*c[1], a*c[2], a*c[3])
# base.:*(c::rbg, a) = rgb(a*c[1], a*c[2], a*c[3])
# base.:+(c::rbg, a) = rgb(a*c[1], a*c[2], a*c[3])