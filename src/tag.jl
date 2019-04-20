# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Tag functions
for T in (Node, Element, Ip, Edge, Face)
    @eval begin
        tag!(object::$T, tag::String) = (object.tag = tag)
        tag!(objects::Array{$T,1}, tag::String) = for object in objects; object.tag = tag end
    end
end

@doc """
    tag!(object, tag)

Modify the field `object.tag` with the value `tag`.

    tag!(itr, tag)
Tag each object in `itr` with the value `tag`.
""" tag!

