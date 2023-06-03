# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Tag functions
for T in (Node, Cell, Block, Element, Ip)
    
    @eval begin
        """
            $(TYPEDSIGNATURES)

        Tag a $($T) `object` by setting the `tag` string.
        """
        tag!(object::$T, tag::String) = (object.tag = tag; nothing)

        """
            $(SIGNATURES)

        Tag all $($T) `objects` in an array using the `tag` string.
        """
        tag!(objects::Array{<:$T,1}, tag::String) = for object in objects; object.tag = tag end
        
    end
    
end

tag!(::Nothing, tag::String) = error("There are no objects to apply tag \"$tag\". The tag! function received 'nothing'.")