# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export XmlDocument, XmlElement

abstract type XmlNode end

mutable struct XmlElement<:XmlNode
    name::String
    attributes::OrderedDict{String,String}
    children::Vector{XmlNode}
    content::Union{AbstractString, Vector{UInt8}}

    # function XmlElement(name::AbstractString, attributes::Union{AbstractDict, Tuple}, children, content::AbstractString="")
        # return new(name, OrderedDict(attributes), children, content)
    # end

    function XmlElement(name::AbstractString; attributes::Union{AbstractDict,Tuple}=Dict(), children=XmlElement[], content::AbstractString="")
        return new(name, OrderedDict(attributes), children, content)
    end
end

addchild!(node::XmlElement, child::XmlNode) = push!(node.children, child)

mutable struct XmlComment<:XmlNode
    content::String
    function XmlComment(content::AbstractString)
        return new(content)
    end
end


haschildren(node::XmlElement) = length(node.children)>0


mutable struct XmlDocument
    attributes::OrderedDict{String,String}
    children::Vector{XmlNode}
    root::XmlElement
    function XmlDocument(attributes::Union{AbstractDict,Tuple})
        return new(OrderedDict(attributes), [])
    end
    function XmlDocument(attributes::Union{AbstractDict,Tuple}, root::XmlElement)
        return new(OrderedDict(attributes), [root], root)
    end
end


function addchild!(node::XmlDocument, child::XmlNode) 
    push!(node.children, child)
    if child isa XmlElement
        node.root = child
    end
end

# Get a list of all nodes with a given attribute
function Base.getindex(doc::XmlDocument, p::Pair{String,String})
    return getindex(doc.root, p)
end


# Get a node from a nested sequence of names
# If among the children there are more than one node with the 
# same name, the last one is considered.
function (node::XmlElement)(args::String...)
    n = node
    for s in args
        found = false
        for child in n.children[end:-1:1]
            if child isa XmlElement && child.name==s
                found = true
                n = child
                break
            end
        end
        found || return nothing
    end
    return n
end

# Get a list of all nodes with a given name
function Base.getindex(node::XmlElement, s::String)
    nodes = XmlElement[]
    for child in node.children
        if child isa XmlElement && child.name==s
            push!(nodes, child)
        end
    end
    return nodes
end

# Get a child according to index
function Base.getindex(node::XmlElement, i::Int)
    i<=length(node.children) && return node.children[i]
    return nothing
end

# Get a child according to s
function getchild(node::XmlElement, s::String)
    for child in node.children
        if child isa XmlElement && child.name==s
            return child
        end
    end
    return nothing
end



# Get a list of all nodes with a given attribute and value (att=>val)
function Base.getindex(node::XmlElement, p::Pair{String,String})
    nodes = XmlElement[]

    att = p.first
    val = p.second

    get(node.attributes, att, nothing) == val && push!(nodes, node)

    for child in node.children
        child isa XmlElement || continue
        append!(nodes, Base.getindex(child, p))
    end

    return nodes
end


# Read a xml node. Return nothing if there is no node.
function readnode(text, pos)
    # check if there is no node
    m = match(r"(\S{2})", text, pos)
    if m===nothing
        return nothing, pos
    elseif m.captures[1]=="</"  # empty content
        return nothing, pos
    elseif m.captures[1][1]!='<' # non empty content
        return nothing, pos
    end
    
    # get node name and attributes range
    rng = findnext(r" *<.+?>", text, pos)
    pos = rng.stop
    str = SubString(text, rng)

    # check if comment
    if str[2]=='!'
        # @show str
        # @show pos
        # rng = findnext("-->", text, pos)
        # pos = rng.stop
        return XmlComment(strip(str[5:end-4])), pos
    end

    name = match(r"<(\w+)", str).captures[1]
    attributes = OrderedDict{String,String}()
    
    # get attributes
    apos = 1
    while true
        m = match(r" (\w+)=\"(.*?)\"", str, apos)
        m===nothing && break
        _,apos = findnext(r" \w+=\".*?\"", str, apos)
        att = m.captures[1]
        val = m.captures[2]
        attributes[att] = val
    end

    # check if self-closing tag
    if str[end-1:end]=="/>"
        return XmlElement(name, attributes=attributes), pos
    end

    # get children or content
    children = XmlNode[]
    content = ""
    while true
        node,pos = readnode(text, pos+1)
        if node!==nothing
            push!(children, node)
        else
            rng = findnext("</$name>", text, pos) 
            content = strip(text[pos:rng.start-1])
            pos = rng.stop
            break
        end
    end

    return XmlElement(name, attributes=attributes, children=children, content=content), pos

end


# read a xml file into a XmlDocument
function XmlDocument(filename::String)
    text = read(filename, String)
    rng = findfirst(r"<\?xml.*?>", text)
    pos = rng.stop

    # get attributes range
    rng = findfirst(r"<\?xml.*?>", text)
    pos = rng.stop
    str = SubString(text, rng)

    # get attributes
    attributes = OrderedDict{String,String}()
    apos = 1
    while true
        m = match(r" (\w+)=\"(.*?)\"", str, apos)
        m===nothing && break
        _,apos = findnext(r" \w+=\".*?\"", str, apos)
        att = m.captures[1]
        val = m.captures[2]
        attributes[att] = val
    end

    # get children or content
    children = XmlNode[]
    while true
        node,pos = readnode(text, pos+1)
        if node!==nothing
            push!(children, node)
        else
            # rng = findnext("</$name>", text, pos) 
            # content = strip(text[pos:rng.start-1])
            # pos = rng.stop
            break
        end
    end

    idx = findfirst(e-> e isa XmlElement, children)
    root = children[idx]

    return XmlDocument(attributes, root)
end


function writenode(f::IOStream, node::XmlComment, level::Int)
    tab = "   "
    println(f, tab^level, "<!-- ", node.content, " -->")
end


function writenode(f::IOStream, node::XmlElement, level::Int)
    tab = "   "
    # Print header
    print(f, tab^level, "<", node.name)
    for (att,val) in node.attributes
        print(f, " $att=\"$val\"")
    end

    # Self-closing tag
    if length(node.children)==0 && node.content==""
        println(f, "/>")
        return
    end

    # Ending header
    println(f, ">")

    # Print content or children
    if length(node.children)==0
        # print content
        if node.name=="DataArray"
            println(f, replace(node.content, r"^ *"m=>tab^(level+1)))
        elseif node.attributes["encoding"]=="raw"
            write(f, "_")
            write(f, node.content)
            println(f)
        else
            println(f, tab^(level+1), node.content)
        end
    else
        for child in node.children
            writenode(f, child, level+1)
        end
    end

    # Print closing tag
    println(f, tab^level, "</", node.name, ">")
end


function save(doc::XmlDocument, filename::String)
    # Open filename
    f = open(filename, "w")

    # Print prolog
    print(f, "<?xml")
    for (att,val) in doc.attributes
        print(f, " $att=\"$val\"")
    end
    println(f, "?>")

    # Print nodes
    for child in doc.children
        writenode(f, child, 0)
    end

    close(f)
end


export to_xml_node
const Xsingletype = Union{Number,Bool,AbstractString,Char,Symbol}


function to_xml_node(arr::Array{<:Any}, name::String="Array", attributes::AbstractDict=Dict{String,String}())
    if eltype(arr) <: Xsingletype
        attributes["format"] = "ascii"
        attributes["type"] = string(eltype(arr))
        attributes["components"] = string(size(arr,2))
        if length(arr)>0
            str = repr("text/plain", arr)
            content = str[findfirst('\n',str)+2 : end]
        else
            content = ""
        end

        return XmlElement(name, attributes=attributes, content=content)
    end

    attributes["type"] = string(eltype(arr))
    children = XmlElement[]
    for item in arr
        push!( children, to_xml_node(item) )
    end

    return XmlElement(name, attributes=attributes, children=children)
end


function to_xml_node(dict::AbstractDict, name::String="Dict", attributes::AbstractDict=Dict{String,String}())
    # vty = valtype(dict)
    s = string(typeof(dict))
    s = split(s,".")[end]
    attributes["type"] = s

    children = [ 
                to_xml_node(collect(keys(dict)), "keys"),
                to_xml_node(collect(values(dict)), "values")
               ]
    #children = XmlElement[]
    #for (k,v) in dict
        #if typeof(v) <: Xsingletype
            #push!(children, XmlElement(string(k), string(v)))
        #else
            #push!(children, XmlElement(string(k), Dict(), [to_xml_node(v)]))
        #end
    #end

    return XmlElement(name, attributes=attributes, children=children)

end


function to_xml_node(obj::Any, name::String=""; exclude::Array{Symbol,1}=Symbol[])
    name=="" && (name=string(typeof(obj)))
    attributes=Dict{String,String}()
    children = XmlElement[]
    ty = typeof(obj)
    fields = fieldnames(ty)
    types  = fieldtypes(ty)

    for (fld,ty) in zip(fields,types)
        fld in exclude && continue
        fld_str = string(fld)
        fld_str[1] == '_' && continue

        value = getfield(obj, fld)

        if ty<:Xsingletype
            attributes[string(fld)] = string(value)
        else
            push!(children, to_xml_node(value, fld_str))
        end
    end

    return XmlElement(name, attributes=attributes, children=children)
    
end


function XmlElement(obj::Any, name::String=""; exclude::Array{Symbol,1}=Symbol[])
    return to_xml_node(obj, name, exclude=exclude)
end
