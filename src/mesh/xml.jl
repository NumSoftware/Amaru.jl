mutable struct Xnode
    name::String
    attributes::OrderedDict{String,String}
    children::Array{Xnode,1}
    content::String

    function Xnode(name::AbstractString, attributes::AbstractDict=Dict{String,String}())
        return new(name, attributes, Xnode[], "")
    end
    function Xnode(name::AbstractString, attributes::AbstractDict, children::Array, content::AbstractString)
        return new(name, attributes, children, content)
    end
end

haschildren(node::Xnode) = length(node.children)>0

# Get node by name
function Base.getindex(node::Xnode, s::String)
    for i=length(node.children):-1:1
        if node.children[i].name==s
            return node.children[i]
        end
    end
    error("child $s not found")
    return node
end


# Get a list of all nodes with a given attribute
function Base.getindex(node::Xnode, p::Pair{String,String})
    nodes = Xnode[]
    att = p.first
    val = p.second

    get(node.attributes, att, nothing) == val && push!(nodes, node)

    for child in node.children
        append!(nodes, getindex(child, p))
    end

    return nodes
end


# Read a xml node. Return nothing if there is no node.
function readnode(text, pos)
    # check if there is no node
    m = match(r"(\S{2})", text, pos)
    if m.captures[1]=="</"  # empty content
        return nothing, pos
    elseif m.captures[1][1]!='<' # non empty content
        return nothing, pos
    end
    
    # get node name and attributes range
    rng = findnext(r" *<.+?>", text, pos)
    pos = rng.stop
    str = SubString(text, rng)
    name = match(r"<(\w+)", str).captures[1]
    attributes = OrderedDict{String,String}()
    
    # get attributes
    apos = 1
    while true
        m = match(r" (\w+)=\"(.*?)\"", str, apos)
        m==nothing && break
        _,apos = findnext(r" \w+=\".*?\"", str, apos)
        att = m.captures[1]
        val = m.captures[2]
        attributes[att] = val
    end

    # get children or content
    children = Xnode[]
    content = ""
    while true
        node,pos = readnode(text, pos+1)
        if node!=nothing
            push!(children, node)
        else
            rng = findnext("</$name>", text, pos) 
            content = strip(text[pos:rng.start-1])
            pos = rng.stop
            break
        end
    end

    return Xnode(name, attributes, children, content), pos

end


mutable struct Xdoc
    attributes::OrderedDict{String,String}
    root::Xnode
end

# Get a list of all nodes with a given attribute
function Base.getindex(doc::Xdoc, p::Pair{String,String})
    return getindex(doc.root, p)
end

export Xdoc

# read a xml file into a Xdoc
function Xdoc(filename::String)
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
        m==nothing && break
        _,apos = findnext(r" \w+=\".*?\"", str, apos)
        att = m.captures[1]
        val = m.captures[2]
        attributes[att] = val
    end

    root, _ = readnode(text, pos+1)

    return Xdoc(attributes, root)
end

function printnode(f::IOStream, node::Xnode, level::Int)
    tab = "   "
    # Print header
    print(f, tab^level, "<", node.name)
    for (att,val) in node.attributes
        print(f, " $att=\"$val\"")
    end
    println(f, ">")

    # Print children or content
    if length(node.children)==0
        println(f, tab^(level+1), node.content)
    else
        for child in node.children
            printnode(f, child, level+1)
        end
    end

    # Print closing tag
    println(f, tab^level, "</", node.name, ">")
end

function save(doc::Xdoc, filename::String)
    # Open filename
    f = open(filename, "w")

    # Print prolog
    print(f, "<?xml")
    for (att,val) in doc.attributes
        print(f, " $att=\"$val\"")
    end
    println(f, "?>")

    # Print nodes
    printnode(f, doc.root, 0)

    close(f)
end
