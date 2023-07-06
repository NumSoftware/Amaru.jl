# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function save_xml(model::Model, filename::String)
    io = IOBuffer()
    nnodes = length(model.nodes)
    nelems = length(model.elems)

    env = model.env
    attributes=OrderedDict{String,String}()

    fields = fieldnames(typeof(env))

    for fld in fields
        fld in (:matparams,:cinc,:cout) && continue
        value = getfield(env, fld)
        value == 0 && continue
        attributes[string(fld)] = string(value)
    end

    # Model
    root = Xnode("Model", attributes=attributes)

    # Materials
    xmats = Xnode("Materials")
    mat_dict = OrderedDict{UInt, MatParams}()
    for elem in model.elems
        hs = hash(elem.matparams)
        haskey(mat_dict, hs) && continue
        mat_dict[hs] = elem.matparams
    end

    mat_idx_dict = OrderedDict{UInt, Int}()
    i = 0
    for (k,v) in mat_dict
        push!(xmats.children, Xnode(v))
        i += 1
        mat_idx_dict[k] = i
    end
    push!(root.children, xmats)

    # Nodes
    xnodes = Xnode("Nodes")
    for node in model.nodes
        atts = OrderedDict(
                           "id"=>string(node.id),
                           "tag"=>string(node.tag),
                           "coord"=>"$(node.coord.x),$(node.coord.y),$(node.coord.z)",
                          )
        xnode = Xnode("Node", attributes=atts)
        for dof in node.dofs
            atts = OrderedDict(
                               "name"=>string(dof.name),
                               "natname"=>string(dof.natname),
                               "eq_id"=>string(dof.eq_id),
                               "prescribed"=>string(dof.prescribed),
                               "keys"=>join(keys(dof.vals), ","),
                               "vals"=>join(values(dof.vals), ","),
                              )
            xdof = Xnode("Dof", attributes=atts)
            push!(xnode.children, xdof)
        end
        push!(xnodes.children, xnode)
    end
    push!(root.children, xnodes)

    # Elements
    xelems = Xnode("Elements")
    for elem in model.elems
        atts = OrderedDict(
                           "id"=>string(elem.id),
                           "shape"=>string(elem.shape.name),
                           "tag"=>string(elem.tag),
                           "material"=>string(mat_idx_dict[hash(elem.matparams)]),
                           "active"=>string(elem.active),
                           "nodes"=>join((n.id for n in elem.nodes), ","),
                           "linked_elems"=>join((e.id for e in elem.linked_elems), ","),
                          )
        elemname = split(string(typeof(elem)), ".")[end]
        xelem = Xnode(elemname, attributes=atts) 
        for ip in elem.ips
            atts = OrderedDict(
                               "id"=>string(ip.id),
                               "tag"=>string(ip.tag),
                              )
            fields = fieldnames(typeof(ip.state))
            keys = []
            vals = []
            for fld in fields
                fld in (:env,) && continue
                push!(keys, fld)
                push!(vals, getfield(ip.state, fld))
            end
            atts["keys"] = join(keys, ",")
            atts["vals"] = join(vals, ",")
            xip = Xnode("Ip", attributes=atts)
            push!(xelem.children, xip)
        end
        push!(xelems.children, xelem)
    end
    push!(root.children, xelems)

    # NodeData
    xnodedata = Xnode("NodeData")
    for (field,D) in model.node_data
        isempty(D) && continue
        isfloat = eltype(D)<:AbstractFloat
        dtype = isfloat ? "Float64" : "Int32"
        ncomps = size(D,2)
        xdata = Xnode("DataArray", attributes=OrderedDict("name"=>string(field), "type"=>dtype, "ncomps"=>string(ncomps)))
        for i in 1:nnodes
            for j in 1:ncomps
                if isfloat
                    @printf io "%17.7e" Float32(D[i,j])
                else
                    print(io, D[i,j], "  ")
                end
            end
        end
        xdata.content = String(take!(io))
        push!(xnodedata.children, xdata)
    end
    push!(root.children, xnodedata)

    # ElemData
    xelemdata = Xnode("ElemData")
    for (field,D) in model.elem_data
        isempty(D) && continue
        isfloat = eltype(D)<:AbstractFloat
        dtype = isfloat ? "Float64" : "Int32"
        ncomps = size(D,2)
        xdata = Xnode("DataArray", attributes=OrderedDict("name"=>string(field), "type"=>dtype, "ncomps"=>string(ncomps)))
        for i in 1:nelems
            for j in 1:ncomps
                if isfloat
                    @printf io "%20.10e" Float32(D[i,j])
                else
                    print(io, D[i,j], "  ")
                end
            end
        end
        xdata.content = String(take!(io))
        push!(xelemdata.children, xdata)
    end
    push!(root.children, xelemdata)

    fileatts = OrderedDict("version"=>"1.0", "encoding"=>"utf-8")
    doc = Xdoc(fileatts, root)
    save(doc, filename)

end

"""
    save(domain, filename)

Saves a domain object into a file. Available formats are vtu, vtk and xml.
"""
function save(domain::Model, filename::String; quiet=false)

    formats = (".vtk", ".vtu", ".xml")
    _, format = splitext(filename)
    format in formats || error("save: Cannot save Model to $filename. Available formats are $formats.")

    if format==".xml"; 
        save_xml(domain, filename)
        quiet || printstyled( "  file $filename written \e[K \n", color=:cyan)
    else
        invoke(save, Tuple{AbstractDomain,String}, domain, filename, quiet=quiet)
    end
end


function save(elems::Array{<:Element,1}, filename::String; quiet=false)
    save(Model(elems), filename, quiet=quiet)
end


function setfields!(obj, keys, vals; exclude::Tuple{Vararg{Symbol}}=())
    for (k,v) in zip(keys,vals)
        field = Symbol(k)
        field in exclude && continue
        ty = fieldtype(typeof(obj), field)
        if ty==Bool
            setfield!(obj, field, parse(Bool, v))
        elseif ty<:Integer
            setfield!(obj, field, parse(Int, v))
        elseif ty<:Real
            setfield!(obj, field, parse(Float64, v))
        elseif ty==Symbol
            setfield!(obj, field, Symbol(v))
        elseif ty<:AbstractString
            setfield!(obj, field, v)
        elseif ty==Vec3
            setfield!(obj, field, Vec3(parse.(Float64, split(v, ","))))
        elseif ty<:AbstractArray
            setfield!(obj, field, parse.(Float64, split(v, (',','[',']'), keepempty=false)))
        else
            # Avoid Meta.parse as much as possible. It spends too much time.
            setfield!(obj, field, eval(Meta.parse(v)))
        end
    end
end

function setfields!(obj, dict; exclude::Tuple{Vararg{Symbol}}=())
    setfields!(obj, keys(dict), values(dict), exclude=exclude)
end


function Model(filename::String; quiet=false)
    suitable_formats = (".xml",)
    

    quiet || printstyled("Loading Model: filename $filename\n", bold=true, color=:cyan)

    basename, format = splitext(filename)

    format in suitable_formats || error("Model: cannot read format \"$format\". Suitable formats are $suitable_formats.")

    domain = Model()
    env = ModelEnv()

    quiet || printstyled("  loading xml file...\r", color=:cyan)
    xdoc = Xdoc(filename)
    xdomain = xdoc.root
    setfields!(env, xdomain.attributes)

    domain.env = env
    domain.ndim = env.ndim

    quiet || printstyled("  setting materials...\e[K\r", color=:cyan)
    materials = MatParams[]
    xmats = xdomain("Materials")
    for xmat in xmats.children
        T = eval(Symbol(xmat.name))
        matparams = ccall(:jl_new_struct_uninit, Any, (Any,), T)
        setfields!(matparams, xmat.attributes)

        push!(materials, matparams)
    end

    quiet || printstyled("  setting nodes...\e[K\r", color=:cyan)
    xnodes = xdomain("Nodes")
    for xnode in xnodes.children
        node = Node()
        setfields!(node, xnode.attributes)

        for xdof in xnode.children
            dof = ccall(:jl_new_struct_uninit, Any, (Any,), Dof)
            setfields!(dof, xdof.attributes, exclude=(:keys, :vals))
            keys = Symbol.(split(xdof.attributes["keys"], ","))
            vals = parse.(Float64, split(xdof.attributes["vals"], ","))
            dof.vals = OrderedDict(keys .=> vals)
            push!(node.dofs, dof)
        end
        
        push!(domain.nodes, node)
    end

    quiet || printstyled("  setting elements...\e[K\r", color=:cyan)
    xelems = xdomain("Elements")
    for xelem in xelems.children
        T = eval(Symbol(xelem.name))

        shape = eval(Symbol(xelem.attributes["shape"]))
        nodesidx = parse.(Int, split(xelem.attributes["nodes"],","))
        nodes = domain.nodes[nodesidx]

        elem = new_element(T, shape, nodes, "", domain.env)
        setfields!(elem, xelem.attributes, exclude=(:shape, :material, :nodes, :linked_elems))

        matidx = parse(Int,xelem.attributes["material"])
        elem.matparams = materials[matidx]
        nips = length(xelem.children)
      
        push!(domain.elems, elem)
    end

    # Setting linked elements
    quiet || printstyled("  setting linked elements...\e[K\r", color=:cyan)
    for (i,xelem) in enumerate(xelems.children)
        linked_str = xelem.attributes["linked_elems"]
        linked_str == "" && continue
        linked_idx = parse.(Int, split(linked_str,","))
        domain.elems[i].linked_elems = domain.elems[linked_idx]
    end

    # Quadrature and initialization
    quiet || printstyled("  setting integration points...\e[K\r", color=:cyan)
    for (i,xelem) in enumerate(xelems.children)
        elem = domain.elems[i]
        nips = length(xelem.children)

        setquadrature!(elem, nips)
        elem_init(elem)

        for (i,xip) in enumerate(xelem.children)
            ip = elem.ips[i]
            setfields!(ip, xip.attributes, exclude=(:keys, :vals))
            keys = Symbol.(split(xip.attributes["keys"], ",")) 
            vals = split(xip.attributes["vals"], r",(?! )")
            setfields!(elem.ips[i].state, keys, vals)
        end
        
        push!(domain.elems, elem)
    end

    domain.faces = get_surface(domain.elems)
    domain.edges = getedges(domain.faces)

    quiet || printstyled("  setting additional data...\e[K\r", color=:cyan)

    TYPES = Dict("Float32"=>Float32, "Float64"=>Float64, "Int32"=>Int32, "Int64"=>Int64)
    nnodes = length(domain.nodes)
    nelems = length(domain.elems)

    xnodedata = xdomain("NodeData")
    for arr in xnodedata.children
        ncomps = parse(Int, arr.attributes["ncomps"])
        dtype = TYPES[arr.attributes["type"]]
        label = arr.attributes["name"]
        if ncomps==1
            domain.node_data[label] = parse.(dtype, split(arr.content))
        else
            domain.node_data[label] = transpose(reshape(parse.(dtype, split(arr.content)), ncomps, nnodes))
        end
    end

    xelemdata = xdomain("ElemData")
    if xelemdata!==nothing
        for arr in xelemdata.children
            ncomps = parse(Int, arr.attributes["ncomps"])
            dtype = TYPES[arr.attributes["type"]]
            label = arr.attributes["name"]
            if ncomps==1
                domain.elem_data[label] = parse.(dtype, split(arr.content))
            else
                domain.elem_data[label] = transpose(reshape(parse.(dtype, split(arr.content)), ncomps, nelems))
            end
        end
    end

    quiet || printstyled( "  file $filename loaded \e[K \n", color=:cyan)

    return domain

end

function get_segment_data(model::Model, X1::Array{<:Real,1}, X2::Array{<:Real,1}, filename::String=""; npoints=50)
    data = model.node_data
    table = DataTable(["s"; collect(keys(data))])
    X1 = [X1; 0.0][1:3]
    X2 = [X2; 0.0][1:3]
    Δ = (X2-X1)/(npoints-1)
    Δs = norm(Δ)
    s1 = 0.0

    for i in 1:npoints
        X = X1 + Δ*(i-1)
        s = s1 + Δs*(i-1)
        cell = find_elem(X, model.elems, model._elempartition, 1e-7, Cell[])
        coords = getcoords(cell)
        R = inverse_map(cell.shape, coords, X)
        N = cell.shape.func(R)
        map = [ p.id for p in cell.nodes ]
        vals = [ s ]
        for (k,V) in data
            val = dot(V[map], N)
            push!(vals, val)
        end
        push!(table, vals)
    end

    filename != "" && save(table, filename)

    return table
end
