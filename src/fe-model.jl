# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    FEModel

A type that represents a finite element model.

# Fields
$(FIELDS)
"""
mutable struct FEModel<:AbstractDomain
    nodes::Array{Node,1}
    elems::Array{Element,1}
    faces::Array{CellFace,1}
    edges::Array{CellEdge,1}
    thickness::Float64
    ndofs   ::Integer
    ctx     ::Context

    # Data
    node_data::OrderedDict{String,Array}
    elem_data::OrderedDict{String,Array}

    # Auxiliary data
    _elempartition::ElemPartition

    function FEModel()
        this = new()

        this.nodes = []
        this.elems = []
        this.faces = []
        this.edges = []
        
        this.thickness = 1.0
        this.ndofs   = 0

        this._elempartition = ElemPartition()
        this.node_data      = OrderedDict()
        this.elem_data      = OrderedDict()
        return this
    end
end

# const FEModel = FEModel


"""
    FEModel(mesh, mats, options...)

Uses a mesh and a list of meterial especifications to construct a finite element `FEModel`.

# Arguments

`mesh` : A finite element mesh

`mats` : Material definitions given as an array of pairs ( tag or location => constitutive model instance )

# Keyword arguments

`stressmodel`
`thickness`
`quiet = false` : If true, provides information of the model construction

"""
function FEModel(
    mesh    :: Mesh,
    matbinds:: Vector{<:Pair},
    ctx :: Context;
    thickness :: Real = 1.0,
    quiet   :: Bool = false,
    # outdir  :: String=".",
    # outkey  :: String="out"
)

    ndim = mesh.ctx.ndim
    @assert ndim>0

    # FEModel and environment data
    model  = FEModel()
    model.thickness = thickness
    model.ctx = ctx
    ctx.ndim = ndim
    # ctx.t = 0.0
    # ctx.outdir = rstrip(outdir, ['/', '\\'])
    # ctx.outkey = outkey

    # check if outdir exists
    # isdir(ctx.outdir) || mkdir(ctx.outdir)

    quiet || printstyled("FE model setup\n", bold=true, color=:cyan)

    # Setting nodes
    model.nodes = copy.(mesh.nodes)

    # Setting new elements
    ncells      = length(mesh.elems)
    model.elems = Array{Element,1}(undef, ncells)
    for matbind in matbinds
        matbind isa Pair{<:Any, <:Pair{DataType, <:Pair{DataType, <:NamedTuple}}} || throw(AmaruException("FEModel: Assigments of element and material models should be specified as: filter => Element => Material => properties"))
        
        filter    = matbind.first
        elem_type = matbind.second.first
        mat_type  = matbind.second.second.first
        args      = matbind.second.second.second
        # Check if material is compatible with the element
        comp_elem_types = [ argtps[2] for argtps in typeofargs(compat_state_type) if length(argtps)>1 && mat_type isa argtps[1] ]

        if !any(isa.(elem_type, comp_elem_types))
            comp_mat_types  = [ argtps[1].parameters[1] for argtps in typeofargs(compat_state_type) if length(argtps)>1 && typeof(argtps[1])!=UnionAll && argtps[2].parameters[1]==elem_type ]

            message = "FEModel: Material model $(mat_type) is not compatible with Element $(elem_type) \n\
            Compatible elements for material $(mat_type): $(join(comp_elem_types, ", ", " and ")) \n\
            Compatible materials for element $(elem_type): $(join(comp_mat_types, ", ", " and "))"
            message = replace(message, r"Amaru\." => "")
            throw(AmaruException(message))
        end


        cells = mesh.elems[filter]
        if !(cells isa Array)
            cells = [ cells ]
        end
        if isempty(cells)
            warn("FEModel: No elements found for expression $(repr(filter))")
        end

        mat   = mat_type(;args...)
        if hasmethod(compat_elem_props, Tuple{Type{elem_type}})
            props = compat_elem_props(elem_type)(;args...)
        else
            props = (;)
        end

        for cell in cells
            if cell.embedded
                elem_t = embedded_type(elem_type)
            else
                elem_t = elem_type
            end

            if compat_shape_family(elem_t) != cell.shape.family
                error("FEModel: Material model $(typeof(mat)) cannot be used with elements for expr. $(repr(filter)) with shape $(cell.shape.name)\n")
            end

            elem = elem_t()
            
            conn = [ p.id for p in cell.nodes ]

            elem.id     = cell.id
            elem.ctx    = ctx
            elem.shape  = cell.shape
            elem.tag    = cell.tag
            elem.nodes  = model.nodes[conn]
            elem.mat    = mat
            elem.active = true
            elem.ips    = [] # jet to be set
            elem.linked_elems = [] # jet to be set
            props!=(;) && (elem.props=props)

            model.elems[cell.id] = elem 
        end
    end

    # Check if all elements have material defined
    undefined_elem_shapes = Set{String}()
    for i in 1:ncells
        if !isassigned(model.elems, i)
            push!(undefined_elem_shapes, mesh.elems[i].shape.name)
        end
    end
    if !isempty(undefined_elem_shapes)
        error("FEModel: missing material definition to allocate elements with shape: $(join(undefined_elem_shapes, ", "))\n")
    end

    # Setting linked elements
    for cell in mesh.elems
        for lcell in cell.linked_elems
            push!(model.elems[cell.id].linked_elems, model.elems[lcell.id])
        end
    end

    # # Setting element references in nodes
    # for elem in mesh.elems
    #     for node in elem.nodes
    #         push!(node.elems, elem)
    #     end
    # end

    # Setting faces
    model.faces = CellFace[]
    for (i,cell) in enumerate(mesh.faces)
        conn = [ p.id for p in cell.nodes ]
        face = CellFace(cell.shape, model.nodes[conn], tag=cell.tag)
        if cell.owner!==nothing
            face.owner = model.elems[cell.owner.id]
        else
            error("FEModel: face $(i) has no owner element")
        end
        face.id = i
        push!(model.faces, face)
    end

    # Setting edges
    model.edges = CellEdge[]
    for (i,cell) in enumerate(mesh.edges)
        conn = [ p.id for p in cell.nodes ]
        edge = CellEdge(cell.shape, model.nodes[conn], tag=cell.tag)
        if cell.owner!==nothing
            edge.owner = model.elems[cell.owner.id]
        end
        edge.id = i
        push!(model.edges, edge)
    end

    # Finishing to configure elements
    ip_id = 0
    for elem in model.elems
        elem_config_dofs(elem)  # dofs
        setquadrature!(elem)    # ips
        for ip in elem.ips      # ip ids
            ip_id += 1
            ip.id = ip_id
        end
    end

    # Initializing elements
    for elem in model.elems
        elem_init(elem)
    end

    if !quiet
        # print("  ", "$stressmodel model      \n")
        @printf "  %5d nodes\n" length(model.nodes)
        @printf "  %5d elements\n" length(model.elems)
    end

    if !quiet
        if model.ctx.ndim==2
            @printf "  %5d edges\n" length(model.faces)
        else
            @printf "  %5d faces\n" length(model.faces)
            @printf "  %5d edges\n" length(model.edges)
        end
        @printf "  %5d materials\n" length(matbinds)
        # @printf "  %5d loggers\n" length(model.loggers)
    end

    # Setting data
    model.node_data["node-id"] = copy(mesh.node_data["node-id"])
    model.elem_data["elem-id"] = copy(mesh.elem_data["elem-id"])
    model.elem_data["cell-type"] = copy(mesh.elem_data["cell-type"])

    return model
end


function FEModel(elems::Array{<:Element,1})
    model            = FEModel()
    ctx              = elems[1].ctx
    model.ctx        = ctx
    model.ctx.ndim   = ctx.ndim

    # Copying nodes
    nodesset = OrderedSet(node for elem in elems for node in elem.nodes)
    model.nodes = copy.(collect(Node, nodesset))
    nodes    = model.nodes

    # Map for nodes
    nodemap = zeros(Int, maximum(node.id for node in nodes))
    for (i,node) in enumerate(nodes)
        nodemap[node.id] = i 
        model.nodes[i].id = i
    end

    # Map for elements
    elemmap = zeros(Int, maximum(elem.id for elem in elems))
    for (i,elem) in enumerate(elems)
        elemmap[elem.id] = i 
    end

    # Get ndim
    ndim = 1
    for node in nodes
        node.coord.y != 0.0 && (ndim=2)
        node.coord.z != 0.0 && (ndim=3; break)
    end
    model.ctx.ndim = ndim

    # Setting elements
    for (i,elem) in enumerate(elems)
        nodeidxs = [ nodemap[node.id] for node in elem.nodes]
        elemnodes = nodes[nodeidxs]
        newelem = new_element(typeof(elem), elem.shape, elemnodes, elem.tag, model.ctx)
        newelem.id = i
        newelem.mat = elem.mat
        newelem.props = elem.props
        push!(model.elems, newelem)
    end
    
    # Setting linked elements
    missing_linked = false
    for (i,elem) in enumerate(elems)
        for lelem in elem.linked_elems
            idx = elemmap[lelem.id]
            idx==0 && (missing_linked=true; continue)
            push!(model.elems[i].linked_elems, model.elems[idx])
        end
    end
    missing_linked && error("FEModel: Missing linked elements while generating a model from a list of elements.")

    # Setting quadrature
    ip_id = 0
    for (newelem, elem) in zip(model.elems, elems)
        setquadrature!(newelem, length(elem.ips))   # ips
        for ip in newelem.ips      # ip ids
            ip_id += 1
            ip.id = ip_id
        end
        elem_init(newelem)
    end

    # Copying states
    for (newelem, elem) in zip(model.elems, elems)
        for i in 1:length(elem.ips)
            copyto!(newelem.ips[i].state, elem.ips[i].state)
        end
    end

    # Setting faces and edges
    model.faces = get_outer_facets(model.elems)
    model.edges = getedges(model.faces)

    # Setting data
    update_output_data!(model)

    return model

end


function update_output_data!(model::FEModel)
    # Updates data arrays in the model
    model.node_data = OrderedDict()
    model.elem_data = OrderedDict()

    # Ids and cell type
    model.node_data["node-id"] = collect(1:length(model.nodes))
    model.elem_data["elem-id"]  = collect(1:length(model.elems))
    model.elem_data["cell-type"] = [ Int(elem.shape.vtk_type) for elem in model.elems ]

    # Nodal values
    nnodes = length(model.nodes)

    # get node field symbols
    node_fields_set = OrderedSet{Symbol}()
    for node in model.nodes
        for dof in node.dofs
            union!(node_fields_set, keys(dof.vals))
        end
    end
    node_fields = collect(node_fields_set)

    # Generate empty lists
    for field in node_fields
        model.node_data[string(field)] = zeros(nnodes)
    end

    # Fill dof values
    for node in model.nodes
        for dof in node.dofs
            for (field,val) in dof.vals
                model.node_data[string(field)][node.id] = val
            end
        end
    end

    # add nodal values from patch recovery (solid elements) : regression + averaging
    V_rec, fields_rec = nodal_patch_recovery(model)
    for (i,field) in enumerate(fields_rec)
        model.node_data[string(field)] = V_rec[:,i]
    end
    append!(node_fields, fields_rec)

    # add nodal values from local recovery (joints) : extrapolation + averaging
    V_rec, fields_rec = nodal_local_recovery(model)
    for (i,field) in enumerate(fields_rec)
        model.node_data[string(field)] = V_rec[:,i]
    end
    append!(node_fields, fields_rec)

    # Nodal vector values
    if :ux in node_fields
        if :uz in node_fields
            model.node_data["U"] = [ model.node_data["ux"] model.node_data["uy"] model.node_data["uz"] ]
        elseif :uy in node_fields
            model.node_data["U"] = [ model.node_data["ux"] model.node_data["uy"] zeros(nnodes) ]
        else
            model.node_data["U"] = [ model.node_data["ux"] zeros(nnodes) zeros(nnodes) ]
        end
    end

    if :vx in node_fields
        if :vz in node_fields
            model.node_data["V"] = [ model.node_data["vx"] model.node_data["vy"] model.node_data["vz"] ]
        elseif :vy in node_fields
            model.node_data["V"] = [ model.node_data["vx"] model.node_data["vy"] zeros(nnodes) ]
        else
            model.node_data["V"] = [ model.node_data["vx"] zeros(nnodes) zeros(nnodes) ]
        end
    end

    # Element values
    nelems = length(model.elems)
    all_elem_vals   = [ elem_vals(elem) for elem in model.elems ]
    elem_fields_set = Set( key for elem in model.elems for key in keys(all_elem_vals[elem.id]) )
    elem_fields     = collect(elem_fields_set)

    # generate empty lists
    for field in elem_fields
        model.elem_data[string(field)] = zeros(nelems)
    end

    # fill elem values
    for elem in model.elems
        for (field,val) in all_elem_vals[elem.id]
            model.elem_data[string(field)][elem.id] = val
        end
    end

end


function get_node_and_elem_vals(model::FEModel)
    # Return symbols and values for nodes and elements
    # Note: nodal ids must be numbered starting from 1

    # nodal values
    nnodes = length(model.nodes)

    # get node field symbols
    node_fields_set = Set{Symbol}()
    for node in model.nodes
        for dof in node.dofs
            union!(node_fields_set, keys(dof.vals))
        end
    end

    # get node field values
    node_fields_idx = OrderedDict( key=>i for (i,key) in enumerate(node_fields_set) )
    nfields = length(node_fields_set)
    NV = zeros(nnodes, nfields)
    for node in model.nodes
        for dof in node.dofs
            for (field,val) in dof.vals
                NV[ node.id, node_fields_idx[field] ] = val
            end
        end
    end

    # add nodal values from patch recovery (solid elements) : regression + averaging
    V_rec, fields_rec = nodal_patch_recovery(model)
    NV = [ NV V_rec ]
    node_fields = [ collect(node_fields_set); fields_rec ]

    # add nodal values from local recovery (joints) : extrapolation + averaging
    V_rec, fields_rec = nodal_local_recovery(model)
    NV = [ NV V_rec ]
    node_fields = [ node_fields; fields_rec ]

    # element values
    nelems = length(model.elems)
    all_elem_vals   = [ elem_vals(elem) for elem in model.elems ]
    elem_fields_set = Set( key for elem in model.elems for key in keys(all_elem_vals[elem.id]) )
    elem_fields_idx = OrderedDict( key=>i for (i,key) in enumerate(elem_fields_set) )
    nfields = length(elem_fields_set)
    EV = zeros(nelems, nfields)
    for elem in model.elems
        for (field,val) in all_elem_vals[elem.id]
            EV[ elem.id, elem_fields_idx[field] ] = val
        end
    end

    elem_fields = collect(elem_fields_set)
    return NV, node_fields, EV, elem_fields

end


function reg_terms(x::Float64, y::Float64, nterms::Int64)
    nterms==6 && return ( 1.0, x, y, x*y, x^2, y^2 )
    nterms==4 && return ( 1.0, x, y, x*y )
    nterms==3 && return ( 1.0, x, y )
    return (1.0,)
end


function reg_terms(x::Float64, y::Float64, z::Float64, nterms::Int64)
    nterms==10 && return ( 1.0, x, y, z, x*y, y*z, x*z, x^2, y^2, z^2 )
    nterms==7  && return ( 1.0, x, y, z, x*y, y*z, x*z )
    nterms==4  && return ( 1.0, x, y, z )
    return (1.0,)
end


function nodal_patch_recovery(model::FEModel)
    # Note: nodal ids must be numbered starting from 1

    ndim = model.ctx.ndim
    nnodes = length(model.nodes)
    length(model.faces)==0 && return zeros(nnodes,0), Symbol[]

    # get node field symbols (to skip if later found at ips)
    node_fields_set = OrderedSet{Symbol}()
    for node in model.nodes
        for dof in node.dofs
            union!(node_fields_set, keys(dof.vals))
        end
    end
    node_fields = collect(node_fields_set)

    # get surface nodes
    bry_nodes_set = Set( node for face in model.faces for node in face.nodes )

    # list for boundary nodes
    at_bound = falses(nnodes)
    for node in bry_nodes_set
        at_bound[node.id] = true
    end

    # generate patches for solid elements
    patches     = [ Element[] for i in 1:nnodes ] # internal patches
    bry_patches = [ Element[] for i in 1:nnodes ] # boundary patches
    for elem in model.elems
        elem.shape.family != BULKCELL && continue
        for node in elem.nodes[1:elem.shape.basic_shape.npoints] # only at corners
            if at_bound[node.id]
                push!(bry_patches[node.id], elem)
            else
                push!(patches[node.id], elem)
            end
        end
    end

    # check if nodes are in at least one internal patch
    npatches = zeros(Int, nnodes)
    haspatch = falses(nnodes)
    for patch in patches
        for elem in patch
            for node in elem.nodes
                haspatch[node.id] = true
            end
        end
    end
    orphan_nodes = [ node for node in model.nodes if (!haspatch[node.id] && at_bound[node.id]) ]

    # add border patches avoiding patches with few elements
    if length(orphan_nodes)>0
        for n=3:-1:1
            # add orphan_nodes patches if patch has more than n elems
            for node in orphan_nodes
                patch = bry_patches[node.id]
                length(patch)>=n && ( patches[node.id] = patch )
            end

            # check for orphan nodes
            for node in orphan_nodes
                patch = patches[node.id]
                for elem in patch
                    for node in elem.nodes
                        haspatch[node.id] = true
                    end
                end
            end
            orphan_nodes = [ node for node in orphan_nodes if !haspatch[node.id] ]
            length(orphan_nodes)==0 && break
        end
    end

    # all data from ips per element and field names
    all_ips_vals   = Array{Array{OrderedDict{Symbol,Float64}},1}()
    all_fields_set = OrderedSet{Symbol}()
    for elem in model.elems
        if elem.shape.family==BULKCELL
            ips_vals = [ ip_state_vals(elem.mat, ip.state) for ip in elem.ips ]
            push!(all_ips_vals, ips_vals)
            union!(all_fields_set, keys(ips_vals[1]))

        else # skip data from non solid elements
            push!(all_ips_vals, [])
        end
    end

    # map field => index
    all_fields_idx = OrderedDict( key=>i for (i,key) in enumerate(all_fields_set) )
    nfields = length(all_fields_set)

    nfields==0 && return zeros(Float64, nnodes, nfields), Symbol[]

    @withthreads begin
        # matrices for all nodal values and repetitions
        V_vals =  zeros(Float64, nnodes, nfields)
        V_reps =  zeros(Int64  , nnodes, nfields)

        # patch recovery
        for patch in patches

            length(patch) == 0 && continue

            # list of fields
            fields = unique( key for elem in patch for key in keys(all_ips_vals[elem.id][1]) )
            setdiff!(fields, node_fields) # remove fields already at nodes if any

            last_subpatch  = [] # elements of a subpatch for a particular field
            invM = Array{Float64,2}(undef,0,0)
            
            local N
            subpatch_ips = Ip[]
            subpatch_nodes = Node[]
            nterms = 0

            for field in fields
                # find subpatch for current field
                subpatch = [ elem for elem in patch if haskey(all_ips_vals[elem.id][1],field) ]

                # get subpatch data
                if subpatch != last_subpatch
                    subpatch_ips   = [ ip for elem in subpatch for ip in elem.ips ]
                    subpatch_nodes = unique( node for elem in subpatch for node in elem.nodes )
                    last_subpatch  = subpatch

                    m = length(subpatch_ips)
                    n = length(subpatch_nodes)

                    # number of polynomial terms
                    if ndim==3
                        # there were problems when using 10 terms
                        nterms = m>=7 ? 7 : m>=4 ? 4 : 1
                    else
                        nterms = m>=6 ? 6 : m>=4 ? 4 : m>=3 ? 3 : 1
                    end

                    # matrix M for regression
                    M = Array{Float64,2}(undef,m,nterms)
                    for (i,ip) in enumerate(subpatch_ips)
                        x, y, z = ip.coord
                        if ndim==3
                            M[i,:] .= reg_terms(x, y, z, nterms)
                        else
                            M[i,:] .= reg_terms(x, y, nterms)
                        end
                    end
                    invM = pinv(M)

                    # find nodal values
                    N = Array{Float64,2}(undef,n,nterms)
                    for (i,node) in enumerate(subpatch_nodes)
                        x, y, z = node.coord
                        if ndim==3
                            N[i,:] .= reg_terms(x, y, z, nterms)
                        else
                            N[i,:] .= reg_terms(x, y, nterms)
                        end
                    end
                end

                # values at ips
                W = Float64[ dict[field] for elem in subpatch for dict in all_ips_vals[elem.id] ]
                # coefficients vector from regression polynomial
                A = invM*W
                # values at nodes
                V = N*A

                # saving for later averaging
                field_idx = all_fields_idx[field]
                for (i,node) in enumerate(subpatch_nodes)
                    V_vals[node.id, field_idx] += V[i]
                    V_reps[node.id, field_idx] += 1
                end

            end
        end

    end

    # average values
    V_vals ./= V_reps
    V_vals[isnan.(V_vals)] .= 0.0

    return V_vals, collect(all_fields_set)

end


function nodal_local_recovery(model::FEModel)
    # Recovers nodal values from non-solid elements such as joints and joint1d elements
    # The element type should implement the elem_recover_nodal_values function
    # Note: nodal ids in the fe model must be numbered starting from 1

    nnodes = length(model.nodes)

    # all local data from elements
    all_node_vals  = Array{OrderedDict{Symbol,Array{Float64,1}}, 1}()
    all_fields_set = OrderedSet{Symbol}()
    rec_elements   = Array{Element, 1}()

    for elem in model.elems
        elem.shape.family == BULKCELL && continue
        node_vals = elem_recover_nodal_values(elem)
        length(node_vals) == 0 && continue

        push!(rec_elements, elem)
        push!(all_node_vals, node_vals)
        union!(all_fields_set, keys(node_vals))
    end

    # map field => index
    all_fields_idx = OrderedDict( key=>i for (i,key) in enumerate(all_fields_set) )
    nfields = length(all_fields_set)

    # matrices for all nodal values and repetitions
    V_vals =  zeros(Float64, nnodes, nfields)
    V_reps =  zeros(Int64  , nnodes, nfields)

    length(rec_elements) == 0 && return V_vals, collect(all_fields_set)

    # local recovery
    for (i,elem) in enumerate(rec_elements)
        node_vals = all_node_vals[i]
        row_idxs  = [ node.id for node in elem.nodes ]

        for (field, vals) in node_vals

            idx = all_fields_idx[field]
            V_vals[row_idxs, idx] .+= vals
            V_reps[row_idxs, idx] .+= 1
        end
    end

    # average values
    V_vals ./= V_reps
    V_vals[isnan.(V_vals)] .= 0.0

    return V_vals, collect(all_fields_set)
end


  