# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Domain


"""
`Domain(mesh, [filekey="out"])`

Creates an `Domain` object based on a Mesh object `mesh` and represents the geometric domain to be analysed by the finite element analysis.

**Fields**

`nodes`: An array of nodes

`elems`: An array of finite elements

`faces`: An array of `Face` objects containing the boundary faces

`edges`: An array of `Edge` objects containing all the boundary faces

`filekey` : An string object that is used as part of the filename of resulting analyses files

"""
mutable struct Domain
    ndim ::Int
    nodes::Array{Node,1}
    elems::Array{Element,1}
    faces::Array{Face,1}
    edges::Array{Edge,1}
    filekey::String

    loggers::Array{AbstractLogger,1}
    nincs::Integer
    nouts::Integer
    ndofs::Integer
    stage::Integer
    shared_data::SharedAnalysisData

    function Domain(;filekey::String="out")
        this = new()

        this.filekey  = filekey

        this.loggers = []
        this.nincs   = 0
        this.ndofs   = 0
        this.nouts   = 0
        this.stage   = 0
        return this
    end
end


function Domain(mesh::Mesh, matbinds::Union{MaterialBind, Array{MaterialBind,1}},
                loggers::Union{AbstractLogger,Array}=[];
                model_type::Symbol=:general, thickness::Real=1.0, filekey::String="out", verbose::Bool=true)

    # Fix arguments
    if typeof(matbinds)==MaterialBind   ; matbinds = [ matbinds ] end
    if typeof(loggers)<:AbstractLogger; loggers = [ loggers ] end

    dom  = Domain(filekey=filekey)

    # Shared analysis data
    ndim = dom.ndim = mesh.ndim
    dom.shared_data = SharedAnalysisData()
    dom.shared_data.ndim = ndim 
    dom.shared_data.model_type = model_type
    dom.shared_data.thickness = thickness
    dom.shared_data.t = 0.0

    if verbose
        print_with_color(:cyan, "Domain setup:", bold=true)
        println()
    end

    # Setting nodes
    dom.nodes = [ Node([p.x, p.y, p.z], tag=p.tag, id=i) for (i,p) in enumerate(mesh.points)]

    # Setting new elements
    verbose && print("  setting elements...\r")
    ncells    = length(mesh.cells)
    dom.elems = Array{Element,1}(ncells)
    Nips      = zeros(Int, ncells)       # list with number of ips per element
    Tips      = Array{TagType,1}(ncells)  # list with the ip tag per element
    Tips     .= ""
    for mb in matbinds
        cells = mesh.cells[mb.expr]
        if isempty(cells)
            warn("Domain: binding material model $(typeof(mb.mat)) to an empty list of cells (expr: $(mb.expr))")
        end
        etype = matching_elem_type(mb.mat)
        for cell in cells
            ty = cell.embedded ? embedded_elem_type(etype) : etype
            if matching_shape_class(ty) != cell.shape.class
                error("Domain: material model $(typeof(mb.mat)) cannot be used with shape $(cell.shape.name) (cell id: $(cell.id))\n")
            end

            conn = [ p.id for p in cell.points ]
            elem = new_element(ty, cell.shape, dom.nodes[conn], dom.shared_data, cell.tag)

            elem.id = cell.id
            elem.mat = mb.mat
            dom.elems[cell.id] = elem
            #Nips[elem.id] = mb.nips
            Nips[elem.id] = cell.nips
            #Tips[elem.id] = mb.iptag
            Tips[elem.id] = cell.iptag
        end
    end

    # Check if all elements have material defined
    undefined_elem_shapes = Set{String}()
    for i=1:ncells
        if !isassigned(dom.elems, i)
            push!(undefined_elem_shapes, mesh.cells[i].shape.name)
        end
    end
    if !isempty(undefined_elem_shapes)
        error("Domain: missing material definition to allocate elements with shape: $(join(undefined_elem_shapes, ", "))\n")
    end

    # Setting linked elements
    for cell in mesh.cells
        for lcell in cell.linked_cells
            push!(dom.elems[cell.id].linked_elems, dom.elems[lcell.id])
        end
    end

    # Setting faces
    dom.faces = Array{Face}(0)
    for (i,cell) in enumerate(mesh.faces)
        conn = [ p.id for p in cell.points ]
        face = Face(cell.shape, dom.nodes[conn], ndim, cell.tag)
        face.oelem = dom.elems[cell.ocell.id]
        face.id = i
        push!(dom.faces, face)
    end

    # Setting edges
    dom.edges = Array{Edge}(0)
    for (i,cell) in enumerate(mesh.edges)
        conn = [ p.id for p in cell.points ]
        edge = Edge(cell.shape, dom.nodes[conn], ndim, cell.tag)
        edge.oelem = dom.elems[cell.ocell.id]
        edge.id = i
        push!(dom.edges, edge)
    end

    # Finishing to configure elements
    ip_id = 0
    for elem in dom.elems
        elem_config_dofs(elem)               # dofs
        elem_config_ips(elem, Nips[elem.id]) # ips
        for ip in elem.ips # updating ip tags
            ip_id += 1
            ip.id = ip_id
            ip.tag = Tips[elem.id]
        end
    end

    # Initializing elements
    verbose && print("  initializing elements...\r")
    for elem in dom.elems
        elem_init(elem)
    end

    # Setting loggers
    for logger in loggers
        setup_logger!(dom, logger)
    end
    dom.loggers = loggers

    if verbose
        print("  ", ndim, "D domain $model_type model      \n")
        @printf "  %5d nodes\n" length(dom.nodes)
        @printf "  %5d elements\n" length(dom.elems)
        if ndim==2
            @printf "  %5d edges\n" length(dom.faces)
        else
            @printf "  %5d faces\n" length(dom.faces)
            @printf "  %5d edges\n" length(dom.edges)
        end
        @printf "  %5d materials\n" length(matbinds)
        @printf "  %5d loggers\n" length(dom.loggers)
        println("  done.")
    end

    return dom
end


# Function for setting loggers
"""
    setlogger!(domain, logger)

Register a new `logger` in `domain`.

    setlogger!(domain, loggers)

Register each logger from the array `loggers` in `domain`.

"""
function setlogger!(dom::Domain, logger::Union{AbstractLogger, Array{<:AbstractLogger,1}})
    loggers = typeof(logger)<:AbstractLogger ? [logger] : logger
    setup_logger!.(dom, loggers)
    dom.loggers = loggers
end


# Function for updating loggers
update_loggers!(domain::Domain) = update_logger!.(domain.loggers)


# Function to reset a domain
#= This is error prone specially with ips in loggers
A new Domain is preferable
function reset!(dom::Domain) 
    dom.nincs = 0
    dom.stage = 0

    # Reconfigure nodes and dofs
    for node in dom.nodes
        empty!(node.dofs)
        empty!(node.dofdict)
    end

    # Reconfigure elements, dofs and ips
    ip_id = 0
    for elem in dom.elems
        elem_config_dofs(elem)
        ip_tags = [ ip.tag for ip in elem.ips]
        elem_config_ips(elem, length(elem.ips)) # ips
        for (i,ip) in enumerate(elem.ips) # updating ip tags
            ip_id += 1
            ip.id = ip_id
            ip.tag = ip_tags[i]
        end
    end

    # Reconfigure loggers
    for logger in dom.loggers
        setup_logger!(dom, logger)
    end
end

function reset!(dom::Domain)
    dom.nincs = 0
    dom.stage = 0

    # Reconfigure nodes and dofs
    for node in dom.nodes
        reset!(node)
    end

    # Reconfigure elements, dofs and ips
    for elem in dom.elems
        reset!(elem)
    end

    # Reconfigure loggers
    for logger in dom.loggers
        reset!(logger)
    end

end

=#


function save(dom::Domain, filename::String; verbose=true, save_ips=false)
    # Saves the dom information in vtk format
    nnodes = length(dom.nodes)
    nelems = length(dom.elems)

    # points, cells and cell types
    points  = [ node.X[i] for node in dom.nodes, i in 1:3]
    cells   = [ [node.id for node in elem.nodes] for elem in dom.elems ]
    cell_tys= [ elem.shape.vtk_type for elem in dom.elems ]
    
    # Node and element data
    point_scalar_data = Dict()
    point_vector_data = Dict()
    cell_scalar_data  = Dict()

    # Get node and elem values
    #node_vals, node_labels, elem_vals, elem_labels = all_node_and_elem_vals(dom.nodes, dom.elems)
    node_vals, node_labels, elem_vals, elem_labels = node_and_elem_vals(dom)
    nncomps = length(node_labels)
    necomps = length(elem_labels)

    # Write vectors
    if :uy in node_labels
        ux_idx = findfirst(node_labels, :ux)
        uy_idx = findfirst(node_labels, :uy)
        uz_idx = findfirst(node_labels, :uz)
        if uz_idx>0
            U = node_vals[:, [ux_idx, uy_idx, uz_idx]]
        else
            U = hcat( node_vals[:, [ux_idx, uy_idx]], zeros(nnodes) )
        end
        point_vector_data["U"] = U
    end

    if :vy in node_labels
        vx_idx = findfirst(node_labels, :vx)
        vy_idx = findfirst(node_labels, :vy)
        vz_idx = findfirst(node_labels, :vz)
        if vz_idx>0
            V = node_vals[:, [vx_idx, vy_idx, vz_idx]]
        else
            V = hcat( node_vals[:, [vx_idx, vy_idx]], zeros(nnodes) )
        end
        point_vector_data["V"] = V
    end

    # Write scalars
    point_scalar_data["node-id"] = [ node.id for node in dom.nodes ]

    # Write nodal scalar data
    for i=1:nncomps
        field = string(node_labels[i])
        point_scalar_data[field] = node_vals[:,i]
    end

    # Write cell data
    cell_scalar_data["elem-id"] = [ elem.id for elem in dom.elems ]
    cell_scalar_data["cell-type"] = [ elem.shape.vtk_type for elem in dom.elems ]

    # Write elem scalar data
    for i=1:necomps
        field = string(elem_labels[i])
        cell_scalar_data[field] = elem_vals[:,i]
    end
    
    vtk_data = VTK_unstructured_grid("Amaru - Finite Element Structures and Tools",
                                     points, cells, cell_tys,
                                     point_scalar_data=point_scalar_data,
                                     point_vector_data=point_vector_data,
                                     cell_scalar_data=cell_scalar_data)

    # Save file
    save_vtk(vtk_data, filename)

    verbose && print_with_color(:green, "  file $filename written (Domain)\n")

    # save ip information as vtk
    if save_ips
        # List of ips
        ips = [ ip for elem in dom.elems for ip in elem.ips ]
        nips = length(ips)

        # points, cells and cell types
        points = [ ips.X[i] for ip in ips, i in 1:3]
        cells  = [ [ip.id] for ip in ips ]
        cell_tys = ones(Int, nips) # VTK_VERTEX

        # get values at ips
        table = DTable()
        for ip in ips
            push!(table, ip_vals(ip))
        end

        # Write ip scalar data
        point_scalar_data = Dict()
        for field in keys(table)
            point_scalar_data[string(field)] = table[field]
        end

        vtk_data = VTK_unstructured_grid("Amaru - Finite Element Structures and Tools",
                                         points, cells, cell_tys,
                                         point_scalar_data=point_scalar_data)

        # Save file
        ips_filename = splitext(filename)[1]*"-ip.vtk"
        save(vtk_data, ips_filename)

        verbose && print_with_color(:green, "  file $ips_filename written (Domain IPs)\n")
    end
end

function node_and_elem_vals(dom::Domain)
    # nodal values
    nnodes = length(dom.nodes)

    node_fields_set = Set{Symbol}()
    #node_fields_set = Set( key for node in dom.nodes for dof in node.dofs for key in keys(dof.vals) )
    for node in dom.nodes
        for dof in node.dofs
            union!(node_fields_set, keys(dof.vals))
        end
    end
    #node_fields_set = Set( key for node in dom.nodes for dof in node.dofs for key in keys(dof.vals) )
    node_fields_idx = OrderedDict( key=>i for (i,key) in enumerate(node_fields_set) )
    nfields = length(node_fields_set)
    NV = zeros(nnodes, nfields)
    for node in dom.nodes
        for dof in node.dofs
            for (field,val) in dof.vals
                NV[ node.id, node_fields_idx[field] ] = val
            end
        end
    end

    # add nodal values from patch recovery (solid elements) : regression + averaging
    V_rec, fields_rec = nodal_patch_recovery(dom)
    NV = [ NV V_rec ]
    node_fields = [ collect(node_fields_set); fields_rec ]

    # add nodal values from local recovery (joints) : extrapolation + averaging
    V_rec, fields_rec = nodal_local_recovery(dom)
    NV = [ NV V_rec ]
    node_fields = [ node_fields; fields_rec ]

    # element values
    nelems = length(dom.elems)
    all_elem_vals   = [ elem_vals(elem) for elem in dom.elems ]
    elem_fields_set = Set( key for elem in dom.elems for key in keys(all_elem_vals[elem.id]) )
    elem_fields_idx = OrderedDict( key=>i for (i,key) in enumerate(elem_fields_set) )
    nfields = length(elem_fields_set)
    EV = zeros(nelems, nfields)
    for elem in dom.elems
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


function nodal_patch_recovery(dom::Domain)
    ndim = dom.shared_data.ndim
    nnodes = length(dom.nodes)
    length(dom.faces)==0 && return zeros(nnodes,0), Symbol[]

    # get surface nodes
    bry_nodes_set = Set( node for face in dom.faces for node in face.nodes )

    # list for boundary nodes
    at_bound = falses(nnodes)
    for node in bry_nodes_set
        at_bound[node.id] = true
    end

    # generate patches for solid elements
    patches     = [ Element[] for i=1:nnodes ] # internal patches
    bry_patches = [ Element[] for i=1:nnodes ] # boundary patches
    for elem in dom.elems
        elem.shape.class != SOLID_SHAPE && continue
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
    orphan_nodes = [ node for node in dom.nodes if (!haspatch[node.id] && at_bound[node.id]) ]

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
    all_fields_set = Set{Symbol}()
    for elem in dom.elems
        if elem.shape.class==SOLID_SHAPE
            ips_vals = [ ip_state_vals(elem.mat, ip.data) for ip in elem.ips ]
            push!(all_ips_vals, ips_vals)
            union!(all_fields_set, keys(ips_vals[1]))
        else # skip data from non solid elements
            push!(all_ips_vals, [])
        end
    end

    # map field => index
    all_fields_idx = OrderedDict( key=>i for (i,key) in enumerate(all_fields_set) )
    nfields = length(all_fields_set)

    # matrices for all nodal values and repetitions
    V_vals =  zeros(Float64, nnodes, nfields)
    V_reps =  zeros(Int64  , nnodes, nfields)

    # patch recovery
    for patch in patches
        length(patch) == 0 && continue

        # list of fields
        fields = unique( key for elem in patch for key in keys(all_ips_vals[elem.id][1]) )

        last_subpatch  = [] # elements of a subpatch for a particular field
        invM = Array{Float64,2}(0,0)
        N    = Array{Int64,2}(0,0)
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
                    nterms = m>=10? 10 : m>=7? 7 : m>=4? 4 : 1
                else
                    nterms = m>=6? 6 : m>=4? 4 : m>=3? 3 : 1
                end

                # matrix M for regression
                M = Array{Float64,2}(m,nterms)
                for (i,ip) in enumerate(subpatch_ips)
                    x, y, z = ip.X
                    if ndim==3
                        M[i,:] .= reg_terms(x, y, z, nterms)
                    else
                        M[i,:] .= reg_terms(x, y, nterms)
                    end
                end
                invM = pinv(M)

                # find nodal values
                N = Array{Float64,2}(n,nterms)
                for (i,node) in enumerate(subpatch_nodes)
                    x, y, z = node.X
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

    # average values
    V_vals ./= V_reps
    V_vals[isnan.(V_vals)] = 0.0

    return V_vals, collect(all_fields_set)

end


function nodal_local_recovery(dom::Domain)
    # Recovers nodal values from non-solid elements as joints and joint1d elements
    # The element type should implement the elem_extrapolated_node_vals function

    ndim = dom.shared_data.ndim
    nnodes = length(dom.nodes)

    # all local data from elements
    all_node_vals  = Array{OrderedDict{Symbol,Array{Float64,1}}, 1}()
    all_fields_set = OrderedSet{Symbol}()
    rec_elements   = Array{Element, 1}()

    for elem in dom.elems
        elem.shape.class == SOLID_SHAPE && continue
        node_vals = elem_extrapolated_node_vals(elem)
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
    V_vals[isnan.(V_vals)] = 0.0

    return V_vals, collect(all_fields_set)
end


function Base.convert(::Type{VTK_unstructured_grid}, dom::Domain)

    # Saves the dom information in vtk format
    nnodes = length(dom.nodes)
    nelems = length(dom.elems)

    # points, cells and cell types
    points  = [ node.X[i] for node in dom.nodes, i in 1:3]
    cells   = [ [node.id for node in elem.nodes] for elem in dom.elems ]
    cell_tys= [ elem.shape.vtk_type for elem in dom.elems ]
    
    # Node and element data
    point_scalar_data = Dict()
    point_vector_data = Dict()
    cell_scalar_data  = Dict()

    # Get node and elem values
    #node_vals, node_labels, elem_vals, elem_labels = all_node_and_elem_vals(dom.nodes, dom.elems)
    node_vals, node_labels, elem_vals, elem_labels = node_and_elem_vals(dom)
    nncomps = length(node_labels)
    necomps = length(elem_labels)

    # Write vectors
    if :uy in node_labels
        ux_idx = findfirst(node_labels, :ux)
        uy_idx = findfirst(node_labels, :uy)
        uz_idx = findfirst(node_labels, :uz)
        if uz_idx>0
            U = node_vals[:, [ux_idx, uy_idx, uz_idx]]
        else
            U = hcat( node_vals[:, [ux_idx, uy_idx]], zeros(nnodes) )
        end
        point_vector_data["U"] = U
    end

    # Write scalars
    point_scalar_data["node-id"] = [ node.id for node in dom.nodes ]

    # Write nodal scalar data
    for i=1:nncomps
        field = string(node_labels[i])
        point_scalar_data[field] = node_vals[:,i]
    end

    # Write cell data
    cell_scalar_data["elem-id"] = [ elem.id for elem in dom.elems ]
    cell_scalar_data["cell-type"] = [ elem.shape.vtk_type for elem in dom.elems ]

    # Write elem scalar data
    for i=1:necomps
        field = string(elem_labels[i])
        cell_scalar_data[field] = elem_vals[:,i]
    end
    
    return VTK_unstructured_grid("Amaru - Finite Element Structures and Tools",
                                 points, cells, cell_tys,
                                 point_scalar_data=point_scalar_data,
                                 point_vector_data=point_vector_data,
                                 cell_scalar_data=cell_scalar_data)
end


import FemMesh.mplot
function mplot(dom::Domain; args...)
    ugrid = convert(VTK_unstructured_grid, dom)
    if dom.ndim==3 # Plot the surface only
        nodes     = get_nodes(dom.faces)
        node_ids  = [ node.id for node in nodes ]
        id_dict   = Dict{Int,Int}( node.id=>i for (i,node) in enumerate(nodes) )
        oelem_ids = [ face.oelem.id for face in dom.faces ]

        ugrid.points = [ node.X[i] for node in nodes, i=1:3 ]
        ugrid.cells  = [ [ id_dict[node.id] for node in face.nodes] for face in dom.faces ]
        ugrid.cell_types = [ face.shape.vtk_type for face in dom.faces ]

        # update data
        for (field, data) in ugrid.point_scalar_data
            ugrid.point_scalar_data[field] = data[node_ids]
        end
        for (field, data) in ugrid.cell_scalar_data
            ugrid.cell_scalar_data[field] = data[oelem_ids]
        end
        for (field, data) in ugrid.point_vector_data
            ugrid.point_vector_data[field] = data[node_ids, :]
        end
    end
    mplot(ugrid; args...)
end
