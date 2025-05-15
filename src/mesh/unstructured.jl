# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function mesh_unstructured(geo::GeoModel; kwargs...)
    args      = checkargs([geo], kwargs, Mesh_Geo_params)
    size      = args.size
    quadratic = args.quadratic
    recombine = args.recombine
    algorithm = args.algorithm
    quiet     = args.quiet

    if !quiet
        println("  Geometry:")
        nsurfs = length(geo.faces)
        nvols  = length(geo.volumes)
        @printf "  %5d surfaces\n" nsurfs
        nvols>0 && @printf "  %5d volumes\n" nvols
    end

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("t1")
    # if geo.size>0
    #     gmsh.option.setNumber("Mesh.CharacteristicLengthMin", geo.size)
    #     gmsh.option.setNumber("Mesh.CharacteristicLengthMax", geo.size)
    # end

    isvolumemesh = length(geo.volumes)>0

    # add points
    for p in geo.points
        sz = p.size==0 ? size : p.size
        gmsh.model.geo.addPoint(p.coord.x, p.coord.y, p.coord.z, sz, p.id)
    end

    # add lines
    for l in geo.edges
        if l isa Line
            p1 = l.points[1].id
            p2 = l.points[2].id
            gmsh.model.geo.addLine(p1, p2, l.id)
        elseif l isa Arc
            p1 = l.points[1].id
            pc = l.points[2].id
            p2 = l.points[3].id

            gmsh.model.geo.addCircleArc(p1, pc, p2, l.id)
        else
            error("Not implemented")
        end
    end

    # add loops
    for loop in geo.loops
        idxs = Int[]
        for dart in loop.darts
            idx = dart.forward ? dart.edge.id : -dart.edge.id
            push!(idxs, idx)
        end

        gmsh.model.geo.addCurveLoop(idxs, loop.id)
    end

    # add surfaces
    for surf in geo.faces
        lo_idxs = [ lo.id for lo in surf.loops ]
        if surf.flat
            gmsh.model.geo.addPlaneSurface(lo_idxs, surf.id) # plane surface
        else
            gmsh.model.geo.addSurfaceFilling(lo_idxs, surf.id) # filling surf
        end
    end
    surf_idxs = [ surf.id for surf in geo.faces ]

    # add volumes
    for vol in geo.volumes
        surf_idxs = [ spin.face.id for spin in vol.spins ]

        gmsh.model.geo.addSurfaceLoop(surf_idxs, vol.id)
        gmsh.model.geo.addVolume([vol.id], vol.id) # not considering volume holes
    end
    vol_idxs = [ vol.id for vol in geo.volumes ]

    gmsh.model.geo.synchronize() # only after geometry entities are defined
    
    for l in geo.edges
        # transfinite
        if l.n>0
            gmsh.model.mesh.set_transfinite_curve(l.id, l.n+1)
        end
    end
    
    # generate mesh
    if algorithm==:delaunay
        gmsh.option.setNumber("Mesh.Algorithm", 5)
    elseif algorithm==:frontal
        gmsh.option.setNumber("Mesh.Algorithm", 6)
    else
        error("Mesh: Wrong algorithm")
    end

    for s in geo.faces
        s.recombine && gmsh.model.mesh.set_recombine(2, s.id)
        s.transfinite && gmsh.model.mesh.set_transfinite_surface(s.id)
    end

    if !isvolumemesh
        tagset = Set([ surf.tag for surf in geo.faces ])
        tagsdict = Dict( tag=>i for (i,tag) in enumerate(tagset) )
        for (tag, gidx) in tagsdict
            surf_idxs = [ surf.id for surf in geo.faces if surf.tag==tag ]
            gmsh.model.addPhysicalGroup(2, surf_idxs, gidx) # ndim, entities, group_id
        end
    else
        tagset = Set([ vol.tag for vol in geo.volumes  ])
        tagsdict = Dict( tag=>i for (i,tag) in enumerate(tagset) )
        for (tag, gidx) in tagsdict
            vol_idxs = [ vol.id for vol in geo.volumes if vol.tag==tag]
            gmsh.model.addPhysicalGroup(3, vol_idxs, gidx) # ndim, entities, group_id
        end
    end

    # embed points
    for p in geo.points
        length(p.edges)==0 || continue

        # search in sub-paths
        found = false
        for spath in geo.subpaths
            if p in spath.path.points
                found = true
                break
            end
        end
        found && continue
    
        # search surfaces
        for s in geo.faces
            s.loops[1].flat || continue
            if inside(p, s.loops[1])
                gmsh.model.mesh.embed(0,[p.id],2,s.id)
            end
        end
    end

    # embed lines
    for l in geo.edges
        l isa Line || continue
        length(l.faces)==0 || continue
        for s in geo.faces
            s.loops[1].flat || continue
            if insidepolygon(l.points, getpoints(s.loops[1]))
                gmsh.model.mesh.embed(1,[l.id],2,s.id)
            end
        end
    end

    tempfile = "_temp.vtk"
    logfile = "_gmsh.log"
    try
        open(logfile, "w") do out
            redirect_stdout(out) do
                # gmsh.write("_temp.geo_unrolled")
                gmsh.model.mesh.generate(isvolumemesh ? 3 : 2)
                quadratic && gmsh.model.mesh.setOrder(2) # quadratic elements
                recombine && gmsh.model.mesh.recombine()
                gmsh.write(tempfile)
            end
        end
    catch err
        error("Error generating unstructured mesh.")
    end
    
    gmsh.finalize()
    mesh = Mesh(tempfile)
    rm(tempfile, force=true)
    rm(logfile, force=true)

    # flip elements
    for elem in mesh.elems
        isinverted(elem) && flip!(elem)
    end

    # set tags for elements
    if !isvolumemesh
        invtagsdict = Dict( i=>tag for (tag,i) in tagsdict )
        for elem in mesh.elems
            elem.tag = invtagsdict[ mesh.elem_data["CellEntityIds"][elem.id] ]
        end
    else
        invtagsdict = Dict( i=>tag for (tag,i) in tagsdict )
        for elem in mesh.elems
            elem.tag = invtagsdict[ mesh.elem_data["CellEntityIds"][elem.id] ]
        end
    end

    delete!(mesh.elem_data, "CellEntityIds")

    # set tag for nodes
    ptagdict = Dict()
    for p in geo.points
        p.tag!="" || continue
        ptagdict[p.coord] = p.tag
    end
    
    for node in mesh.nodes
        tag = get(ptagdict, node.coord, "")
        tag != "" || continue
        node.tag = tag
    end

    synchronize!(mesh, sortnodes=true)

    return mesh

end
