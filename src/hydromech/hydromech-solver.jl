# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export HydroAnalysis, HydromechAnalysis

HydromechAnalysis_params = [
    FunInfo(:HydromechAnalysis, "Hydromechanical analysis properties."),
    KwArgInfo(:stressmodel, "Stress model", :d3, values=(:planestress, :planestrain, :axisymmetric, :d3)),
    KwArgInfo(:thickness, "Thickness for 2d analyses", 1.0, cond=:(thickness>0)),
    KwArgInfo(:g, "Gravity acceleration", 0.0, cond=:(g>=0)),
    KwArgInfo(:gammaw, "Water unit weight", 0.0, cond=:(gammaw>0)),
]
@doc docstring(HydromechAnalysis_params) HydromechAnalysis()

mutable struct HydromechAnalysisProps<:TransientAnalysis
    stressmodel::Symbol # plane stress, plane strain, etc.
    thickness::Float64  # thickness for 2d analyses
    g::Float64 # gravity acceleration
    γw::Float64 # water unit weight
    
    function HydromechAnalysisProps(; kwargs...)
        args = checkargs(kwargs, HydromechAnalysis_params)
        this = new(args.stressmodel, args.thickness, args.g, args.gammaw)
        return this
    end
end

HydroAnalysis = HydromechAnalysis = HydromechAnalysisProps


# Assemble the global stiffness matrix
function hm_mount_global_matrices(model::FEModel,
                                  ndofs::Int,
                                  Δt::Float64,
                                  )
    # Assembling matrix G

    α = 1.0 # time integration factor

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]
        RHS = zeros(ndofs)

        for elem in model.elems
            ty                         = typeof(elem)
            has_stiffness_matrix       = hasmethod(elem_stiffness, (ty,))
            has_coupling_matrix        = hasmethod(elem_coupling_matrix, (ty,))
            has_conductivity_matrix    = hasmethod(elem_conductivity_matrix, (ty,))
            has_compressibility_matrix = hasmethod(elem_compressibility_matrix, (ty,))
            has_RHS_vector             = hasmethod(elem_RHS_vector, (ty,))

            # Assemble the stiffness matrix
            if has_stiffness_matrix
                K, rmap, cmap = elem_stiffness(elem)
                nr, nc = size(K)
                for i in 1:nr
                    for j in 1:nc
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, K[i,j])
                    end
                end
            end

            # Assemble the coupling matrices
            if has_coupling_matrix
                Cup, rmap, cmap = elem_coupling_matrix(elem)
                nr, nc = size(Cup)
                for i in 1:nr
                    for j in 1:nc
                        # matrix Cup
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, Cup[i,j])

                        # matrix Cup'
                        push!(R, cmap[j])
                        push!(C, rmap[i])
                        push!(V, Cup[i,j])
                    end
                end
            end

            # Assemble the conductivity matrix
            if has_conductivity_matrix
                H, rmap, cmap =  elem_conductivity_matrix(elem)
                nr, nc = size(H)
                for i in 1:nr
                    for j in 1:nc
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, α*Δt*H[i,j])
                    end
                end

                # Assembling RHS components
                Uw = [ dof.vals[:uw] for node in elem.nodes for dof in node.dofs if dof.name==:uw ]
                RHS[rmap] -= Δt*(H*Uw)
            end

            # Assemble the compressibility matrix
            if has_compressibility_matrix
                Cpp, rmap, cmap =  elem_compressibility_matrix(elem)
                nr, nc = size(Cpp)
                for i in 1:nr
                    for j in 1:nc
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, Cpp[i,j])
                    end
                end
            end

            # Assemble remaining RHS vectors
            if has_RHS_vector
                Q, map = elem_RHS_vector(elem)
                RHS[map] += Δt*Q
            end
        end
    end

    # generating sparse matrix G
    local G
    try
        G = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show ndofs
        @show err
    end

    yield()

    return G, RHS
end


function complete_uw_h(model::FEModel)
    haskey(model.node_data, "uw") || return
    Uw = model.node_data["uw"]
    H  = model.node_data["h"]

    for elem in model.elems
        elem.shape.family==BULKCELL || continue
        elem.shape==elem.shape.basic_shape && continue
        npoints  = elem.shape.npoints
        nbpoints = elem.shape.basic_shape.npoints
        map = [ elem.nodes[i].id for i in 1:nbpoints ]
        Ue = Uw[map]
        He = H[map]
        C = elem.shape.nat_coords
        for i in nbpoints+1:npoints
            id = elem.nodes[i].id
            R = C[i,:]
            N = elem.shape.basic_shape.func(R)
            Uw[id] = dot(N,Ue)
            H[id] = dot(N,He)
        end
    end
end


function solve!(model::FEModel, ana::HydromechAnalysis; args...)
    name = "Solver for seepage and hydromechanical analyses"
    status = stage_iterator!(name, hm_stage_solver!, model; args...)
    return status
end

hm_stage_solver_params = [
    FunInfo(:mech_stage_solver!, "Solves a load stage of a hydromechanical analysis."),
    ArgInfo(:model, "FEModel object"),
    ArgInfo(:stage, "Stage object"),
    KwArgInfo(:tol, "Force tolerance", 0.01, cond=:(tol>0)),
    KwArgInfo(:dTmin, "Relative minimum increment size", 1e-7, cond=:(0<dTmin<1) ),
    KwArgInfo(:dTmax, "Relative maximum increment size", 0.1, cond=:(0<dTmax<1) ),
    KwArgInfo(:rspan, "Relative span to residue reapplication", 0.01, cond=:(0<rspan<1) ),
    KwArgInfo(:scheme, "Global solving scheme", :FE, values=(:FE, :ME, :BE, :Ralston) ),
    KwArgInfo(:maxits, "Maximum number of NR iterations", 5, cond=:(1<=maxits<=10)),
    KwArgInfo(:autoinc, "Flag to set auto-increments", false),
    KwArgInfo(:quiet, "Flat to set silent mode", false),
]
@doc docstring(hm_stage_solver_params) hm_stage_solver!()

function hm_stage_solver!(model::FEModel, stage::Stage; args...)
    args = checkargs(args, mech_stage_solver_params)
    
    tol     = args.tol      
    ΔTmin   = args.dTmin    
    ΔTmax   = args.dTmax   
    rspan   = args.rspan    
    scheme  = args.scheme   
    maxits  = args.maxits 
    autoinc = args.autoinc  
    quiet   = args.quiet 

    ctx = model.ctx
    println(ctx.log, "Hydromech FE analysis: Stage $(stage.id)")

    solstatus = success()

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    tspan     = stage.tspan
    ctx       = model.ctx
    saveouts = stage.nouts > 0

    stressmodel = ctx.stressmodel
    ctx.ndim==3 && @check stressmodel==:d3

    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, stage.bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown bcs
    pmap  = nu+1:ndofs   # map for prescribed bcs
    model.ndofs = length(dofs)

    println(ctx.info, "unknown dofs: $nu")
    println(ctx.log, "unknown dofs: $nu")

    quiet || nu==ndofs && println(ctx.alerts, "solve_system!: No essential boundary conditions")

    if stage.id == 1
        # Setup quantities at dofs
        for dof in dofs
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
            if dof.name==:uw
                dof.vals[:h] = 0.0 # water head
            end
        end

        # Save initial file and loggers
        update_records!(model, force=true)
        complete_uw_h(model)
    end

    # get elevation Z for all Dofs
    Z = zeros(ndofs)
    for node in model.nodes
        for dof in node.dofs
            Z[dof.eq_id] = node.coord[ctx.ndim]
        end
    end

    # Get global parameters
    gammaw = model.ctx.γw
    isnan(gammaw) && error("solve!: gammaw parameter was not set in FEModel")
    gammaw > 0 || error("hm_solve: invalid value for gammaw: $gammaw")

    # Get the domain current state and backup
    State = [ ip.state for elem in model.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    t    = ctx.t     # current time
    ΔTbk = 0.0
    ΔTcheck = saveouts ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT,0.01))

    inc  = 0             # increment counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration
    Rc   = zeros(ndofs)  # vector of cumulated residues
    sysstatus = ReturnStatus()

    # Get boundary conditions
    Uex, Fex = get_bc_vals(model, bcs, t) # get values at time t  #TODO pick internal forces and displacements instead!

    # Get unbalanced forces
    Fin = zeros(ndofs)
    for elem in model.elems
        elem_internal_forces(elem, Fin)
    end
    Fex .-= Fin # add negative forces to external forces vector

    # Get global vectors from values at dofs
    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}

    while T < 1.0-ΔTmin
        ctx.ΔT = ΔT

        # Update counters
        inc += 1
        ctx.inc = inc

        println(ctx.log, "  inc $inc")

        # Get forces and displacements from boundary conditions
        Δt = tspan*ΔT
        ctx.t = t + Δt
        UexN, FexN = get_bc_vals(model, bcs, t+Δt) # get values at time t+Δt

        ΔUex = UexN - U
        ΔFex = FexN - F

        ΔTcr = min(rspan, 1-T)    # time span to apply cumulated residues
        αcr  = min(ΔT/ΔTcr, 1.0)  # fraction of cumulated residues to apply
        T<1-rspan && (ΔFex .+= αcr.*Rc) # addition of residuals

        ΔUex[umap] .= 0.0
        ΔFex[pmap] .= 0.0

        R   .= ΔFex
        ΔUa .= 0.0
        ΔUi .= ΔUex    # essential values at iteration i

        # Newton Rapshon iterations
        residue   = 0.0
        maxfails  = 3  # maximum number of it. fails with residual change less than 90%
        nfails    = 0  # counter for iteration fails
        nits      = 0
        residue1  = 0.0
        converged = false
        syserror  = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Predictor step for FE, ME and BE
            if scheme in (:FE, :ME, :BE)
                G, RHS = hm_mount_global_matrices(model, ndofs, Δt)
                R .+= RHS

                # Solve
                status = solve_system!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                println(ctx.log, status.message)


                failed(status) && (syserror=true; break)

                copyto!.(State, StateBk)
                ΔUt    = ΔUa + ΔUi
                ΔFin, status = update_state!(model.elems, ΔUt, Δt)
                failed(status) && (syserror=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap] )
            end

            # Corrector step for ME and BE
            if residue > tol && scheme in (:ME, :BE)
                G2, RHS = hm_mount_global_matrices(model, ndofs, Δt)
                if scheme=="ME"
                    G = 0.5*(G + G2)
                elseif scheme=="BE"
                    G = G2
                end
                status = solve_system!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(status) && (syserror=true; break)

                copyto!.(State, StateBk)
                ΔUt    = ΔUa + ΔUi
                ΔFin, status = update_state!(model.elems, ΔUt, Δt)
                failed(status) && (syserror=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions

            @printf(ctx.log, "    it %d  residue: %-10.4e\n", it, residue)

            it==1 && (residue1=residue)
            residue < tol  && (converged=true; break)
            isnan(residue) && break
            it>maxits      && break
            it>1 && res>lastres && break

        end

        ctx.residue = residue
        q = 0.0 # increment size factor for autoinc

        if syserror
            println(ctx.alerts, sysstatus.message)
            println(ctx.log, sysstatus.message)
            converged = false
        end
        # quiet || sysstatus.message!="" && println(ctx.alerts, sysstatus.message)

        if converged
            # Update nodal natural and essential values for the current stage
            U .+= ΔUa
            F .+= ΔFin
            Uex .= UexN
            Fex .= FexN
            Rc .= (1.0-αcr).*Rc .+ R  # update cumulated residue

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
                if dof.name==:uw
                    dof.vals[:h] = Z[i] + U[i]/gammaw
                end
            end

            # Update time
            t += Δt
            T += ΔT
            ctx.t = t
            ctx.T = T

            # Check for saving output file
            checkpoint = T>Tcheck-ΔTmin
            if checkpoint
                ctx.out += 1
                Tcheck += ΔTcheck # find the next output time
            end

            rstatus = update_records!(model, checkpoint=checkpoint)
            if failed(rstatus)
                println(ctx.alerts, rstatus.message)
                return rstatus
            end

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    q = 1+tanh(log10(ftol/(res1+eps())))
                    q = max(q, 1.1)

                    ΔTtr = min(q*ΔT, ΔTmax, 1-T)
                    if T+ΔTtr>Tcheck-ΔTmin
                        ΔTbk = ΔT
                        ΔT = Tcheck-T
                        @assert ΔT>=0.0
                    else
                        ΔT = ΔTtr
                        ΔTbk = 0.0
                    end
                end
            end
        else
            # Restore counters
            inc -= 1
            ctx.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(ctx.log, "      increment failed")
                q = (1+tanh(log10(tol/(residue1+eps()))))
                q = clamp(q, 0.2, 0.9)
                syserror && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < ΔTmin
                    solstatus = failure("Solver did not converge.")
                    break
                end
            else
                solstatus = failure("Residue is greater than the tolerance. Try `autoinc=true`. ")
                break
            end
        end
    end

    complete_uw_h(model)
    failed(solstatus) && update_records!(model, force=true)

    return solstatus

end