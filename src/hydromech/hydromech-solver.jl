# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export HydroAnalysis, HydromechAnalysis

mutable struct HydromechAnalysisProps<:Analysis
    stressmodel::String # plane stress, plane strain, etc.
    thickness::Float64  # thickness for 2d analyses
    g::Float64 # gravity acceleration
    γw::Float64 # water unit weight
    
    function HydromechAnalysisProps(;stressmodel="3d", thickness=1.0, g=0.0, gammaw=0)
        @check stressmodel in ("plane-stress", "plane-strain", "axisymmetric", "3d")
        @check thickness>0
        @check g>=0
        @check gammaw>0
        return new(stressmodel, thickness, g, gammaw)
    end
end

HydroAnalysis = HydromechAnalysis = HydromechAnalysisProps


# Assemble the global stiffness matrix
function hm_mount_global_matrices(model::Model,
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


function complete_uw_h(model::Model)
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


# function update_state!(model::Model, ΔUt::Vect, ΔFin::Vect, Δt::Float64, verbosity::Int)
#     # Update
#     verbosity>1 && print("    updating... \r")

#     # Get internal forces and update data at integration points (update ΔFin)
#     ΔFin .= 0.0
#     for elem in model.elems
#         status = update_elem!(elem, ΔUt, ΔFin, Δt)
#         failed(status) && return status
#     end
#     return success()
# end


"""
    solve!(domain, bcs, options...) -> Bool

Performs one stage finite element hydro-mechanical analysis of a `domain`
subjected to a list of boundary conditions `bcs`.

# Arguments

`model` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)

# Keyword arguments

`tspan = NaN` : Time lapse for the transient analysis in the current stage

`end_time = NaN` : Final time for the transient analysis in the current stage

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxincs = 1000000` : Maximum number of increments

`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`Ttol     = 1e-9` : Pseudo-time tolerance

`scheme  = :FE` : Predictor-corrector scheme at each increment. Available scheme is `:FE`

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`verbose = true` : If true, provides information of the analysis steps

`silent = false` : If true, no information is printed
"""

function solve!(model::Model, ana::HydromechAnalysis; args...)
    name = "Solver for seepage and hydromechanical analyses"
    status = stage_iterator!(name, hm_stage_solver!, model; args...)
    return status
end


function hm_stage_solver!(model::Model, stage::Stage, logfile::IOStream, sline::StatusLine; 
    tol     :: Number  = 1e-2,
    Ttol    :: Number  = 1e-8,
    rspan   :: Number  = 1e-2,
    scheme  :: String  = "FE",
    maxits  :: Int     = 5,
    autoinc :: Bool    = false,
    outdir  :: String  = ".",
    outkey  :: String  = "out",
    quiet  :: Bool    = false
    )

    println(logfile, "Hydromech FE analysis: Stage $(stage.id)")
    stage.status = :solving

    solstatus = success()
    scheme in ("FE", "ME", "BE") || error("solve! : invalid scheme \"$(scheme)\"")

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    tspan     = stage.tspan
    env       = model.env
    save_outs = stage.nouts > 0

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
    println(logfile, "unknown dofs: $nu")
    message(sline, "  unknown dofs: $nu")

    quiet || nu==ndofs && message(sline, "solve_system!: No essential boundary conditions", Base.warn_color)

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
        update_output_data!(model)
        update_single_loggers!(model)
        update_multiloggers!(model)
        update_monitors!(model)
        complete_uw_h(model)
        save_outs && save(model, "$outdir/$outkey-0.vtu", quiet=true)
    end

    # get elevation Z for all Dofs
    Z = zeros(ndofs)
    for node in model.nodes
        for dof in node.dofs
            Z[dof.eq_id] = node.coord[env.ndim]
        end
    end

    # Get global parameters
    gammaw = model.env.ana.γw
    isnan(gammaw) && error("solve!: gammaw parameter was not set in Model")
    gammaw > 0 || error("hm_solve: invalid value for gammaw: $gammaw")

    # Get the domain current state and backup
    State = [ ip.state for elem in model.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    t    = env.t     # current time

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT,0.01))

    ΔTbk = 0.0
    ΔTcheck = save_outs ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    inc  = 0             # increment counter
    iout = env.out # file output counter
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

    while T < 1.0-Ttol
        env.ΔT = ΔT

        # Update counters
        inc += 1
        env.inc = inc

        println(logfile, "  inc $inc")

        # Get forces and displacements from boundary conditions
        Δt = tspan*ΔT
        env.t = t + Δt
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
        errored   = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Predictor step for FE, ME and BE
            if scheme in ("FE", "ME", "BE")
                G, RHS = hm_mount_global_matrices(model, ndofs, Δt)
                R .+= RHS

                # Solve
                status = solve_system!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                println(logfile, status.message)


                failed(status) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt    = ΔUa + ΔUi
                ΔFin, status = update_state!(model.elems, ΔUt, Δt)
                failed(status) && (errored=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap] )
            end

            # Corrector step for ME and BE
            if residue > tol && scheme in ("ME", "BE")
                G2, RHS = hm_mount_global_matrices(model, ndofs, Δt)
                if scheme=="ME"
                    G = 0.5*(G + G2)
                elseif scheme=="BE"
                    G = G2
                end
                status = solve_system!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(status) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt    = ΔUa + ΔUi
                ΔFin, status = update_state!(model.elems, ΔUt, Δt)
                failed(status) && (errored=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions

            @printf(logfile, "    it %d  residue: %-10.4e\n", it, residue)

            it==1 && (residue1=residue)
            residue < tol  && (converged=true; break)
            isnan(residue) && break
            it>maxits      && break
            #it>1 && residue>lastres && break
            residue>0.9*lastres && (nfails+=1)
            nfails==maxfails    && break
        end

        q = 0.0 # increment size factor for autoinc

        if errored
            println(logfile, sysstatus.message)
            converged = false
        end
        quiet || sysstatus.message!="" && message(sline, sysstatus.message, Base.default_color_warn)

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
            env.t = t
            env.T = T
            env.residue = residue

            # Check for saving output file
            if T>Tcheck-Ttol && save_outs
                env.out += 1
                iout = env.out
                rm.(glob("*conflicted*.dat", "$outdir/"), force=true)
                
                update_output_data!(model)
                # update_embedded_disps!(active_elems, model.node_data["U"])

                update_multiloggers!(model)
                save(model, "$outdir/$outkey-$iout.vtu", quiet=true)

                Tcheck += ΔTcheck # find the next output time
            end

            update_single_loggers!(model)
            update_monitors!(model)
            flush(logfile)

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    if nits==1
                        q = (1+tanh(log10(tol/residue1)))^1
                    else
                        q = 1.0
                    end

                    ΔTtr = min(q*ΔT, 1/nincs, 1-T)
                    if T+ΔTtr>Tcheck-Ttol
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
            env.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(logfile, "      increment failed")
                q = (1+tanh(log10(tol/residue1)))
                q = clamp(q, 0.2, 0.9)
                errored && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < Ttol
                    solstatus = failure("solver did not converge")
                    break
                end
            else
                solstatus = failure("solver did not converge")
                break
            end
        end
    end

    if !save_outs
        update_output_data!(model)
        complete_uw_h(model)
        update_multiloggers!(model)
    end

    return solstatus

end