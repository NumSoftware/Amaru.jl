# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ThermoAnalysis, ThermomechAnalysis

ThermomechAnalysis_params = [
    FunInfo(:ThermomechAnalysis, "Thermomechanical analysis properties."),
    KwArgInfo(:stressmodel, "Stress model", :d3, values=(:planestress, :planestrain, :axisymmetric, :d3)),
    KwArgInfo(:thickness, "Thickness for 2d analyses", 1.0, cond=:(thickness>0)),
    KwArgInfo(:g, "Gravity acceleration", 0.0, cond=:(g>=0)),
    KwArgInfo(:T0, "Reference temperature", 0.0, cond=:(T0>=-273.15)),
]
@doc docstring(ThermomechAnalysis_params) ThermomechAnalysis()

mutable struct ThermomechAnalysisProps<:TransientAnalysis
    stressmodel::Symbol # plane stress, plane strain, etc.
    thickness::Float64  # thickness for 2d analyses
    g::Float64 # gravity acceleration
    T0::Float64 # reference temperature
    
    function ThermomechAnalysisProps(; kwargs...)
        args = checkargs(kwargs, ThermomechAnalysis_params)
        this = new(args.stressmodel, args.thickness, args.g, args.T0)
        return this
    end
end

ThermoAnalysis = ThermomechAnalysis = ThermomechAnalysisProps

# Assemble the global stiffness matrix
function tm_mount_global_matrices(model::Model,
                                  ndofs::Int,
                                  Δt::Float64,
                                 )
    # Assembling matrix G

    α = 1.0 # time integration factor
    T0k = model.env.ana.T0 + 273.15

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]
        RHS = zeros(ndofs)

        for elem in model.elems
            ty                      = typeof(elem)
            has_stiffness_matrix    = hasmethod(elem_stiffness, (ty,))
            has_coupling_matrix     = hasmethod(elem_coupling_matrix, (ty,))
            has_conductivity_matrix = hasmethod(elem_conductivity_matrix, (ty,))
            has_mass_matrix         = hasmethod(elem_mass_matrix, (ty,))

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
                Cut, rmap, cmap = elem_coupling_matrix(elem)
                nr, nc = size(Cut)
                for i in 1:nr
                    for j in 1:nc
                        # matrix Cut
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, Cut[i,j])

                        # matrix Cut'
                        push!(R, cmap[j])
                        push!(C, rmap[i])
                        push!(V, T0k*Cut[i,j]) # transposed and multiplied by T0k
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
                Ut = [ dof.vals[:ut] for node in elem.nodes for dof in node.dofs if dof.name==:ut ]
                RHS[rmap] -= Δt*(H*Ut)
            end

            # Assemble the mass matrix
            if has_mass_matrix
                M, rmap, cmap =  elem_mass_matrix(elem)
                nr, nc = size(M)
                for i in 1:nr
                    for j in 1:nc
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, M[i,j])
                    end
                end
            end
        end
    end

    # generating sparse matrix G
    local G
    try
        G = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show err
    end

    yield()

    return G, RHS
end


function complete_ut_T(model::Model)
    haskey(model.node_data, "ut") || return
    Ut = model.node_data["ut"]
    T0 = model.env.ana.T0

    for elem in model.elems
        elem.shape.family==BULKCELL || continue
        elem.shape==elem.shape.basic_shape && continue
        npoints  = elem.shape.npoints
        nbpoints = elem.shape.basic_shape.npoints
        map = [ elem.nodes[i].id for i in 1:nbpoints ]
        Ute = Ut[map]
        C = elem.shape.nat_coords
        for i in nbpoints+1:npoints
            id = elem.nodes[i].id
            R = C[i,:]
            N = elem.shape.basic_shape.func(R)
            Ut[id] = dot(N,Ute)
        end
    end

    model.node_data["ut"] = Ut .+ T0
end


function solve!(model::Model, ana::ThermomechAnalysis; args...)
    name = "Solver for thermal and thermomechanical analyses"
    status = stage_iterator!(name, tm_stage_solver!, model; args...)
    return status
end

tm_stage_solver_params = [
    FunInfo(:tm_stage_solver!, "Solves a load stage of a thremomechanical analysis."),
    ArgInfo(:model, "Model object"),
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
@doc docstring(tm_stage_solver_params) tm_stage_solver!()

function tm_stage_solver!(model::Model, stage::Stage; args...)
    args = checkargs(args, tm_stage_solver_params)
    
    tol     = args.tol      
    ΔTmin   = args.dTmin    
    ΔTmax   = args.dTmax   
    rspan   = args.rspan    
    scheme  = args.scheme   
    maxits  = args.maxits 
    autoinc = args.autoinc  
    quiet   = args.quiet 

    env = model.env
    println(env.log, "Thermomech FE analysis: Stage $(stage.id)")

    solstatus = success()

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    tspan     = stage.tspan
    env       = model.env
    saveouts = stage.nouts > 0
    T0        = env.ana.T0
    ftol      = tol

    stressmodel = env.ana.stressmodel
    env.ndim==3 && @check stressmodel==:d3
    

    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, stage.bcs) # unknown dofs first
    ndofs    = length(dofs)
    umap     = 1:nu         # map for unknown bcs
    pmap     = nu+1:ndofs   # map for prescribed bcs
    model.ndofs = length(dofs)

    println(env.info,"unknown dofs: $nu")
    println(env.log, "unknown dofs: $nu")

    quiet || nu==ndofs && println(env.alerts, "No essential boundary conditions")

    if stage.id == 1
        # Setup quantities at dofs
        for dof in dofs
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
            if dof.name==:ut
                dof.vals[:T] = T0 # real temperature
            end
        end

        update_records!(model, force=true)
        complete_ut_T(model)
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in active_elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    t       = env.t     # current time
    ΔTbk    = 0.0
    ΔTcheck = saveouts ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT, ΔTmax, ΔTcheck))

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
    Uex, Fex = get_bc_vals(model, bcs, t) # get values at time t  #TODO check for incremental or target values

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

    if scheme==:FE
        p1=1.0; q11=1.0
    elseif scheme==:ME
        p1=1.0; q11=1.0; a1=0.5; a2=0.5
    elseif scheme==:BE
        p1=1.0; q11=1.0; a1=0.0; a2=1.0
    elseif scheme==:Ralston
        p1=2/3; q11=2/3; a1=1/4; a2=3/4
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}

    while T < 1.0-ΔTmin
        env.ΔT = ΔT

        # Update counters
        inc += 1
        env.inc = inc

        println(env.log, "  inc $inc")

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
        ΔUi .= ΔUex  # essential values at iteration i

        # Newton Rapshon iterations
        nits      = 0
        res       = 0.0
        res1      = 0.0
        converged = false
        syserror  = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = res # residue from last iteration

            # Predictor step
            G, RHS = tm_mount_global_matrices(model, ndofs, Δt)
            ΔUitr = p1*ΔUi
            Rtr   = q11*(R + RHS)
            sysstatus = solve_system!(G, ΔUitr, Rtr, nu)   # Changes unknown positions in ΔUi and R
            failed(sysstatus) && (syserror=true; break)
            copyto!.(State, StateBk)
            ΔUt = ΔUa + ΔUitr
            ΔFin, sysstatus = update_state!(active_elems, ΔUt, Δt)
            failed(sysstatus) && (syserror=true; break)

            # Corrector step
            if scheme==:FE
                ΔUi = ΔUitr
            else
                G2, RHS = tm_mount_global_matrices(model, ndofs, Δt)
                R .+= RHS
                G = a1*G + a2*G2
                R = R + RHS
                sysstatus = solve_system!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (syserror=true; break)
                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                ΔFin, sysstatus = update_state!(active_elems, ΔUt, Δt)
                failed(sysstatus) && (syserror=true; break)
            end

            # Update accumulated essential values
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions
            res = norm(R, Inf)

            @printf(env.log, "    it %d  residue: %-10.4e\n", it, res)

            it==1 && (res1=res)
            res<ftol && (converged=true; break)

            isnan(res) && break
            it>maxits  && break

            it>1 && res>lastres && break
        end

        env.residue = res
        q = 0.0 # increment size factor for autoinc

        if syserror
            println(env.alerts, sysstatus.message)
            println(env.log, sysstatus.message)
            converged = false
        end

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
                if dof.name==:ut
                    dof.vals[:T] = U[i] + T0
                end
            end

            # Update time
            t += Δt
            T += ΔT
            env.t = t
            env.T = T

            # Check for saving output file
            checkpoint = T>Tcheck-ΔTmin
            if checkpoint
                env.out += 1
                Tcheck += ΔTcheck # find the next output time
            end

            complete_ut_T(model)
            rstatus = update_records!(model, checkpoint=checkpoint)
            if failed(rstatus)
                println(env.alerts, rstatus.message)
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
            env.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(env.log, "      increment failed")

                q = 1+tanh(log10(ftol/(res1+eps())))
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

    complete_ut_T(model)
    failed(solstatus) && update_records!(model, force=true)

    return solstatus

end
