# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export AcousticAnalysis, AcousticMechAnalysis, AcousticContext, AcousticMechContext


AcousticMechContext_params = [
    FunInfo(:MechContext, "Context for mechanical analysis"),
    KwArgInfo(:ndim, "Analysis dimension", 0),
    KwArgInfo(:stressmodel, "Stress model", :d3, values=(:planestress, :planestrain, :axisymmetric, :d3, :none)),
    KwArgInfo(:g, "Gravity acceleration", 0.0, cond=:(g>=0))
]
@doc docstring(AcousticMechContext_params) AcousticMechContext

mutable struct AcousticMechContext<:Context
    stressmodel::Symbol # plane stress, plane strain, etc.
    g::Float64 # gravity acceleration

    ndim     ::Int       # Analysis dimension
    thickness::Float64 # to be set after FEModel creation

    function AcousticMechContext(; kwargs...)
        args = checkargs(kwargs, AcousticMechContext_params)
        this = new()
        
        # Analysis related
        this.stressmodel = args.stressmodel
        this.ndim        = args.ndim
        this.g           = args.g
        return this
    end
end

const AcousticContext = AcousticMechContext


AcousticMechAnalysis_params = [
    FunInfo(:AcousticMechAnalysis, "AcousticMechanical analysis."),
]
@doc docstring(AcousticMechAnalysis_params) AcousticMechAnalysis()

mutable struct AcousticMechAnalysis<:TransientAnalysis
    model ::FEModel
    ctx   ::AcousticMechContext
    sctx  ::SolverContext

    stages  ::Array{Stage}
    loggers ::Array{AbstractLogger,1}
    monitors::Array{AbstractMonitor,1}

    function AcousticMechAnalysis(model::FEModel; outdir=".", outkey="out")
        this = new(model, model.ctx)
        this.stages = []
        this.loggers = []
        this.monitors = []  
        this.sctx = SolverContext()
        
        this.sctx.outkey = outkey
        this.sctx.outdir = rstrip(outdir, ['/', '\\'])
        isdir(this.sctx.outdir) || mkdir(this.sctx.outdir) # create output directory if it does not exist

        model.ctx.thickness = model.thickness
        if model.ctx.stressmodel==:none
            if model.ctx.ndim==2
                model.ctx.stressmodel = :planestrain
            else
                model.ctx.stressmodel = :d3
            end
        end

        return this
    end
    
end

const AcousticAnalysis = AcousticMechAnalysis

function am_mount_M(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_mass, (ty,)) || continue
            
            Me, rmap, cmap = elem_acoustic_mass(elem)
            
            nr, nc = size(Me)
            for i in 1:nr
                for j in 1:nc
                    val = Me[i,j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    M = sparse(R, C, V, ndofs, ndofs)
    yield()

    return M
end


function am_mount_K(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_stiffness, (ty,)) || continue

            Ke, rmap, cmap = elem_acoustic_stiffness(elem)

            nr, nc = size(Ke)
            for i in 1:nr
                for j in 1:nc
                    val = Ke[i,j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    K = sparse(R, C, V, ndofs, ndofs)
    yield()

    return K
end


function solve!(model::FEModel, ana::AcousticMechAnalysis; args...)
    name = "Solver for acoustic mechanical analyses"
    status = stage_iterator!(name, am_stage_solver!, model; args...)
    return status
end

am_stage_solver_params = [
    FunInfo(:tm_stage_solver!, "Solves a load stage of a thremomechanical analysis."),
    ArgInfo(:model, "FEModel object"),
    KwArgInfo(:tol, "Force tolerance", 0.01, cond=:(tol>0)),
    KwArgInfo(:dTmin, "Relative minimum increment size", 1e-7, cond=:(0<dTmin<1) ),
    KwArgInfo(:dTmax, "Relative maximum increment size", 0.1, cond=:(0<dTmax<1) ),
    KwArgInfo(:rspan, "Relative span to residue reapplication", 0.01, cond=:(0<rspan<1) ),
    KwArgInfo(:scheme, "Global solving scheme", :FE, values=(:FE, :ME, :BE, :Ralston) ),
    KwArgInfo(:maxits, "Maximum number of NR iterations", 5, cond=:(1<=maxits<=10)),
    KwArgInfo(:autoinc, "Flag to set auto-increments", false),
    KwArgInfo(:quiet, "Flat to set silent mode", false),
]
@doc docstring(tm_solver_params) solve!(::AcousticMechAnalysis; args...)


function solve!(ana::AcousticMechAnalysis; args...)
    args = checkargs(args, tm_solver_params)
    if !args.quiet
        printstyled("Solver for acoustic-mechanical analyses", "\n", bold=true, color=:cyan)
        println("  stress model: ", ana.ctx.stressmodel)
    end

    status = stage_iterator!(tm_stage_solver!, ana; args...)
    return status
end


function tm_stage_solver!(ana::AcousticMechAnalysis, stage::Stage; args...)
    args = NamedTuple(args)
    
    tol     = args.tol      
    ΔTmin   = args.dTmin    
    ΔTmax   = args.dTmax   
    # rspan   = args.rspan    
    # scheme  = args.scheme   
    maxits  = args.maxits 
    autoinc = args.autoinc  
    quiet   = args.quiet 

    model = ana.model
    ctx = model.ctx
    sctx = ana.sctx
    println(sctx.log, "ThermoMech FE analysis: Stage $(stage.id)")

    solstatus = success()

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    tspan     = stage.tspan
    saveouts  = stage.nouts > 0

    stressmodel = ctx.stressmodel
    ctx.ndim==3 && @check stressmodel==:d3

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

    println(sctx.info, "unknown dofs: $nu")
    println(sctx.log, "unknown dofs: $nu")

    # quiet || nu==ndofs && println(sctx.alerts, "solve_system!: No essential boundary conditions", Base.warn_color) #TODO

    # Dictionary of data keys related with a dof
    components_dict = Dict( 
                            :up => (:up, :fq, :vp, :ap),
                            :ux => (:ux, :fx, :vx, :ax),
                            :uy => (:uy, :fy, :vy, :ay),
                            :uz => (:uz, :fz, :vz, :az),
                            :rx => (:rx, :mx, :vrx, :arx),
                            :ry => (:ry, :my, :vry, :ary),
                            :rz => (:rz, :mz, :vrz, :arz))
    
    # Setup quantities at dofs
    if stage.id == 1
        for dof in dofs
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = 0.0
            dof.vals[vs] = 0.0
            dof.vals[as] = 0.0
        end

        # Initial accelerations
        A = zeros(ndofs)
        V = zeros(ndofs)
        Uex, Fexa = get_bc_vals(model, bcs) # get values at time t
        
        M = am_mount_M(model.elems, ndofs)

        # initial acceleration
        sysstatus = solve_system!(M, A, Fexa, nu)
        # failed(sysstatus) && return failure("solver: Reduced mass matrix is singular")


        # Initial values at nodes
        for (i,dof) in enumerate(dofs)
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[vs] = V[i]
            dof.vals[as] = A[i]
            dof.vals[fs] = Fexa[i]
        end

        # Save initial file and loggers
        update_records!(ana, force=true)
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in active_elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    t       = sctx.t     # current time
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
    ΔUi  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUk  = zeros(ndofs)  # vector of essential values for current iteration
    
    ΔFexi = zeros(ndofs)  # increment of external natural bc
    ΔUexi = zeros(ndofs)  # increment of external essential bc
    Fexi! = zeros(ndofs)  # last value of external natural bc
    Uexi! = zeros(ndofs)  # last value of external essential bc

    sysstatus = ReturnStatus()

    # Get boundary conditions
    # Uex, Fexa = get_bc_vals(model, bcs, t) # get values at time t  #TODO pick internal forces and displacements instead!

    # for elem in model.elems
    #     elem_internal_forces(elem, Fin)
    # end
    # Fexa .-= Fin # add negative forces to external forces vector

    # Get global vectors from values at dofs
    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}
    local At, Vt

    while T < 1.0-ΔTmin
        sctx.ΔT = ΔT
        Δt     = tspan*ΔT
        sctx.t  = t + Δt

        # Update counters
        inc += 1
        sctx.inc = inc 

        println(sctx.log, "  inc $inc")

        # Get forces and displacements from boundary conditions
        Uexi, Fexi = get_bc_vals(model, bcs, t+Δt) # get values at time t+Δt
        ΔUexi .= Uexi - Uexi!
        ΔF = Fexi - Fexi!

        ΔUexi[umap] .= 0.0
        ΔF[pmap] .= 0.0

        ΔF[pmap] .= 0.0  # Zero at prescribed positions

        ΔUi .= 0.0
        ΔUk .= ΔUexi    # essential values at iteration i

        # Newton Rapshon iterations
        res   = 0.0
        nits  = 0
        res1  = 0.0
        converged = false
        syserror   = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUk.=0.0) # essential values are applied only at first iteration
            lastres = res # res from last iteration

            K = am_mount_K(model.elems, ndofs)
            M = am_mount_M(model.elems, ndofs)

            K′ = K + 4/Δt^2*M   # pseudo-stiffness matrix
            ΔF′ = ΔF + M*(A + 4*V/Δt - 4*ΔUi/Δt^2)
            # Solve
            solve_system!(K′, ΔUk, ΔF′, nu)

            # Restore the state to last converged increment
            copyto!.(State, StateBk)

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin .= 0.0
            ΔUit   = ΔUi + ΔUk
            ΔFin, sysstatus = update_state!(model.elems, ΔUit, Δt)
            failed(sysstatus) && (syserror=true; break)

            Vt = -V + 2/Δt*ΔUit
            At = -A + 4/Δt^2*(ΔUit - V*Δt)

            R = ΔF - (ΔFin + M*At)
            res = maximum(abs, R[umap] )

            # Update accumulated displacement
            ΔUi .+= ΔUk

            @printf(sctx.log, "    it %d  res: %-10.4e\n", it, res)

            it==1 && (res1=res)
            res < tol  && (converged=true; break)
            isnan(res) && break
            it>maxits  && break
            it>1 && res>lastres && break
        end

        q = 0.0 # increment size factor for autoinc

        if syserror
            println(sctx.log, sysstatus.message)
            converged = false
        end
        quiet || sysstatus.message!="" && println(sctx.alerts, sysstatus.message)

        if converged
            # Update nodal natural and essential values for the current stage
            U  .+= ΔUi
            F  .+= ΔFin

            Uexi! = Uexi
            Fexi! .= Fexi

            A .= At
            V .= Vt

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                us, fs, vs, as = components_dict[dof.name]
                dof.vals[us] = U[i]
                dof.vals[fs] = F[i]
                dof.vals[vs] = V[i] 
                dof.vals[as] = A[i]
            end

            # Update time
            t += Δt
            T += ΔT
            sctx.t = t
            sctx.T = T
            sctx.residue = res

            # Check for saving output file
            checkpoint = T>Tcheck-ΔTmin
            if checkpoint
                sctx.out += 1
                Tcheck += ΔTcheck # find the next output time
            end

            rstatus = update_records!(ana, checkpoint=checkpoint)
            if failed(rstatus)
                println(sctx.alerts, rstatus.message)
                return rstatus
            end

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    if nits==1
                        q = 1.0 + tanh(log10(tol/res1))
                    else
                        q = 1.0
                    end

                    ΔTtr = min(q*ΔT, 1/nincs, 1-T)
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
            sctx.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(sctx.log, "      increment failed")
                q = (1+tanh(log10(tol/res1)))
                q = clamp(q, 0.2, 0.9)
                syserror && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < ΔTmin
                    solstatus = failure("solver did not converge")
                    break
                end
            else
                solstatus = failure("solver did not converge")
                break
            end
        end
    end

    failed(solstatus) && update_records!(ana, force=true)

    return solstatus

end
