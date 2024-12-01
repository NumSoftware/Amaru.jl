# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export AcousticModalAnalysis


AcousticModalAnalysis_params = [
    FunInfo(:DynAnalysis, "Dynamic analysis"),
    ArgInfo(:model, "Finite element model"),
]
@doc docstring(AcousticModalAnalysis_params) AcousticModalAnalysis


mutable struct AcousticModalAnalysis<:Analysis
    model ::FEModel
    ctx   ::AcousticMechContext
    sctx  ::SolverContext

    stages  ::Array{Stage}
    loggers ::Array{AbstractLogger,1}
    monitors::Array{AbstractMonitor,1}

    freqs::Array{Float64,1} # frequencies
    modes::Array{Float64,2} # modes
    
    function AcousticModalAnalysis(model::FEModel; outdir=".", outkey="out")
        this = new(model, model.ctx)
        this.stages = []
        this.loggers = []
        this.monitors = []  

        this.freqs = zeros(0)
        this.modes = zeros(0,0)

        this.sctx = SolverContext()
        this.sctx.outdir = outdir
        this.sctx.outkey = outkey

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


function solve!(model::FEModel, ana::AcousticModalAnalysis; args...)
    name = "Solver for dynamic modal analyses"
    status = stage_iterator!(name, acoustic_modal_solver!, model; args...)
    return status
end


acoustic_modal_solver_params = [
    FunInfo( :acoustic_modal_solver!, "Finds the frequencies and vibration modes of an acoustic system."),
    ArgInfo( :model, "FEModel object"),
    ArgInfo( :stage, "Stage object"),
    KwArgInfo( :nmodes, "Number of modes to be computed", 5),
    KwArgInfo( :quiet, "Flag to set silent mode", false),
]
@doc docstring(acoustic_modal_solver_params) solve!(::AcousticModalAnalysis; args...)


function solve!(ana::AcousticModalAnalysis; args...)
    args = checkargs(args, acoustic_modal_solver_params)
    if !args.quiet
        printstyled("Solver for acousticanical modal analyses", "\n", bold=true, color=:cyan)
        println("  stress model: ", ana.ctx.stressmodel)
    end

    status = stage_iterator!(acoustic_modal_solver!, ana; args...)
    return status
end


function acoustic_modal_solver!(ana::AcousticModalAnalysis, stage::Stage; kwargs...)
    args   = checkargs(kwargs, acoustic_modal_solver_params)
    nmodes = args.nmodes
    quiet  = args.quiet

    model = ana.model
    ctx = model.ctx
    sctx = ana.sctx

    quiet || println(sctx.log, "Modal analysis for acoustic systems")

    stressmodel = ctx.stressmodel
    ctx.ndim==3 && @check stressmodel==:d3

    
    # todo: check there are not force boundary conditions

    # get only bulk elements
    # model = FEModel(model.elems.bulks) #todo: update surface/edges

    # check density
    for elem in model.elems
        elem.props.ρ==0 && error("mech_modal_solver: density should not be zero")
    end

    # get dofs organized according to boundary conditions
    dofs, nu    = configure_dofs!(model, stage.bcs)

    ndofs       = length(dofs)
    model.ndofs = length(dofs)
    if !quiet
        println(sctx.log, "unknown dofs: $nu")
        println(sctx.info, "unknown dofs: $nu")
    end

    # setup quantities at dofs
    for dof in dofs
        dof.vals[dof.name]    = 0.0
        dof.vals[dof.natname] = 0.0
    end

    K11 = am_mount_K(model.elems, ndofs)[1:nu, 1:nu]
    M11 = am_mount_M(model.elems, ndofs)[1:nu, 1:nu]

    L = zeros(nu) # vector for the inverse of the lumped mass matrix
    for i in 1:size(M11,2)
        L[i] = 1.0/M11[i,i]
    end

    L = diag(M11) # lumped matrix in vector form
    P = K11./L  # inv(M11)*K11  mass normalized stiffness matrix

    # invM = inv(Matrix(M11))
    # P = invM*K11

    # eingenvalues and eingenvectors
    λ, V = eigs(P, nev=nmodes, which=:SR) # find λ with the smallest real parts / :SM for smallest magnitude

    # filter positive real eigenvals
    filter = [ i for i in eachindex(λ) if isreal(λ[i]) && real(λ[i])>0 ]
    perm   = sortperm(real(λ[filter]))[1:nmodes] # sorting from smalles to largest
    filter = filter[perm]

    λ = λ[filter]
    V = V[:, filter]

    ω = λ.^0.5  # true frequencies

    update_output_data!(model)
    save(model, joinpath(sctx.outdir, "$(sctx.outkey)-0.vtu"), quiet=true)

    sctx.inc = 1
    sctx.ΔT  = 1.0
    
    # save modes
    for i in 1:nmodes
        U = V[:,i] # modal displacements
        
        for (k,dof) in enumerate(dofs[1:nu])
            dof.vals[dof.name] = U[k]
        end
        sctx.T = i/nmodes
        sctx.out = i
        update_output_data!(model)
        save(model, joinpath(sctx.outdir, "$(sctx.outkey)-$i.vtu"), quiet=true)
    end

    # reset displacement values
    for dof in dofs
        dof.vals[dof.name] = 0.0
    end

    # save frequencies and modes
    ana.freqs = ω
    ana.modes = V

    return success()
end

