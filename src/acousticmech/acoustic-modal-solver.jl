
export AcousticModalAnalysis

AcousticModalAnalysis_params = [
    FunInfo(:AcousticModalAnalysis, "Mechanical modal analysis"),
    KwArgInfo(:stressmodel, "Stress model", :d3, values=(:planestress, :planestrain, :axisymmetric, :d3)),
    KwArgInfo(:thickness, "Thickness for 2d analyses", 1.0, cond=:(thickness>0)),
    KwArgInfo(:g, "Gravity acceleration", 0.0, cond=:(g>=0))
]
@doc docstring(AcousticModalAnalysis_params) AcousticModalAnalysis


mutable struct AcousticModalAnalysis<:Analysis
    stressmodel::Symbol # plane stress, plane strain, etc.
    thickness::Float64  # thickness for 2d analyses
    g::Float64 # gravity acceleration
    freqs::Vector{Float64}
    modes::Matrix{Float64}
    
    function AcousticModalAnalysis(; kwargs...)
        args = checkargs(kwargs, AcousticModalAnalysis_params)
        this = new(args.stressmodel, args.thickness, args.g, [], zeros(0,0))
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
@doc docstring(acoustic_modal_solver_params) acoustic_modal_solver!()


function acoustic_modal_solver!(model::FEModel, stage::Stage; kwargs...)
    args     = checkargs(kwargs, acoustic_modal_solver_params)
    nmodes    = args.nmodes
    quiet    = args.quiet

    ctx = model.ctx
    quiet || println(ctx.log, "Modal analysis for acoustic systems")

    # get dofs organized according to boundary conditions
    dofs, nu    = configure_dofs!(model, stage.bcs)

    ndofs       = length(dofs)
    model.ndofs = length(dofs)
    if !quiet
        println(ctx.log, "unknown dofs: $nu")
        println(ctx.info, "unknown dofs: $nu")
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
    save(model, joinpath(ctx.outdir, "$(ctx.outkey)-0.vtu"), quiet=true)

    model.ctx.inc = 1
    model.ctx.ΔT  = 1.0
    
    # save modes
    for i in 1:nmodes
        U = V[:,i] # modal displacements
        
        for (k,dof) in enumerate(dofs[1:nu])
            dof.vals[dof.name] = U[k]
        end
        model.ctx.T = i/nmodes
        model.ctx.out = i
        update_output_data!(model)
        save(model, joinpath(ctx.outdir, "$(ctx.outkey)-$i.vtu"), quiet=true)
    end

    # reset displacement values
    for dof in dofs
        dof.vals[dof.name] = 0.0
    end

    model.ctx.freqs = ω
    model.ctx.modes = V

    return success()

end

