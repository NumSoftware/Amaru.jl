# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Mazars

mutable struct MazarsIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Tensor2
    ε::Tensor2
    lastΔε::Tensor2
    φt::Float64
    φc::Float64
    φ::Float64  # damage
    ε̅max::Float64
    D::Array{Float64,2}
    function MazarsIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.lastΔε = zeros(6)
        this.φ = 0.0
        this.ε̅max = 0.0
        this.D = zeros(6,6)
        return this
    end
end

mutable struct Mazars<:Material
    E ::Float64
    nu::Float64
    ε̅0::Float64
    At::Float64
    Bt::Float64
    Ac::Float64
    Bc::Float64
    ρ::Float64
    De::Tensor4
    invDe::Tensor4

    function Mazars(prms::Dict{Symbol,Float64})
        return Mazars(;prms...)
    end

    function Mazars(;E=NaN, nu=0.0, eps0=NaN, At=NaN, Bt=NaN, Ac=NaN, Bc=NaN, rho=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert At>0.0
        @assert Ac>0.0
        @assert Bt>0.0
        @assert Bc>0.0
        @assert rho>=0.0

        this     = new(E, nu, eps0, At, Bt, Ac, Bc, rho=0.0)
        this.De  = calcDe(E, nu, :general)
        this.invDe  = inv(this.De)
        return this 
    end
end

# Returns the element type that works with this material model
matching_elem_type(::Mazars) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::Mazars, shared_data::SharedAnalysisData) = MazarsIpState(shared_data)

function set_state(ipd::MazarsIpState; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig.*V2M
    else
        if length(sig)!=0; error("Mazars: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*V2M
    else
        if length(eps)!=0; error("Mazars: Wrong size for strain array: $eps") end
    end
end


# special functions
pos(x) = (abs(x)+x)/2.0
neg(x) = (-abs(x)+x)/2.0


function calcD2(mat::Mazars, ipd::MazarsIpState)
    # There is something wrong with the derivatives here

    if ipd.φ <= 0.0
        # Elastic constitutive matrix
        return mat.De
    else
        # Principal stresses and principal directions
        σp = principal(ipd.σ)
        σp = [ σp; zeros(3) ]
        # Eigen vectors
        p1 = V[:, 1]
        p2 = V[:, 2]
        p3 = V[:, 3]
        P  = [ matrix2Mandel(p1*p1'), matrix2Mandel(p2*p2'), matrix2Mandel(p3*p3') ]

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        Dei= inv(mat.De)
        εt = Dei*σt
        εc = Dei*σc
        εv = sum(pos.(εt)) + sum(pos.(εc))

        # Equivalent strain scalar
        εp = principal(ipd.ε)
        ε̅ = √sum( pos(εp[i])^2 for i=1:3 )

        # Tensile and compression damage weights
        αt =  sum(pos.(εt))/εv
        αc =  sum(pos.(εc))/εv

        # Constitutive matrix calculation
        dφtdε̅ = (1.0-mat.At)*mat.ε̅0*ε̅^-2 + mat.At*mat.Bt*exp( -mat.Bt*(ε̅-mat.ε̅0) )
        dφcdε̅ = (1.0-mat.Ac)*mat.ε̅0*ε̅^-2 + mat.Ac*mat.Bc*exp( -mat.Bc*(ε̅-mat.ε̅0) )
        dε̅dε  = sum( pos(εp[i])*(1+sign(εp[i]))*P[i] for i=1:3 ) / (2*ε̅)
        #dε̅dε  = [ dε̅dε; zeros(3) ]
        dφdε  = (αt*dφtdε̅ + αc*dφcdε̅) * dε̅dε

        #@show dφdε'*mat.De
        D     = (1.0 - ipd.φ)*mat.De - (dφdε'*mat.De)'*ipd.ε'

        return D
    end
end


function calcDsecante(mat::Mazars, ipd::MazarsIpState)
    if ipd.φ <= 0.0
        return mat.De
    else
        D = (1.0 - ipd.φ)*mat.De 
        return D
    end
end

function calcD(mat::Mazars, ipd::MazarsIpState)
    if ipd.φ <= 0.0
        return mat.De
    else
        D = zeros(6,6)
        δ = mat.ε̅0/20  # 20 seems to work better
        κ = 1.0e-15    # important because some components may be zero
        δε = sign.(ipd.lastΔε .+ κ)*δ
        #δε = ones(6)*δ

        Δε = zeros(6)
        for i=1:6
            Δε[i] = δε[i]
            if i>1; Δε[i-1] = 0.0 end
            Δσ = mazars_Δσ(mat, ipd, Δε)
            Di = Δσ/Δε[i]
            D[:,i] = Di
        end
        return D
    end
end


function mazars_Δσ(mat::Mazars, ipd::MazarsIpState, Δε::Vect)::Vect
    # Auxiliary function to aid the calculation of numerical D matrix
    ε = ipd.ε + Δε    # Total strain
    εp = principal(ε) # Principal stresses tensor

    # Equivalent strain scalar
    ε̅ = √sum( pos(εp[i])^2 for i=1:3 )
    ε̅max = max(ipd.ε̅max, mat.ε̅0)

    if ε̅ < ε̅max  # linear-elastic increment
        φ = ipd.φ
    else # increment with damage
        σp = principal(ipd.σ)
        σp = [ σp; zeros(3) ]

        # Damage calculation
        φt = 1.0 - (1-mat.At)*mat.ε̅0/ε̅ - mat.At/exp(mat.Bt*(ε̅-mat.ε̅0))
        φc = 1.0 - (1-mat.Ac)*mat.ε̅0/ε̅ - mat.Ac/exp(mat.Bc*(ε̅-mat.ε̅0))
        #φt = max(ipd.φt, φt)
        #φc = max(ipd.φc, φc)

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        εt = mat.invDe*σt
        εc = mat.invDe*σc

        εv = sum(pos(εt[i]) + pos(εc[i]) for i=1:3)

        # Tensile and compressive damage weights
        αt =  sum(pos(εt[i]) for i=1:3)/εv
        αc =  sum(pos(εc[i]) for i=1:3)/εv
        if εv==0.0
            αt = αc = 0.5
        end

        # Damage variable
        φ = αt*φt + αc*φc
    end

    # Total stress and stress increment
    σ = (1.0 - φ)*mat.De*ε

    Δσ    = σ - ipd.σ
    return Δσ 
    
end


function stress_update(mat::Mazars, ipd::MazarsIpState, Δε::Array{Float64,1})
    if !isnan(Δε[1])
        ipd.lastΔε = copy(Δε)
    end
    σini  = ipd.σ
    ipd.ε = ipd.ε + Δε

    E  = mat.E
    nu = mat.nu

    # Principal stresses tensor
    εp = principal(ipd.ε)

    # Equivalent strain scalar
    ε̅ = √sum( pos(εp[i])^2 for i=1:3 )
    ipd.ε̅max = max(ipd.ε̅max, mat.ε̅0)

    if ε̅ < ipd.ε̅max  # linear-elastic increment
        ipd.σ = (1.0 - ipd.φ)*mat.De*ipd.ε
    else # increment with damage: A previous elastic step is desired
        ipd.ε̅max = ε̅

        # Principal stresses and principal directions
        σp = principal(ipd.σ)
        σp = [ σp; zeros(3) ]

        # Damage calculation
        φt = 1.0 - (1-mat.At)*mat.ε̅0/ε̅ - mat.At/exp(mat.Bt*(ε̅-mat.ε̅0))
        φc = 1.0 - (1-mat.Ac)*mat.ε̅0/ε̅ - mat.Ac/exp(mat.Bc*(ε̅-mat.ε̅0))

        #ipd.φt = max(ipd.φt, φt)
        #ipd.φc = max(ipd.φc, φc)
        ipd.φt = φt
        ipd.φc = φc

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        εt = mat.invDe*σt
        εc = mat.invDe*σc

        εv = sum(pos(εt[i]) + pos(εc[i]) for i=1:3)

        # Tensile and compressive damage weights
        αt =  sum(pos(εt[i]) for i=1:3)/εv
        αc =  sum(pos(εc[i]) for i=1:3)/εv
        if εv==0.0
            αt = αc = 0.5
        end

        # Damage variable
        φ = αt*ipd.φt + αc*ipd.φc
        ipd.φ = max(φ, ipd.φ)
        #ipd.φ = φ

        # Total stress and stress increment
        ipd.σ = (1.0 - ipd.φ)*mat.De*ipd.ε
    end

    Δσ    = ipd.σ - σini
    return Δσ 
end

function ip_state_vals(mat::Mazars, ipd::MazarsIpState)
    σ  = ipd.σ
    ε  = ipd.ε
    ndim = ipd.shared_data.ndim
    sr2  = √2.

    if ndim==2;
        return Dict(
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :sxy => σ[4]/sr2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :exy => ε[4]/sr2,
          :p   => sum(σ[1:3])/3.0,
          :dam => ipd.φ,
          :damt => ipd.φt,
          :damc => ipd.φc,
          :eq => ipd.ε̅max )
      else
        return Dict(
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :sxy => σ[4]/sr2,
          :syz => σ[5]/sr2,
          :sxz => σ[6]/sr2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :exy => ε[4]/sr2,
          :eyz => ε[5]/sr2,
          :exz => ε[6]/sr2,
          :ev  => trace(ε),
          :p   => trace(σ)/3.0,
          :dam => ipd.φ,
          :damt => ipd.φt,
          :damc => ipd.φc,
          :eq => ipd.ε̅max
          )
      end
end

function elem_and_node_vals(mat::Mazars, elem::Element)
    # call generic method from AbsSolid
    node_vals, e_vals = invoke(elem_and_node_vals, Tuple{AbsSolid, Element}, mat, elem)

    # vector with nodal damage values
    D = node_vals[:dam] 
    Dc = node_vals[:damc] 
    Dt = node_vals[:damt] 

    # fix nodal values outside damage range
    for i=1:length(D)
        if D[i]>1.0; D[i] = 1.0 end
        if D[i]<0.0; D[i] = 0.0 end
        if Dc[i]>1.0; Dc[i] = 1.0 end
        if Dc[i]<0.0; Dc[i] = 0.0 end
        if Dt[i]>1.0; Dt[i] = 1.0 end
        if Dt[i]<0.0; Dt[i] = 0.0 end
    end

    node_vals[:dam] = D
    node_vals[:damc] = Dc
    node_vals[:damt] = Dt

    return node_vals, e_vals
end
