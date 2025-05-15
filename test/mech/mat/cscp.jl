using Amaru
using Test

# mesh
bls = [
       Block( [0 0 0; 0.1 0.1 0.1], shape=HEX8, nx=1, ny=1, nz=1, tag="solids"),
      ]
msh= Mesh(bls)

fc = -30.87e3
ft = 2.95e3

fc_fun = PathFunction(:M, 0.0, fc)
fc_fun = PathFunction(:M, 0.0, 0.5*fc, :Q, 0, 1*fc, 0.002, fc, :C, 0.004, fc, 0.005, 0.2*fc, 0.01, 0.1*fc )

ft_fun = PathFunction(:M, 0.0, ft)
# ft_fun = PathFunction(:M, 0.0, ft, :C, 0.0001, 0.2*ft, 0.0002, 0.0, 0.0003, 0.0)
# ft_fun = PathFunction(:M, 0.0, ft, :C, 0.0004, 0.2*ft, 0.0006, 0.0, 0.0009, 0.0)
ft_fun = PathFunction(:M, 0.0, ft, :C, 0.00002, 0.2*ft, 0.00005, 0.0, 0.0002, 0.0)

p_fun = PathFunction(:M, 0.0, 2.5*fc)
p_fun = PathFunction(:M, 0.0, 1.5*fc, :Q, 0, 1.9*fc, 0.004, 1.9*fc, :Q, 0.04, 2*fc, 0.04, 2.5*fc)
# p_fun = PathFunction(:M, 0.0, 1.61*fc, :L, 0.004, 1.9*fc, :Q, 0.04, 2*fc, 0.04, 2.5*fc)


# fem domain
mats = [
       "solids" => MechSolid => CSCP => (E=30e6, nu=0.2, alpha=1.55, beta=1.14, fc=fc_fun, ft=ft_fun, pmin=p_fun)
      ]

ana = MechAnalysis()
model = FEModel(msh, mats, ana)

loggers = [
           [0.05, 0.05, 0.05] => IpLogger("cscp.table")
          ]
setloggers!(model, loggers)

mons = [
    [0.05, 0.05, 0.05] => IpMonitor(:(sxx, syy), stop=:( rho<0.3*rho_max))
    ]
setmonitors!(model, mons)


θ = 225
α = θ

U      = 0.0014
ΔT     = 0.01
T      = 0.0
factor = 1.1
dT0    = 0.001

while T<1.0
    ux = cosd(α)*U*ΔT
    uy = sind(α)*U*ΔT
    # boundary conditions
    bcs = [
        x==0 => NodeBC(ux=0),
        y==0 => NodeBC(uy=0),
        z==0 => NodeBC(uz=0),
    ]

    if cosd(θ)!=0
        bcs = [ bcs
            x==0.1 => NodeBC(ux=ux)
        ]
    end
    if sind(θ)!=0
        bcs = [ bcs
            y==0.1 => NodeBC(uy=uy)
        ]
    end

    addstage!(model, bcs)
    res = solve!(model, autoinc=true, tol=0.01, rspan=0.02, dT0=dT0)
    global dT0 = 0.05

    sx = model.loggers[1].table[:σxx][end]
    sy = model.loggers[1].table[:σyy][end]
    β  = (atand(sy, sx) + 360) % 360

    global T = T+ΔT
    global ΔT *= factor

    global α = α + 2*(θ-β)

    @show θ
    @show β

    contains(res.message, "Stop") && break
    contains(res.message, "not converge") && break
    
end


include("plot_cscp.jl")