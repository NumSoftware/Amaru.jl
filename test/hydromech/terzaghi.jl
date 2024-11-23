using Amaru
using Test

# Mesh generation


blocks = [
    Block( [0 0 0; 1 1 10], nx=1, ny=1, nz=10, cellshape=HEX8, tag="solids"),
]

msh = Mesh(blocks)

# Finite element analysis


# Analysis data
load   = 10.0
k      = 1.0E-5  # permeability
E      = 5000.0  # Young modulus
nu     = 0.25    # Poisson
gw     = 10.0    # water specific weight
hd     = 10.0    # drainage height
mv     = (1+nu)*(1-2*nu)/(E*(1-nu))
cv     = k/(mv*gw)   # consolidation coefficient
T      = [ 0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0 ]
times  = T*(hd^2/cv)
lapses = diff(times)
t1     = lapses[1]/10

materials = [
    # "solids" => HMSolid => LinearElasticSeep => (E=E, nu=nu, k=k)
    "solids" => HMSolid => HMCombined{ConstPermeability,LinearElastic} => (E=E, nu=nu, k=k)
]

ana = HydromechAnalysis(gammaw=gw)
model = FEModel(msh, materials, ana)

log1 = NodeGroupLogger()
loggers = [
    :(x==0 && y==0) => log1
]
setloggers!(model, loggers)


# Stage 1: hydrostatic pore-pressure


tlong = 10000*hd^2/cv
bcs = [
    :(z==0)  => NodeBC(ux=0, uy=0, uz=0),
    :(x>=0)  => NodeBC(ux=0, uy=0),
    :(z==10) => NodeBC(uw=0.),
]
addstage!(model, bcs, tspan=tlong, nincs=2, nouts=1)
solve!(model, tol=1e-2)
model.ctx.t = 0.0

# Stage 2: loading
bcs = [
    :(z==0)  => NodeBC(ux=0, uy=0, uz=0),
    :(x>=0)  => NodeBC(ux=0, uy=0),
    :(z==10) => SurfaceBC(tz=-load),
    :(z==10) => NodeBC(uw=0.),
]
addstage!(model, bcs, tspan=t1, nincs=4, nouts=1)
solve!(model, tol=1e-2)

# Stage 3: draining
bcs = [
    :(z==0)  => NodeBC(ux=0, uy=0, uz=0),
    :(x>=0)  => NodeBC(ux=0, uy=0),
    :(z==10) => SurfaceBC(tz=-load),
    :(z==10) => NodeBC(uw=0.),
]

Uw_vals = [] # A list with porepressure vectors
for tspan in lapses
    addstage!(model, bcs, tspan=tspan, nincs=20, nouts=1)
    solve!(model, tol=1e-2)
    push!( Uw_vals, log1.book[end][:uw] )
end


# Output
if @isdefined(makeplots) && makeplots
    # Terzaghi's 1-d consolidation
    function calc_Ue(Z, T)
        sum = 0.0
        for i in 0:4
            M = pi/2*(2*i+1)
            sum = sum + 2/M*sin(M*Z)*exp(-M^2*T)
        end
        return sum
    end

    using PyPlot


    # numerical curves
    book  = log1.book
    Uwini = book.tables[2][:uw]       # hydrostatic porepressure
    Z     = 1 .- book.tables[2][:z]/hd # normalized depth
    for Uw in Uw_vals
        dUw = (Uw - Uwini)/load       # excess of porepressure
        plot(dUw, Z, "-o")
    end

    # analytical curves
    for Ti in T
        Ue = calc_Ue.(Z, Ti)
        plot(Ue, Z, "k")
    end

    ylim(1,0)
    show()
end
