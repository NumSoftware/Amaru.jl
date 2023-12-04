using Amaru
#using Test

# Mesh generation
blocks = [Block([0 0; 1 2], nx = 4, ny = 4, tag = "solids")]

msh = Mesh(blocks)

# Finite element analysis

# Analysis data
k     = 0.0502   # thermal conductivity kW/m/K
rho   = 7.8   # water specific weight Ton/m3
cv    = 510.7896  # specific heat (capacity) kJ/Ton/K
E     = 200e6 # kN/m2 = kPa
nu    = 0.3
alpha = 1.2e-5 #  thermal expansion coefficient  1/K or 1/°C
materials = ["solids" => LinearElasticThermo(
    E = E,
    nu = nu,
    k = k,
    rho = rho,
    cv = cv,
    alpha = alpha,
)]
model = FEModel(msh, materials)

log1 = NodeGroupLogger()
loggers = [:(x == 0.5 && y==1.0) => log1]
setloggers!(model, loggers)

bcs = [
     :(y == 0) => NodeBC(ut = 100.0),
     :(y == 2) => NodeBC(ut = 20.0),
     #:(y >= 0) => NodeBC(ux = 0, uy = 0),
     :(y == 0) => NodeBC(ux = 0, uy = 0),
     :(y == 2) => EdgeBC(ty = 100), # kPa
]

tm_solve!(
    model,
    bcs,
    end_time = 20000.0,
    tol = 0.01,
    nincs = 20,
    nouts = 10,
)

