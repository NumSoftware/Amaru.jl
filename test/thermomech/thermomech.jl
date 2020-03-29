using Amaru
#using Test

# Mesh generation
blocks = [Block([0 0; 1 2], nx = 2, ny = 4, tag = "solids")]

msh = Mesh(blocks, verbose = true)

# Finite element analysis

# Analysis data
k     = 0.0502  # thermal conductivity kW/m/K
rho   = 7.8   # water specific weight Ton/m3
cv    = 486.0  # specific heat (capacity) kJ/Ton/K
E     = 200e6 # kPa
nu    = 0.3
alpha = 1.2e-5 #  thermal expansion coefficient  1/K or 1/°C

materials = ["solids" => ElasticSolidThermo(
    E     = E,
    nu    = nu,
    k     = k,
    rho   = rho,
    cv    = cv,
    alpha = alpha,
)]
dom = Domain(msh, materials)

loggers = [:(y == 1) => NodeGroupLogger("book3.dat")]
setloggers!(dom, loggers)

bcs = [
    :(x == 0) => NodeBC(ut = 10.0),
    :(y == 2) => NodeBC(ut = 20.0),
    :(y == 0) => NodeBC(ux = 0, uy = 0),
    :(x == 1) => NodeBC(fx = 100.0),
]

tm_solve!(
    dom,
    bcs,
    end_time = 300000.0,
    tol      = 0.1,
    nincs    = 10,
    verbose  = true,
    nouts    = 10,
)
# Output
#save(dom, "dom3.vtk")
#save(log1, "book3.dat")


