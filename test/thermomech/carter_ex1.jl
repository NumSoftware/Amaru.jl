using Amaru
using Test

# Mesh generation

blocks = [Block([0 0.0; 0.5 0.5], nx = 2, ny = 2, cellshape = QUAD4, tag = "solids"),
          Block([0 0.5; 0.5 0.75], nx = 2, ny = 2, cellshape = QUAD4, tag = "solids"),
          Block([0 0.75; 0.5 0.875], nx = 2, ny = 2, cellshape = QUAD4, tag = "solids"),
          Block([0 0.875; 0.5 1], nx = 2, ny = 2, cellshape = QUAD4, tag = "solids"),
          ]

msh = Mesh(blocks, printlog=true)

# Finite element analysis

# Analysis data
k = 0.0502   # thermal conductivity kW/m/K
rho = 6.0   # water specific weight Ton/m3
cv = 151.227*10^3 # specific heat (capacity) kJ/Ton/K
E = 200e6 # kN/m2
nu = 0.3
alpha = 5e-5 #  thermal expansion coefficient  1/K or 1/°C

materials = ["solids" => ElasticSolidThermo(
    E = E,
    nu = nu,
    k = k,
    rho = rho,
    cv = cv,
    alpha = alpha,
)]

dom = Domain(msh, materials)

log1 = NodeGroupLogger()
loggers = [:(x == 0.25) => log1]
setloggers!(dom, loggers)


bcs = [
     :(x == 0) =>  NodeBC(k = 0),
     :(x == 0.5) =>  NodeBC(k = 0),
     :(y == 0) =>  NodeBC(k = 0),
     :(y == 0) => NodeBC(ux = 0, uy =0),
     :(y == 1) => FaceBC(ut = 100),
]

tm_solve!(
    dom,
    bcs,
    end_time = 13547.80876*1.0,
    tol = 0.01,
    nincs = 50,
    printlog=true,
    #nouts = 10,
)
