using Amaru
using Test

# Mesh generation
coord1 = [0 0; 0.1 0; 0.0905325 0.0405325 ;0.07071  0.07071]
coord2 = [0 0; 0.07071  0.07071; 0.0405325 0.0905325 ;0 0.1]
coord3 = [0.1 0 ;0.2 0; 0 0.2; 0 0.1; 0.15 0; 0.14142 0.14142; 0 0.15; 0.07071  0.07071 ]
coord4 = [0.2 0; 1.0 0; 0 1.0; 0 0.2; 0.6 0; 0.7071 0.7071; 0 0.6; 0.14142 0.14142]

blocks = [Block(coord1, nx = 1, ny = 1, cellshape = TRI3, tag = "solids"),
          Block(coord2, nx = 1, ny = 1, cellshape = TRI3, tag = "solids"),
          Block(coord3, nx = 4, ny = 4, cellshape = QUAD4, tag = "solids"),
          Block(coord4, nx = 4, ny = 4, cellshape = QUAD4, tag = "solids"),
    ]

msh = Mesh(blocks)

# Finite element analysis

# Analysis data
k = 0.0502   # thermal conductivity kW/m/K
rho = 6.0   # water specific weight Ton/m3
cv = 140.79 # specific heat (capacity) kJ/Ton/K
E = 200e6 # kN/m2
nu = 0.3
alpha = 5e-5 #  thermal expansion coefficient  1/K or 1/°C

materials = ["solids" << LinearElasticThermo(
    E = E,
    nu = nu,
    k = k,
    rho = rho,
    cv = cv,
    alpha = alpha,
)]
model = FEModel(msh, materials)

bcs = [
     :(x^2 + y^2 >= $0.99^2) << NodeBC(ut =200.0),
     :(y == 0 ) << NodeBC(uy = 0),
     :(x == 0 ) << NodeBC(ux = 0),
     :(x == 0 ) << NodeBC(k = 0),
     :(y == 0 ) << NodeBC(k = 0),
]

tm_solve!(
    model,
    bcs,
    end_time = 16827.49004*1.0,
    tol = 0.1,
    nincs = 50,
    #nouts = 10,
)
