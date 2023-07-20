using Amaru
using Test

R = 22.27
zmin = R*sind(26.67)

# Finite element model
a = 4
b = 16

node = Node(R, 0, 0)
mesh = revolve(node, minangle=26.67, maxangle=90, base=[0,0,0], axis=[0,-1,0], n=a)
mesh = revolve(mesh, base=[0,0,0], axis=[0,0,1], angle=360, n=b)
tag!(mesh.elems, "shell")

# Finite element analysis

# Analysis data
k     = 0.0502 # thermal conductivity kW/m/Ka
rho   = 7.8    # material specific weight Ton/m3
cv    = 486.0  # specific heat (capacity) kJ/Ton/K
E     = 200_000_000  # kPa
nu    = 0.3
alpha = 1.2e-5 # thermal expansion coefficient  1/K or 1/°C
th    = 0.1

# materials = ["shell"<< TMShell << ElasticShellThermo << (E=E, nu=nu, k=k, alpha=alpha,thickness=th, rho=rho, cv=cv) ]
materials = ["shell"<< TMShell << TMCombined{ConstConductivity, LinearElastic} << (E=E, nu=nu, k=k, alpha=alpha, thickness=th, rho=rho, cv=cv) ]
# materials = ["shell"<< TMShell << TMCombined{ConstConductivity, VonMises} << (E=E, nu=nu, k=k, alpha=alpha, thickness=th, rho=rho, cv=cv, H=0.0, fy=100000.0) ]

ana = ThermomechAnalysis(T0=0.0)
model = FEModel(mesh, materials, ana)
addmonitor!(model, :(x==0 && y==0 && z==$R) => NodeMonitor(:uz))
addmonitor!(model, :(x==0 && y==0 && z==$R) => NodeMonitor(:ut))

bcs = [

	:(z==$zmin) => NodeBC(ux=0, uy=0, uz=0, rx=0, ry=0, rz=0),
    :(x==0 && y==0 && z==$R) => NodeBC(fz=10000),
    :(z==$zmin) => NodeBC(ut = 50),
]

addstage!(model, bcs, tspan=20_000_000, nincs=10)
solve!(model, autoinc=true, tol=0.1,)