using Amaru
#using Revise
nouts=20

th3 = 0.019
th2 = 0.006
th1 = 0.003
a = 0.100
b = 0.033
c = 0.250
L  = 2.0

# mesh
bl1 = Block( [a/2 a/2 0; a/2 -a/2 0], n=4, cellshape=LIN3, tag="shell1")
bl2 = Block( [a/2 -a/2 0; -a/2 -a/2 0], n=4, cellshape=LIN3, tag="shell1")
bl3 = Block( [-a/2 -a/2 0; -a/2 a/2 0], n=4, cellshape=LIN3, tag="shell1")
bl4 = Block( [-a/2 a/2 0; a/2 a/2 0], n=4, cellshape=LIN3, tag="shell1")

bl5 = Block( [0 b 0; 0 a/2 0], n=4, cellshape=LIN3, tag="shell2")
bl6 = Block( [0 -b 0; 0 -a/2 0], n=4, cellshape=LIN3, tag="shell2")

bls = [bl1,bl2,bl3,bl4,bl5,bl6]

mesh = Mesh(bls,ndim=3)
mesh = extrude(mesh, length=L, n=40, axis=[0,0,1])

 
plate1 = Block([-c/2 -c/2; c/2 c/2], nx=10, ny=10, cellshape=QUAD8, tag="placa")
mesh2 = Mesh(plate1, ndim=3)

plate2 = Block([-c/2 -c/2 L; c/2 -c/2 L; c/2 c/2 L; -c/2 c/2 L], nx=10, ny=10, cellshape=QUAD8, tag="placa")
mesh3 = Mesh(plate2, ndim=3)

mesh = Mesh(mesh,mesh2,mesh3)

node1 = nearest(mesh.nodes, [0,0.025,2])
node1.coord = [0,0.033,2]

node2 = nearest(mesh.nodes, [0,-0.025,2])
node2.coord = [0,-0.033,2]

node3 = nearest(mesh.nodes, [0,0.025,0])
node3.coord = [0,0.033,0]

node4 = nearest(mesh.nodes, [0,-0.025,0])
node4.coord = [0,-0.033,0]

#save(mesh4,"test.vtu")
#error()

E = 210_000_000_000
fy =500_000_000
nu = 0.3

k     = 54 # 45 a 50 thermal conductivity W/m/K
rho   = 7850    # material specific weight kg/m3
cv    = 480  # specific heat (capacity) kJ/Ton/K
alpha = 1.2e-5 # thermal expansion coefficient  1/K or 1/°C

#k = PathFunction(:M, 0, 14.2, :L, 1000, 24.5, :L, 1500, 24.5)
#cv = PathFunction(:M, 0, 500, :L, 800, 580, :L, 1500, 24.5)
#alpha = PathFunction(:M, 0, 1.2e-5, :L, 750, 28)


mats = ["shell1" << TMShell << TMCombined{ConstConductivity,LinearElastic} <<  (E=E, fy=fy, nu=nu, thickness = th1, k=k, alpha = alpha, rho=rho, cv=cv),
"shell2" << TMShell << TMCombined{ConstConductivity,LinearElastic} <<  (E=E, fy=fy, nu=nu, thickness = th2, k=k, alpha = alpha, rho=rho, cv=cv),
"placa"<< TMShell << TMCombined{ConstConductivity,LinearElastic} <<  (E=E, fy=fy, nu=nu, thickness = th3, k=k, alpha = alpha, rho=rho, cv=cv) ]

#mats = ["shell1" << TMShell << TMCombined{ConstConductivity,VonMises} <<  (E=E, fy=fy, nu=nu, thickness = th1, k=k, alpha = alpha, rho=rho, cv=cv),
#       "shell2" << TMShell << TMCombined{ConstConductivity,VonMises} <<  (E=E, fy=fy, nu=nu, thickness = th2, k=k, alpha = alpha, rho=rho, cv=cv) ]


#ana = MechAnalysis(stressmodel="3d")
#ana = ThermomechAnalysis(T0=0.0, stressmodel="3d")
#save(mesh, "test.vtu")
#error()

# Finite element model
#addmonitor!(model, [L/2, 0, 0] => IpMonitor(:svm))
#addmonitor!(model, [L, 0, 0] => IpMonitor(:svm))

  ana = ThermomechAnalysis(T0=0.0)
  model = FEModel(mesh, mats, ana, outdir="shell/temp-100")

addlogger!(model, and(y==a/2, x==0, z==L-0.1)=> NodeLogger("column-data.dat"))
addmonitor!(model, and(y==a/2, x==0, z==L-0.1) => NodeMonitor(:uz))

ut1 =:20*t/1000
fz = -110000*t/10
ut2 =:(20 + 345*log10(8*t/60+1))

bcs1 =
[
    z>=0 => EdgeBC(uy=0, uz=0, ux=0, rx=0, ry=0, rz=0, ut=ut1),
]
addstage!(model, bcs1, nouts=1, tspan=1000)
solve!(model, tol=0.1, rspan=0.01, autoinc=true, quiet=false)

bcs2 =
[
    z==0 => EdgeBC(uy=0, uz=0, ux=0),
    #z==L => EdgeBC(uz=0),
    and(z==L, x==0, y==0) => NodeBC(fz=fz),
    z>=0 => EdgeBC(ut=20),
]
addstage!(model, bcs2, nouts=1, tspan=10)
solve!(model, tol=0.1, rspan=0.01, autoinc=true, quiet=false)


bcs3 =
[
    z==0 => EdgeBC(uy=0, uz=0, ux=0),
    #z==L => EdgeBC(uz=0),
    and(z==L, x==0, y==0) => NodeBC(fz=-110000),
    #z>=0 => EdgeBC(ut=20),
    and(z>=L/4, z<=3*L/4) => EdgeBC(ut=ut2),
    #x==0 => EdgeBC(ut=ut),
]
addstage!(model, bcs3, nouts=1, tspan=7*60)
solve!(model, tol=0.1, rspan=0.01, autoinc=true, quiet=false)
