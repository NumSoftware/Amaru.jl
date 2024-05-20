using Amaru
using Revise

r = 0.40
L = 2.0
th = 0.02
n = 20
nouts = 20

E = 200_000_000_000
fy = 500_000_000
nu = 0.3

#k     = 54.108 # thermal conductivity W/m/K
rho   = 7800    # material specific weight kg/m3
cv    = 439.20  # specific heat (capacity) kJ/Ton/K
alpha = 1.2e-5 # thermal expansion coefficient  1/K or 1/°C


k = [
    20.64624298101154  54.10821643286573
    798.203337367805   27.414829659318634
    1199.193436377706  27.414829659318634    
]

k_fun = PathFunction(:M, 20.64, 54.108, :L, 798.20, 27.41, :L, 1199.19, 27.41)


bl = Block( [0 r 0; L r 0], nx=n, cellshape=LIN3)
mesh = Mesh(bl)
mesh = revolve(mesh, base=[0,0,0], axis=[1,0,0], minangle=0, maxangle=360, n=n)

mats = [:all => TMShell => TMCombined{NLConductivity,LinearElastic} =>  (E=E, fy=fy, nu=nu, thickness=th, k=k_fun, alpha=alpha, rho=rho, cv=cv) ]

ana = ThermomechAnalysis(T0=0.0)

model = FEModel(mesh, mats, ana, outdir="shell")

m = 8
# Finite element model
#addmonitor!(model, [L/2, 0, 0] => IpMonitor(:svm))
#addmonitor!(model, [L, 0, 0] => IpMonitor(:svm))

# addlogger!(model, and(x>L-L/m/2, y^2+z^2>r^2) => IpGroupLogger("pipe-vm-tmshell1-160.dat"))
# addlogger!(model, :(y>$r-$r/$m/2) => IpGroupLogger("pipe-vm-tmshell2-160.dat"))
# addlogger!(model, :(y==$r)=> NodeGroupLogger("pipe-vm-tmshell3-160.dat"))

uz = -0.001
#uz = :(-0.00005*t)
#ux = :(-0.0001*t)
ut = :(20+0.5*t)
#tn = :(-85650000*$th)

bcs =
[
    x==0 => NodeBC(ux=0, uy=0, uz=0, rx=0, ry=0, rz=0),
    x==L => EdgeBC(ux=0, uy=0, uz=uz, rx=0, ry=0, rz=0),
    #x>=0 => SurfaceBC(tn=tn),
    x>=0 => NodeBC(ut=ut),
]

addstage!(model, bcs, nouts=nouts, tspan=150)
#solve!(model, tol=1, dTmax = 0.01, autoinc=true, quiet=false)
solve!(model, tol=1, autoinc=true, quiet=false)