
using Amaru
using Test

lx = 0.1
nx = 6
ny = 3
nz = 3

# mesh
bls = [
       Block3D( [0 0 0; lx 0.1 0.1], nx=nx, ny=ny, nz=nz, tag="solids"),
      ]
msh= Mesh(bls)

mats = [
        "solids" << MechSolid() <<  DamageConcrete(E=27e6, nu=0.2, ft=2.4e3, GF=73e-3, fc=-30e3, epsc=-0.002),
       ]

tol = 1.0
nincs = 50

T = Amaru.DTable([:s1,:s2])

ana = MechAnalysis()
model = FEModel(msh, mats, ana)
tag!(model.elems[1].ips[1], "ip")
logs = [
        "ip" << IpLogger(),
        :(x==$lx) << FaceLogger(),
       ]

setloggers!(model, logs)

bcs = [
       :(x==0) << NodeBC(ux=0),
       :(y==0) << NodeBC(uy=0),
       :(z==0) << NodeBC(uz=0),
       #:(x==0 && y==0 && z==0) << NodeBC(uy=0, uz=0),
       :(x==$lx) << SurfaceBC(ux=0.0001),
       #:(x==0.1 && y==0.1 && z==0.1) << NodeBC(ux=0.0001),
      ]

solve!(model, autoinc=true, scheme=:ME, nincs=nincs, nouts=5, tol=tol, maxits=2).success

t = get_segment_data(model, [0,0,0], [0.1,0.0,0.0], "data.dat")
using PyPlot
plot(t["s"], t["sxx"], "-")
error()

if @isdefined(makeplots) && makeplots
    using PyPlot
    tab1 = logs[1].second.table
    tab2 = logs[2].second.table

    #plot(tab2[:ux], tab1[:sxx], "-o")
    #plot(tab1[:w], tab1[:sxx], "-o")
    #plot(tab2[:ux], tab2[:fx], "-o")
    plot(tab1[:w], tab2[:fx], "-o")
end


