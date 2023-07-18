using Amaru
using Test

# mesh generation

bl  = Block( [0 0; 0.2 0.1], nx=3, ny=1, cellshape=QUAD4)
msh = Mesh(bl)
iptag!(msh.elems[2], 100)
tag!(msh.elems[2], 100)
#tag!(msh.elems[5], 100)
#tag!(msh.elems, 100)
tag!(msh.faces[:(x>0 && x<0.2 && y==0.1)], 10)
tag!(msh.nodes[3], 3)
tag!(msh.nodes[5], 5)

# finite element analysis

E = 27.e6

mats = [
    MaterialBind(:solids, LinearElastic << (E=E, nu=0.0)),
    MaterialBind(100, SmearedCrack << (E=E, nu=0.2, ft=2.4e3, mu=1.4, alpha=1.0, wc=1.7e-4, ws=1.85e-5, softcurve="hordijk" ) ),
]

# Loggers
log_ip = IpLogger(100)
log_face = FaceLogger(10)
log_n3 = NodeLogger(3)
log_n5 = NodeLogger(5)
loggers = [log_ip, log_face, log_n3, log_n5]

#model = FEModel(msh, mats, stressmodel="plane-stress", thickness=1.0)
model = FEModel(msh, mats, loggers)

# Boundary conditions
bcs = [
       NodeBC(:(x==0 && y==0 ), :(ux=0, uy=0 )),
       NodeBC(:(x==0.2 && y==0 ), :(uy=0 )),
       SurfaceBC(10, :(uy=8*1e-4)),
       #BC(:node, :(x<0.1), :(ux=0, uy=0 )),
       #BC(:node, :(x>0.1), :(ux=1e-5)),
       #BC(:node, :(x>0.1), :(ux=0, uy=8.0*1.7e-4)),
       #BC(:node, :(x>0.1), :(ux=8.0*1.7e-4, uy=8.0*1.7e-4)),

       #BC(:face, :(x==0), :(ux=0, uy=0 )),

       #BC(:face, :(x==0.2), :(ux=8.0*1.7e-4, uy=0)),
       #BC(:face, :(x==0.2), :(uy=8*1.7e-4, ux=0)),

       #BC(:face, :(x==0.2), :(ux=8.0*1.7e-4, uy=8.0*1.7e-4)),
      ]


@test solve!(model, autoinc=true, nincs=100, maxits=4, tol=0.1, quiet=true, scheme=:FE, nouts=50, maxincs=0).success


if @isdefined(makeplots) && makeplots
    using PyPlot
    tab = log_face.table
    plot(log_n5.table[:ux] - log_n3.table[:ux], log_face.table[:fy],  marker="o", color="blue")
    show()

    #tab = log_ip.table
    #plot(tab[:w], tab[:s1], marker="o", color="blue")
    #show()
    #plot(tab[:up], tab[:s3], marker="o", color="blue")
    #show()
    #plot(tab[:s1], tab[:s3], marker="o", color="blue")
    #show()
end
