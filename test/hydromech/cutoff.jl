using Amaru

# Geometria 
Lrocha    = 30.00;
Hrocha    = 15.00;
Lgap      = 0.50;
Hgap      = 8.00;

blocks = [
    Block( [0 0; (Lrocha/2-Lgap/2) Hrocha], nx=14, ny=15, cellshape=QUAD8, tag="solids"),
    Block( [(Lrocha/2-Lgap/2) 0; (Lrocha/2+Lgap/2) (Hrocha-Hgap)], nx=1, ny=7, cellshape=QUAD8, tag="solids"),
    Block( [(Lrocha/2+Lgap/2) 0; Lrocha Hrocha], nx=14, ny=15, cellshape=QUAD8, tag="solids"),
]

mesh  = Mesh(blocks, verbose=true) 

# Analysis data
k    = 1.0e-8  # permeability
gw   = 10.0    # water specific weight

tag!(mesh.cells[:solids], "solids")

# FEM analysis
mats = [
    "solids" => LinSeep(k=k, gammaw=gw, S=0.0),
]

dom = Domain(mesh, mats)

# Stage 1: pore-pressure stabilization
bcs = [
       :(y==$Hrocha) => NodeBC(uw=0),
      ]

hm_solve!(dom, bcs, end_time=1.0e2, nincs=1, tol=1e-2, nouts=1, verbose=true)

dom.env.t = 0.0

#tag!(dom.elems[:solids][:nodes][:(x==1.0 && y==1.0 && z==1.0)], "input")

times = [ 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1 ]


# Stage 2: volume application
bcs = [
       :(x<=$Lrocha/2 && y==$Hrocha) => FaceBC(uw=10*gw),
       :(x>=$Lrocha/2 && y==$Hrocha) => FaceBC(uw=0),
      ]

for t in times[1:end]
    hm_solve!(dom, bcs, end_time=t, nincs=10, tol=1e-2, nouts=1, verbose=true)
end

save(dom, "dom1.vtk")
