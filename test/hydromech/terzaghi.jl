using Amaru

# Mesh generation

blocks = [
    Block3D( [0 0 0; 1 1 10], nx=1, ny=1, nz=10, shape=HEX8),
]

mesh = Mesh(blocks, verbose=true)

# Finite element analysis

# Analysis data
load = 10.0   
k    = 1.0E-6  # permeability
E    = 5000.0  # Young modulus
nu   = 0.25    # Poisson
gw   = 10.0    # water specific weight
hd   = 10.0    # drainage height
mv   = (1+nu)*(1-2*nu)/(E*(1-nu))
cv   = k/(mv*gw)   # consolidation coefficient
T = [ 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0 ]
times = T*(hd^2/cv)*2
@show times
t1 = times[1]/10
times .+= t1
@show times

materials = [
    MaterialBind(:solids, ElasticSolidLinSeep(E=E, nu=nu, k=k, gw=gw) ),
]

logger = [
    GroupLogger(:nodes, :(x==0 && y==0) ),
]

dom = Domain(mesh, materials, logger)

# Stage 0: hydrostatic pore-pressure
tlong = 1000000.0

bcs = [
    BC(:node, :all, :(ux=0, uy=0) ),
    BC(:node, :(z==0), :(ux=0, uy=0, uz=0) ),
    BC(:node, :(z==10), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=tlong, nincs=1, tol=1, saveincs=true, verbose=true)

#exit()
dom.shared_data.t = 0.0

# Stage 1: loading
#pt(t) = t>=t1? -load : -load/t1*t
pt(t) = t>=t1? -load : 0.0

bcs = [
    BC(:node, :all, :(ux=0, uy=0) ),
    BC(:node, :(z==0), :(ux=0, uy=0, uz=0) ),
    BC(:face, :(z==10), :(tz=$pt(t)) ),
    BC(:node, :(z==10), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=t1, saveincs=true, verbose=true)


# Stage 2: draining

bcs = [
    BC(:node, :all, :(ux=0, uy=0) ),
    BC(:node, :(z==0), :(ux=0, uy=0, uz=0) ),
    BC(:face, :(z==10), :(tz=$pt(t)) ),
    BC(:node, :(z==10), :(uw=0.) ),
]

for t in times[1:4]
    @show t
    hm_solve!(dom, bcs, end_time=t, nincs=1, tol=1, nouts=1, saveincs=true, verbose=true)
end

# Output

function calc_Ue(Z, T)
    sum = 0.0
    for i=0:4
		M = pi/2*(2*i+1)
		sum = sum + 2/M*sin(M*Z)*exp(-M^2*T)
    end
    return sum
end

using PyPlot
save(logger[1], "book.dat")

book = logger[1].book

Uwini = book.tables[2][:uw]
Z     = 1 - book.tables[2][:z]/hd
@show Uwini
for (i,table) in enumerate(book.tables[4:end])
    dUw = (table[:uw] - Uwini)/load
    plot(dUw, Z, "-o")
end

for Ti in T
    Ue = calc_Ue.(Z, Ti)
    plot(Ue, Z, "k")
end

@show book.tables[end][:uw]

ylim(1,0)
show()

#save(dom, "dom.vtk")

