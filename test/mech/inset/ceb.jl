using Amaru
using Test

# Mesh:

bls = [
       Block3D( [0 0 0; 1.0 6.0 1.0], nx=1, ny=5, nz=1),
       BlockInset( [0.5 3 0.5; 0.5 6.0 0.5], curvetype="polyline"),
      ]

msh = Mesh(bls, verbose=true)

bar_points  = sort(msh.cells[:lines][:points], by=get_y)
cell_points = msh.cells[:solids][:points]

#tip = msh.cells[:lines][:nodes][:(x==0.5)][1]

tag!(bar_points[end], "tip")
tag!(cell_points, "fixed_points")
iptag!(msh.cells[:joints1D], "joint_ips")

# Finite elements:

phi = 30*pi/180; dm=0.15138; c=20.0
mats = [
        MaterialBind(:solids,   ElasticSolid(E=24e3, nu=0.2)),
        #MaterialBind(:lines,    ElasticRod(E=1.e8, A=0.005)),
        MaterialBind(:lines,    ElasticRod(E=200e6, A=0.00011)),
        MaterialBind(:joints1D, CEBJoint1D(TauM=12, TauR=3, s1=0.001, s2=0.0011, s3=0.004, alpha=0.5, beta=0.5,
                                           ks=(12/0.001)*5, kn=50000, A=0.005))
       ]

log_tip     = NodeLogger("tip")
log_jnt_ip  = IpLogger("joint_ips")
log_jnt_ips = IpGroupLogger("joint_ips", by=get_y)
loggers = [ log_tip, log_jnt_ip, log_jnt_ips ]

dom = Domain(msh, mats, loggers)

bc1 = NodeBC("fixed_points", :(ux=0, uy=0, uz=0))
bc2 = NodeBC("tip", :(uy=+0.0003))

@test solve!(dom, [bc1, bc2], nincs=10, autoinc=true, verbose=true)

bc2 = NodeBC("tip", :(uy=-0.0001))
@test solve!(dom, [bc1, bc2], nincs=10, autoinc=true, verbose=true)

bc2 = NodeBC("tip", :(uy=+0.0006))
@test solve!(dom, [bc1, bc2], nincs=10, autoinc=true, verbose=true)

bc2 = NodeBC("tip", :(uy=-0.0005))
@test solve!(dom, [bc1, bc2], nincs=10, autoinc=true, verbose=true)

bc2 = NodeBC("tip", :(uy=+0.005))
@test solve!(dom, [bc1, bc2], nincs=30, autoinc=true, verbose=true)

if Amaru.Debug.makeplots
    using PyPlot
    tab = log_jnt_ip.table
    plot(tab[:ur], tab[:tau], marker="o", color="blue")
    show()

    #book = log_joint_ips.book
    #save(log_joint_ips, "book.dat")
    #for tab in book.tables
        #plot(tab[:y], tab[:tau], marker="o")
    #end
    #show()
end
