using Amaru
using Test

# Mesh:

bls = [
       Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=5, nz=1, tag="solids"),
       BlockInset( [0.5 3 0.5; 0.5 6.0 0.5], curvetype="polyline", tag="bars", jointtag="joints"),
      ]

msh = Mesh(bls, verbose=true)

bar_points  = sort(msh.cells["bars"].points, by=get_y)

tag!(bar_points[end], "tip")
tag!(msh.cells["solids"].points, "fixed_points")


# Finite elements:

mats = [
        "solids" => ElasticSolid(E=24e3, nu=0.2),
        "bars"   => ElasticRod(E=200e6, A=0.00011),
        "joints" => CEBJoint1D(TauM=12, TauR=3, s1=0.001, s2=0.0011, s3=0.004, alpha=0.5, beta=0.5,
                                 ks=(12/0.001)*5, kn=50000, A=0.005)
       ]

dom = Domain(msh, mats)
tag!(dom.elems["joints"].ips, "joint_ips")

loggers = [
           "tip"       => NodeLogger(),
           "joint_ips" => IpLogger(),
           "joint_ips" => IpGroupLogger()
           ]

setloggers!(dom, loggers)

bcs = [
       "fixed_points" => NodeBC(ux=0, uy=0, uz=0),
       "tip"          => NodeBC(uy=0.0003),
      ]


@test solve!(dom, bcs, nincs=10, autoinc=true, verbose=false)

bcs[2] = "tip" => NodeBC(uy=-0.0001)
@test solve!(dom, bcs, nincs=10, autoinc=true, verbose=false)

bcs[2] = "tip" => NodeBC(uy=+0.0006)
@test solve!(dom, bcs, nincs=10, autoinc=true, verbose=false)

bcs[2] = "tip" => NodeBC(uy=-0.0005)
@test solve!(dom, bcs, nincs=10, autoinc=true, verbose=false)

bcs[2] = "tip" => NodeBC(uy=+0.005)
@test solve!(dom, bcs, nincs=20, autoinc=true, verbose=false)

if Amaru.config.makeplots
    using PyPlot
    log_jnt_ip = loggers[2].second
    tab = log_jnt_ip.table
    plot(tab[:ur], tab[:tau], marker="o", color="blue")
    show()

    #log_joint_ips = loggers[3].second
    #book = log_joint_ips.book
    #save(log_joint_ips, "book.dat")
    #for tab in book.tables
        #plot(tab[:y], tab[:tau], marker="o")
    #end
    #show()
end
