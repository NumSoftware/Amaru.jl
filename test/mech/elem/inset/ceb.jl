using Amaru
using Test

# Mesh:

bls = [
       Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=1, tag="solids"),
       BlockInset( [0.5 3 0.5; 0.5 6.0 0.5], curvetype="polyline", tag="bars", jointtag="joints"),
      ]

msh = Mesh(bls, verbosity=0)
save(msh, "mesh.vtk")

bar_points  = sort(msh.elems["bars"].nodes, by=get_y)

tag!(bar_points[end], "tip")
tag!(msh.elems["solids"].nodes, "fixed_points")


# Finite elements:

mats = [
        "solids" => ElasticSolid(E=24e3, nu=0.2),
        "bars"   => ElasticRod(E=200e6, A=0.00011),
        "joints" => CebRSJoint(TauM=12, TauR=3, s1=0.001, s2=0.0011, s3=0.004, alpha=0.5, beta=0.5,
                                 ks=(12/0.001)*5, kn=5000, A=0.005)
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

tol = 0.01
# tol = 0.002
scheme = "BE"
scheme = "FE"
# scheme = "ME"
scheme = "Ralston"
nincs=20
verbosity=1
#verbosity=0
#maxits=5
maxits=3
# color="blue"
color="green"
color="red"
#color="black"

@test solve!(dom, bcs, nincs=nincs, autoinc=true, scheme=scheme, tol=tol, maxits=maxits, verbosity=0).success

bcs[2] = "tip" => NodeBC(uy=-0.0001)
@test solve!(dom, bcs, nincs=nincs, autoinc=true, scheme=scheme, tol=tol, maxits=maxits, verbosity=0).success

bcs[2] = "tip" => NodeBC(uy=+0.0006)
@test solve!(dom, bcs, nincs=nincs, autoinc=true, scheme=scheme, tol=tol, maxits=maxits, verbosity=0).success

bcs[2] = "tip" => NodeBC(uy=-0.0005)
@test solve!(dom, bcs, nincs=nincs, autoinc=true, scheme=scheme, tol=tol, maxits=maxits, verbosity=0).success

try
bcs[2] = "tip" => NodeBC(uy=+0.005)
@test solve!(dom, bcs, nincs=nincs, autoinc=true, scheme=scheme, tol=tol, maxits=maxits, verbosity=0).success
catch
end


if Amaru.config.makeplots
    using PyPlot
    log_jnt_ip = loggers[2].second
    tab = log_jnt_ip.table

    cplot([
           (x= tab[:ur], y=tab[:tau], marker="o", color=color)
          ],
          #"test.pdf"
         )

    #log_joint_ips = loggers[3].second
    #book = log_joint_ips.book
    #save(log_joint_ips, "book.dat")
    #for tab in book.tables
        #plot(tab[:y], tab[:tau], marker="o")
    #end
    #show()
end
