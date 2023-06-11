using Amaru
using Test

# Mesh:

bls = [
       Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=1, tag="solids"),
       BlockInset( [0.5 3 0.5; 0.5 6.0 0.5], curvetype="polyline", tag="bars", jointtag="joints"),
      ]

msh = Mesh(bls)
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

model = FEModel(msh, mats)
tag!(model.elems["joints"].ips, "joint_ips")

loggers = [
           "tip"       => NodeLogger(),
           "joint_ips" => IpLogger(),
           "joint_ips" => IpGroupLogger()
           ]

setloggers!(model, loggers)

bcs = [
       "fixed_points" => NodeBC(ux=0, uy=0, uz=0),
       "tip"          => NodeBC(uy=0.0003),
      ]

tol = 0.01
scheme = "Ralston"
nincs=20
maxits=3
color="red"

addstage!(model, bcs, nincs=nincs)
bcs[2] = "tip" => NodeBC(uy=-0.0001)
addstage!(model, bcs, nincs=nincs)
bcs[2] = "tip" => NodeBC(uy=+0.0006)
addstage!(model, bcs, nincs=nincs)
bcs[2] = "tip" => NodeBC(uy=-0.0005)
addstage!(model, bcs, nincs=nincs)
bcs[2] = "tip" => NodeBC(uy=+0.005)
addstage!(model, bcs, nincs=nincs)

@test solve!(model, autoinc=true, scheme=scheme, tol=tol, maxits=maxits).success

if @isdefined(makeplots) && makeplots
    using PyPlot
    log_jnt_ip = loggers[2].second
    tab = log_jnt_ip.table

    cplot(
           (x= tab[:ur], y=tab[:tau], marker="o", color=color)
         )
end
