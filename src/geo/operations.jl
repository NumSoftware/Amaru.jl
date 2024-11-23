
export array!

function move!(subpath::SubPath; dx::Real=0.0, dy::Real=0.0, dz::Real=0.0)
    for p in subpath.path.points
        p.coord = p.coord + Vec3(dx, dy, dz)
    end
end

function array!(geo::GeoModel, subpath::SubPath; nx=1, ny=1, nz=1, dx=0.0, dy=0.0, dz=0.0)
    for k in 0:nz-1
        for j in 0:ny-1
            for i in 0:nx-1
                i==j==k==0 && continue
                cp = copy(subpath)
                move!(cp, dx=i*dx, dy=j*dy, dz=k*dz)
                addsubpath!(geo, cp)
            end
        end
    end
end
