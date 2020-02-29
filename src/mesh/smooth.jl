# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

# This file contains the code for smoothing meshes

# Returns a matrix with the cell coordinates
function cellcoords(c::Cell)
    ndim = c.shape.ndim
    C = Array{Float64}(undef, length(c.points), ndim)
    for (i,point) in enumerate(c.points)
        C[i,1] = point.x
        C[i,2] = point.y
        if ndim>2
            C[i,3] = point.z
        end
    end
    return C
end

# Basic coordinates are defined considering an area/volume equal to 1.0
# TODO: change to reference_coords
function basic_coords(shape::ShapeType) #check
    if shape == TRI3
        #return √2*[0.0 0; 1 0 ; 0 1]
        a = 2.0/3.0^0.25
        return [ 0.0 0.0; a 0.0; a/2.0 a/2*√3 ]
    end
    if shape == TRI6
        #return √2*[0.0 0; 1 0 ; 0 1; 0.5 0; 0.5 0.5; 0 0.5]
        a = 2.0/3.0^0.25
        h = a/2*√3
        return [ 0.0 0; a 0; a/2.0 a/2*√3; a/2.0 0; 0.75*a h/2.0; 0.25*a h/2.0 ]
    end
    if shape == QUAD4
        return [ 0.0 0.0; 1.0 0.0; 1.0 1.0 ; 0.0 1.0 ]
    end
    if shape == TET4
        a = (6.0*√2.0)^(1.0/3.0)
        return [   0.0        0.0          0.0 ;
                   a         0.0          0.0 ;
                 a/2.0  √3.0/2.0*a          0.0 ;
                 a/2.0  √3.0/6.0*a  1.0/3.0*√6*a ]
    end
    if shape == QUAD8
        return [ 0.0 0; 1 0; 1 1 ; 0 1; 0.5 0; 1 0.5; 0.5 1; 0 0.5 ]
    end
    if shape == QUAD9
        return [ 0.0 0; 1 0; 1 1 ; 0 1; 0.5 0; 1 0.5; 0.5 1; 0 0.5; 0.5 0.5 ]
    end
    if shape == QUAD12
        return [ 0.0 0; 1 0; 1 1 ; 0 1;
        1/3 0.0; 1.0 1/3; 2/3 1.0; 0.0 2/3;
        2/3 0.0; 1.0 2/3; 1/3 1.0; 0.0 1/3 ]
    end
    if shape == HEX8
        return [ 0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 1.0; 1.0 1.0 1.0; 0.0 1.0 1.0 ]
    end
    if shape == HEX20
        return [ 0.0 0.0 0.0;
                 1.0 0.0 0.0;
                 1.0 1.0 0.0;
                 0.0 1.0 0.0;
                 0.0 0.0 1.0;
                 1.0 0.0 1.0;
                 1.0 1.0 1.0;
                 0.0 1.0 1.0;
                 0.5 0.0 0.0;
                 1.0 0.5 0.0;
                 0.5 1.0 0.0;
                 0.0 0.5 0.0;
                 0.5 0.0 1.0;
                 1.0 0.5 1.0;
                 0.5 1.0 1.0;
                 0.0 0.5 1.0;
                 0.0 0.0 0.5;
                 1.0 0.0 0.5;
                 1.0 1.0 0.5;
                 0.0 1.0 0.5 ]
    end
    if shape == WED6
        a = (4.0/√3.0)^(1.0/3.0)
        return [ 0.0  0.0  0.0; a  0.0  0.0; a/2.0  a/2*√3  0.0 ;
                 0.0  0.0   a; a  0.0   a; a/2.0  a/2*√3   a ]
    end
    if shape == WED15
        a = (4.0/√3.0)^(1.0/3.0)
        b = a/2*√3
        return [   0  0    0; a        0     0; a/2    b    0;
                   0  0    a; a        0     a; a/2    b    a;
                 a/2  0    0; a*3/4  b/2     0; a/4  b/2    0;
                 a/2  0    a; a*3/4  b/2     a; a/4  b/2    a;
                   0  0  a/2; a        0   a/2; a/2    b  a/2 ]
    end

    if shape == PYR5
        a = (√2/3)*(1/3)
        return [ 0.0   0.0      0.0;
                   a   0.0      0.0;
                   a     a      0.0;
                 0.0     a      0.0;
                 a/2   a/2   √2/2*a;
                ]
    end

    error("No basic coordinates for shape $(shape.name)")
end

# Returns a rotation matrix for a cell based in their first points
function cell_orientation(cell::Cell)
    shape   = cell.shape
    geo_dim = shape.ndim

    if geo_dim == 2
        p1 = cell.points[1]
        p2 = cell.points[2]

        l1 = p2.x - p1.x
        m1 = p2.y - p1.y
        T1 = [l1, m1]
        l1, m1 = T1/norm(T1)
        return [ l1  -m1; m1  l1 ]
    end

    if geo_dim == 3
        p1 = cell.points[1]
        p2 = cell.points[2]
        if shape == HEX8 || shape == HEX20
            p3 = cell.points[4]
        else
            p3 = cell.points[3]
        end

        T1 = [ p2.x - p1.x, p2.y - p1.y, p2.z - p1.z ]
        T2 = [ p3.x - p1.x, p3.y - p1.y, p3.z - p1.z ]
        T3 = cross(T1, T2)
        T2 = cross(T3, T1) # redefine T2

        T1 = T1/norm(T1)
        T2 = T2/norm(T2)
        T3 = T3/norm(T3)

        return [T1 T2 T3]

    end
    error("cell_orientation: Not implemented for shape ", cell.shape)
end


# Matrix D for the simplified FEM analysis
function matrixD(E::Float64, nu::Float64)
    c = E/((1.0+nu)*(1.0-2.0*nu))
    [ c*(1.0-nu)      c*nu        c*nu             0.0             0.0             0.0
          c*nu   c*(1.0-nu)       c*nu             0.0             0.0             0.0
          c*nu       c*nu    c*(1.0-nu)            0.0             0.0             0.0
           0.0        0.0         0.0   c*(1.0-2.0*nu)            0.0             0.0
           0.0        0.0         0.0             0.0   c*(1.0-2.0*nu)            0.0
           0.0        0.0         0.0             0.0             0.0   c*(1.0-2.0*nu) ]
end


#include("../tools/linalg.jl")

# Matrix B for the simplified FEM analysis
function matrixB(ndim::Int, dNdX::Matx, detJ::Float64, B::Matx)
    nnodes = size(dNdX,2)
    sqr2 = √2.0
    B .= 0.0
    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[1,i]
            B[2,2+j*ndim] = dNdX[2,i]
            B[4,1+j*ndim] = dNdX[2,i]/sqr2; B[4,2+j*ndim] = dNdX[1,i]/sqr2
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[1,i]
            dNdy = dNdX[2,i]
            dNdz = dNdX[3,i]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,1+j*ndim] = dNdy/sqr2;   B[4,2+j*ndim] = dNdx/sqr2
            B[5,2+j*ndim] = dNdz/sqr2;   B[5,3+j*ndim] = dNdy/sqr2
            B[6,1+j*ndim] = dNdz/sqr2;   B[6,3+j*ndim] = dNdx/sqr2
        end
    end

    return detJ
end

function matrixK(cell::Cell, ndim::Int64, E::Float64, nu::Float64)
    nnodes = length(cell.points)

    C = cellcoords(cell)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)

    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    IP = get_ip_coords(cell.shape)
    D = matrixD(E, nu)

    for i=1:size(IP,1)
        R    = vec(IP[i,1:3])
        w    = IP[i,4]

        # compute B matrix
        dNdR = cell.shape.deriv(R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        if detJ < 0.0
            #@error "Negative jacobian determinant in cell" cell=cell.id ip=i coords=C shape=cell.shape.name
            #error()
        end
        matrixB(ndim, dNdX, detJ, B)

        # compute K
        coef = detJ*w
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end
    return K
end

function matrixK2(cell::Cell, ndim::Int64, E::Float64, nu::Float64)
    nnodes = length(cell.points)

    C = cellcoords(cell)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)
    DB = Array{Float64}(undef, 6, nnodes*ndim)

    IP = get_ip_coords(cell.shape)

    D = matrixD(E, nu)
    for i=1:size(IP,1)
        R    = vec(IP[i,1:3])
        w    = IP[i,4]
        detJ = matrixB(cell, ndim, R, C, B)

        @gemm DB = D*B
        coef = detJ*w
        @gemm K += coef*B'*DB
        #K   += B'*D*B*detJ*w
    end
    return K
end

function get_map(c::Cell)
    ndim = c.shape.ndim

    map = Int[]
    for p in c.points
        for i=1:ndim
            push!(map, (p.id-1)*ndim + i)
        end
    end

    return map
end

# Mount global stiffness matrix
function mountKg(mesh::Mesh, E::Float64, nu::Float64, A)
    ndim = mesh.ndim
    ndof = ndim*length(mesh.points)
    R, C, V = Int64[], Int64[], Float64[]

    for c in mesh.cells
        Ke  = matrixK(c, ndim, E, nu)
        map = get_map(c)
        nr, nc = size(Ke)
        for i=1:nr
            for j=1:nc
                push!(R, map[i])
                push!(C, map[j])
                push!(V, Ke[i,j])
            end
        end
    end

    # mount A and A'
    nbc = size(A,1)
    for (i,j,val) in zip(findnz(A)...)
            push!(R, ndof+i)
            push!(C, j)
            push!(V, val)
            push!(R, j)
            push!(C, ndof+i)
            push!(V, val)
    end

    return sparse(R, C, V, ndof+nbc, ndof+nbc)
end

# Check if an array of faces are coplanar
# Returns the normal unitary vector of a face
# If the face is not flat returns [0,0,0]
function normal_to_faces(faces::Array{Cell, 1})
    ndim = 1 + faces[1].shape.ndim
    points = Point[]

    for f in faces
        for p in f.points
            push!(points, p)
        end
    end

    # mount coordinates matrix
    C = Array{Float64}(undef, length(points), ndim)
    for (i,p) in enumerate(points)
        C[i,1] = p.x
        C[i,2] = p.y
        if ndim>2
            C[i,3] = p.z
        end
    end

    C += 1.e-5
    I = ones(size(C,1))
    N = pinv(C)*I # best fit normal

    if norm(C*N - I)<1e-5
        return N/norm(N)
    end

    return zeros(ndim)
end


function faces_normal(faces::Array{Cell,1}, facetol)
    ndim = 1 + faces[1].shape.ndim

    #normals = Set{Array{Float64,1}}()
    normals = Array{Float64,1}[]

    for f in faces

        C = getcoords(f, ndim)

        # move the coordinates to avoid singular case
        # when the regression line/planes crosses the origin
        if ndim==2
            C .+= [pi pi^1.1]
        else
            C .+= [pi pi^1.1 pi^1.2]
        end

        I = ones(size(C,1))
        N = pinv(C)*I # best fit normal

        #if norm(C*N - I) < facetol # TODO: check if points in face are coplanar
            #N = N/norm(N)
            normalize!(N)
            #check if the normal already exists
            if all( [ norm(N-NN)>facetol for NN in normals ])
                push!(normals, N)
            end
        #else
            #return Array{Float64,1}[]
        #end
    end

    return normals
end

# Auxiliary structure for a surface node
mutable struct sNode
    point::Point
    faces::Array{Cell}
    normals
end

# Mount global matrix A
function mountA(mesh::Mesh, fixed::Bool, conds, facetol)
    # get border faces
    ndim = mesh.ndim
    npoints = length(mesh.points)
    surf_cells = get_surface(mesh.cells)
    surf_points, surf_patches = get_patches(surf_cells)
    border_nodes = [ sNode(point, patch, nothing) for (point,patch) in zip(surf_points,surf_patches)]

    local fconds=Function[]
    if conds!= nothing
        for c in conds
            ff = quote
                (x,y,z) -> ($c)
            end
            push!(fconds, eval(ff))
        end
    end

    # find the number of bcs
    n = 0  # number of bcs
    for node in border_nodes
        if conds!= nothing
            p = node.point
            if any( Bool[ ff(p.x, p.y, p.z) for ff in fconds ] )
                node.normals = Array{Float64,1}[]
                n += ndim
                continue
            end
        end

        node.normals = faces_normal(node.faces, facetol)
        nnorm        = length(node.normals)
        if nnorm==1 || nnorm==2
            n += nnorm
        else
            n += ndim
        end
    end


    #  for fixed boundary
    if fixed
        n = length(border_nodes)
        #A = zeros(n*ndim, length(mesh.points)*ndim)
        R, C, V = Int64[], Int64[], Float64[]

        for i=1:n
            for j=1:ndim
                #A[ (i-1)*ndim+j, (border_nodes[i].point.id-1)*ndim+j ] = 1.0
                push!(R, (i-1)*ndim+j )
                push!(C, (border_nodes[i].point.id-1)*ndim+j )
                push!(V, 1.0)
            end
        end
        A = sparse(R, C, V, n*ndim, npoints*ndim)
    else
        # mount matrix A (Lagrange multipliers) according to bcs
        #A = zeros(n, length(mesh.points)*ndim)
        R, C, V = Int64[], Int64[], Float64[]

        baserow = 0
        for node in border_nodes
            basecol = (node.point.id-1)*ndim
            normals = node.normals
            nnorm   = length(normals)

            if nnorm==1 || nnorm==2  # all faces in up to 2 planes
                for i=1:nnorm
                    for j=1:ndim
                        push!(R, baserow+1)
                        push!(C, basecol+j)
                        push!(V, normals[i][j])
                        #A[ baserow+1, basecol+j ] = normals[i][j]
                    end
                    baserow += 1
                end
            else # zero or more than two non coplanar faces
                for j=1:ndim
                    #A[ baserow+j, basecol+j ] = 1.0
                    push!(R, baserow+j)
                    push!(C, basecol+j)
                    push!(V, 1.0)
                end
                baserow += ndim
            end

        end

        A = sparse(R, C, V, n, npoints*ndim)
    end

    return A
end



function rigid_transform(source::Array{Float64,2}, target::Array{Float64,2}, pindexes::Array{Int64,1}=Int[])
    # source: Source coordingates
    # target: Destination reference
    # pindexes: Indexes for prioriy matching
    # C: If provided, will represent the center of rotation
    # e.g. R, T = rigid_transform(BC, C0)
    A, B = copy(source), copy(target)
    #@test size(A) == size(B)

    n = size(A,1)
    ndim = size(A,2)

    # Centralizing both sets of points
    #pindexes = []
    #if length(pindexes)==0
        #cA = mean(A, dims=1)
        #cB = mean(B, dims=1)
    #else
        #cA = mean(A[pindexes,:], dims=1)
        #cB = mean(B[pindexes,:], dims=1)
    #end
    cA = mean(A, dims=1)
    cB = mean(B, dims=1)

    A .-= cA
    B .-= cB

    for i in pindexes
        A[i,:] .*= 10
        #B[i,:] .*= 1000
    end


    # Singular value decomposition
    U, S, V = svd(A'*B)

    # Rotation matrix
    R = V*U'

    # special reflection case
    if det(R) < 0.0
       #println("Reflection detected")
       R[:2] *= -1
    end

    T = cB - cA*R'
    #T = -cA*R'

    return R, T

end


# Mount a vector with nodal forces
function force_bc(mesh::Mesh, E::Float64, nu::Float64, α::Float64, extended::Bool)
    n    = length(mesh.points)
    ndim = mesh.ndim
    Fbc  = zeros(n*ndim)
    if extended
        qmin = minimum(c.quality for c in mesh.cells)
    end

    for c in mesh.cells
        # get coordinates matrix
        np = length(c.points)
        C0 = cellcoords(c)
        V  = abs(cell_extent(c)) # area or volume

        #extended || (V = α*V)

        #R0 = cell_orientation(c)
        s = V^(1.0/ndim) # scale factor

        extended && (s *= α)

        BC = basic_coords(c.shape)*s

        # initial alignment
        #C = BC*R0'
        C = BC

        # align C with cell orientation
        R, d = rigid_transform(C, C0)
        D = repeat( d , np, 1)
        C1 = C*R' + D
        U  = C1 - C0  # displacements matrix

        U  = vec(U')  # displacements vector

        K = matrixK(c, ndim, E, nu)

        if extended
            if qmin<1.0
                F = K*U*((1-c.quality)/(1-qmin))
            else
                F = zeros(length(U))
            end
        else
            F = K*U
        end

        # add forces to Fbc
        for (i,point) in enumerate(c.points)
            pid = point.id
            for j = 1:ndim
                Fbc[(pid-1)*ndim+j] += F[(i-1)*ndim+j]
            end
        end
    end

    return Fbc
end



# Mount a matrix with nodal displacements
function find_disps(mesh::Mesh, patches, extended, α, qmin, in_border)
    n    = length(mesh.points)
    m    = length(mesh.cells)
    ndim = mesh.ndim
    UU = zeros(n,ndim)
    W  = ones(m) # qualities

    for c in mesh.cells
        # get coordinates matrix
        np = length(c.points)
        C0 = cellcoords(c)
        v  = abs(cell_extent(c)) # area or volume
        W[c.id] = v
        #W[c.id] = 1.0
        #W[c.id] = (2.0 - c.quality)
        #W[c.id] = c.quality^0.5

        #R0 = cell_orientation(c)
        s = v^(1.0/ndim) # scale factor

        extended && (s *= α)

        BC = basic_coords(c.shape)*s
        C = BC

        # align C with cell orientation
        pindexes = [ i for i=1:np if in_border[ c.points[i].id] ]
        #@show pindexes
        R, D = rigid_transform(C, C0, pindexes)
        C1 = C*R' .+ D
        U  = C1 - C0  # displacements matrix

        #if extended && qmin<1.0
            #U .= (1-c.quality)/(1-qmin)*U
        #end
        #U .*= 0.666
        U .*= 0.5

        # add to UU
        dmap = [ p.id for p in c.points ]
        UU[dmap, :] .+= W[c.id].*U
    end

    # weighted average
    for node in mesh.points
        patch = patches[node.id]
        sumW = sum( W[c.id] for c in patch )
        UU[node.id, :] ./= sumW
    end

    return UU
end


function str_histogram(hist::Array{Int64,1})
    m = maximum(hist)
    H = round.(Int, hist./m*7)
    chars = [" ","_","▁","▂","▃","▄","▅","▆","▇","█"]
    return "["*join( hist[i]==0 ? " " : chars[H[i]+2] for i=1:length(H) )*"]"
end


function fast_smooth!(mesh::Mesh; verbose=true, alpha::Float64=1.0, target::Float64=0.97,
                 fixed::Bool=false, maxit::Int64=30, mintol::Float64=1e-3, tol::Float64=1e-4,
                 facetol=1e-4, savesteps::Bool=false, savedata::Bool=false, bin::Float64=0.05,
                 filekey::String="smooth", conds=nothing, extended=false, smart=false)

    # tol   : tolerance in change of mesh quality for succesive iterations
    # mintol: tolerance in change of worst cell quality in a mesh for succesive iterations

    #verbose && printstyled("Mesh smoothing:\n", bold=true, color=:cyan)
    verbose && printstyled("Mesh fast-$(smart ? "smart-" : "")cells-fitting smoothing:\n", bold=true, color=:cyan)

    # check for not allowed cells
    for c in mesh.cells
        if c.shape.family != SOLID_SHAPE
            error("smooth!: cells of family $(c.shape.family) are not allowed for smoothing: $(c.shape.name)")
        end
    end

    ndim = mesh.ndim
    npoints = length(mesh.points)

    points, patches = get_patches(mesh)  # key points and corresponding patches

    # get a list of surface nodes (sNodes) that include a list of faces per point
    surf_cells = get_surface(mesh.cells)
    surf_points, surf_patches = get_patches(surf_cells)
    border_nodes = [ sNode(point, patch, nothing) for (point,patch) in zip(surf_points,surf_patches)]

    # find normals for border nodes
    for node in border_nodes
        node.normals = faces_normal(node.faces, facetol)
    end

    # arrays of flags
    in_border = falses(length(mesh.points))
    border_idxs = [ n.point.id for n in border_nodes ]
    in_border[border_idxs] .= true

    # map vector for border nodes
    map_pn = zeros(Int, length(mesh.points)) # map point-node
    for (i,node) in enumerate(border_nodes)
        map_pn[ node.point.id ] = i
    end

    # Stats
    Q = Float64[ c.quality for c in mesh.cells]
    q    = mean(Q)
    qmin = minimum(Q)
    qmax = maximum(Q)
    dev  = stdm(Q, q)
    q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])

    stats = DTable()
    hists = DTable()
    push!(stats, OrderedDict(:qavg=>q, :qmin=>qmin, :qmax=>qmax, :dev=>dev))

    hist  = fit(Histogram, Q, 0.0:bin:1.0, closed=:right).weights
    push!(hists, OrderedDict(Symbol(r) => v for (r,v) in zip(0.0:bin:1-bin,hist)))

    verbose && @printf("%4s  %5s  %5s  %5s  %5s  %5s  %5s  %7s  %9s  %10s\n", "it", "qmin", "q1", "q2", "q3", "qmax", "qavg", "sdev", "time", "histogram (0:$bin:1]")
    verbose && @printf("%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s", 0, qmin, q1, q2, q3, qmax, q, dev, "-")
    verbose && println("  ", str_histogram(hist))

    #nits = 0
    mesh.cell_scalar_data["quality"] = Q
    savesteps && save(mesh, "$filekey-0.vtk", verbose=false)

    for i=1:maxit

        sw = StopWatch()

        # Forces vector needed for smoothing
        D = find_disps(mesh, patches, extended, alpha, qmin, in_border)

        # Save last step file with current forces
        #savesteps && save(mesh, "$filekey-$(i-1).vtk", verbose=false)

        # Update mesh
        for point in mesh.points
            id = point.id
            X0 = [point.x, point.y, point.z][1:ndim]
            X  = X0 + vec(D[id,:])

            # skip surface nodes over non-flat locations
            if in_border[id]
                node = border_nodes[ map_pn[id] ]
                normals = node.normals
                nnorm   = length(normals)

                (nnorm == 1 && ndim==2) || (nnorm in (1,2) && ndim==3) || continue

                ΔX = X - X0
                if nnorm==1
                    n1 = normals[1]
                    X = X0 + ΔX - dot(ΔX,n1)*n1 # projection to surface plane
                else
                    n3 = normalize(cross(normals[1], normals[2]))
                    X = X0 + dot(ΔX,n3)*n3 # projection to surface edge
                end
            end

            # update key point coordinates
            point.x = X[1]
            point.y = X[2]
            if ndim==3 point.z = X[3] end

            if smart
                patch = patches[point.id]

                patch_qmin0 = minimum( c.quality for c in patch )

                # get patch new quality values
                patch_q = [ cell_quality(c) for c in patch ]

                patch_qmin = minimum(patch_q)
                #γ = 0.8
                γ = 1.0
                if patch_qmin < γ*patch_qmin0
                #if patch_qmin < patch_qmin0
                    # restore point coordinates if no improvement
                    point.x = X0[1]
                    point.y = X0[2]
                    if ndim==3; point.z = X0[3] end
                else
                    # update quality values: important
                    for (c,q) in zip(patch, patch_q)
                        c.quality = q
                    end
                end
            end

        end

        for c in mesh.cells
            c.quality = cell_quality(c)
        end

        Q = Float64[ c.quality for c in mesh.cells]
        new_q = mean(Q)
        new_qmin = minimum(Q)
        new_qmin <= 0.0 && error("smooth!: got negative quality value (qmin=$new_qmin).")

        mesh.cell_scalar_data["quality"] = Q
        savesteps && save(mesh, "$filekey-$i.vtk", verbose=false)

        Δq    = abs(q - new_q)
        Δqmin = new_qmin - qmin

        q    = new_q
        qmin = new_qmin
        qmax = maximum(Q)
        dev  = stdm(Q, q)
        q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])

        push!(stats, OrderedDict(:qavg=>q, :qmin=>qmin, :qmax=>qmax, :dev=>dev))

        hist = fit(Histogram, Q, 0.0:bin:1.0, closed=:right).weights
        push!(hists, OrderedDict(Symbol(r) => v for (r,v) in zip(0.0:bin:1-bin,hist)))

        verbose && @printf("%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s", i, qmin, q1, q2, q3, qmax, q, dev, see(sw, format=:ms))
        verbose && println("  ", str_histogram(hist))

        if Δq<tol && Δqmin<mintol && i>1
            break
        end

        #nits = i
    end
    # Set forces to zero for the last step
    #mesh.point_vector_data["forces"] = zeros(length(mesh.points), 3)
    #savesteps && save(mesh, "$filekey-$nits.vtk", verbose=false)

    savedata && save(stats, "$filekey-stats.dat")
    savedata && save(hists, "$filekey-hists.dat")
    verbose && println("  done.")

    return nothing
end



function smooth!(mesh::Mesh; verbose=true, alpha::Float64=1.0, target::Float64=0.97,
                 fixed::Bool=false, maxit::Int64=10, mintol::Float64=2e-2, tol::Float64=1e-3,
                 facetol=1e-4, savesteps::Bool=false, savedata::Bool=false, bin::Float64=0.05,
                 filekey::String="smooth", conds=nothing, extended=false, smart=false)

    # tol   : tolerance in change of mesh quality for succesive iterations
    # mintol: tolerance in change of worst cell quality in a mesh for succesive iterations

    verbose && printstyled("Mesh $(smart ? "smart-" : "")cells-fitting smoothing:\n", bold=true, color=:cyan)

    # check for not allowed cells
    for c in mesh.cells
        if c.shape.family != SOLID_SHAPE
            error("smooth!: cells of family $(c.shape.family) are not allowed for smoothing: $(c.shape.name)")
        end
    end

    # Elastic constants
    E  = 1.0
    nu = 0.0

    ndim = mesh.ndim
    npoints = length(mesh.points)
    points, patches = get_patches(mesh)  # key points and corresponding patches

    # Stats
    Q = Float64[ c.quality for c in mesh.cells]
    q    = mean(Q)
    qmin = minimum(Q)
    qmax = maximum(Q)
    dev  = stdm(Q, q)
    q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])

    stats = DTable()
    hists = DTable()
    push!(stats, OrderedDict(:qavg=>q, :qmin=>qmin, :qmax=>qmax, :dev=>dev))

    hist  = fit(Histogram, Q, 0.0:bin:1.0, closed=:right).weights
    push!(hists, OrderedDict(Symbol(r) => v for (r,v) in zip(0.0:bin:1-bin,hist)))

    verbose && @printf("%4s  %5s  %5s  %5s  %5s  %5s  %5s  %7s  %9s  %10s\n", "it", "qmin", "q1", "q2", "q3", "qmax", "qavg", "sdev", "time", "histogram (0:$bin:1]")
    verbose && @printf("%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s", 0, qmin, q1, q2, q3, qmax, q, dev, "-")
    verbose && println("  ", str_histogram(hist))


    # Lagrange multipliers matrix
    A   = mountA(mesh, fixed, conds, facetol)
    nbc = size(A,1)

    nits = 0

    for i=1:maxit

        #if i>2; extended=true end
        sw = StopWatch()

        # Forces vector needed for smoothing
        #verbose && print("calculating forces...")
        F   = force_bc(mesh, E, nu, alpha, extended)

        if ndim==2
            mesh.point_vector_data["forces"] = [ reshape(F, ndim, npoints)' zeros(npoints)]
        else
            mesh.point_vector_data["forces"] = reshape(F, ndim, npoints)'
        end

        # Save last step file with current forces
        savesteps && save(mesh, "$filekey-$(i-1).vtk", verbose=false)

        # Augmented forces vector
        F   = vcat( F, zeros(nbc) )

        # global stiffness plus LM
        #verbose && print("\rmounting stiffness matrix...")
        K = mountKg(mesh, E, nu, A)
        #verbose && println("\r")

        # Solve system
        #verbose && print("\rsolving system...         ")
        LUf = lu(K)
        U = LUf\F

        # Update mesh
        #verbose && print("\rupdating mesh...  ")
        for point in mesh.points
            id = point.id
            X0 = [point.x, point.y, point.z][1:ndim]
            pos = (point.id-1)*ndim+1
            X   = X0 + vec( U[pos:pos+ndim-1] )

            # update key point coordinates
            point.x = X[1]
            point.y = X[2]
            if ndim==3 point.z = X[3] end

            if smart
                patch = patches[point.id]

                patch_qmin0 = minimum( c.quality for c in patch )
                patch_qavg0 = mean( c.quality for c in patch )

                # get patch new quality values
                patch_q = [ cell_quality(c) for c in patch ]

                patch_qmin = minimum(patch_q)
                patch_qavg = mean(patch_q)

                #γ = 1.0
                #γ = 0.9
                γ = 0.8
                if patch_qmin < γ*patch_qmin0
                #if patch_qavg < γ*patch_qavg0
                    # restore point coordinates if no improvement
                    point.x = X0[1]
                    point.y = X0[2]
                    if ndim==3; point.z = X0[3] end
                else
                    # update quality values
                    for (c,q) in zip(patch, patch_q)
                        c.quality = q
                    end
                end
            end

        end

        if !smart
            for c in mesh.cells
                c.quality = cell_quality(c)
            end
        end

        Q = Float64[ c.quality for c in mesh.cells]
        new_q = mean(Q)
        new_qmin = minimum(Q)
        mesh.cell_scalar_data["quality"] = Q

        Δq    = abs(q - new_q)
        Δqmin = new_qmin - qmin

        q    = new_q
        qmin = new_qmin
        qmax = maximum(Q)
        dev  = stdm(Q, q)
        q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])

        push!(stats, OrderedDict(:qavg=>q, :qmin=>qmin, :qmax=>qmax, :dev=>dev))

        hist = fit(Histogram, Q, 0.0:bin:1.0, closed=:right).weights
        push!(hists, OrderedDict(Symbol(r) => v for (r,v) in zip(0.0:bin:1-bin,hist)))

        verbose && @printf("%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s", i, qmin, q1, q2, q3, qmax, q, dev, see(sw, format=:ms))
        verbose && println("  ", str_histogram(hist))

        nits = i

        if Δq<tol && Δqmin<mintol && i>1
            break
        end
    end
    @show nits

    n_bad_cells = count(q -> q<=0, Q)
    n_bad_cells>0 && @warn "Invalid cells found: " n=n_bad_cells

    # Set forces to zero for the last step
    mesh.point_vector_data["forces"] = zeros(length(mesh.points), 3)
    savesteps && save(mesh, "$filekey-$nits.vtk", verbose=false)

    savedata && save(stats, "$filekey-stats.dat")
    savedata && save(hists, "$filekey-hists.dat")

    # update data at current mesh structure
    #mesh.cell_scalar_data["quality"] = Q

    verbose && println("  done.")

    return nothing
end





"""
laplacian_smooth!(mesh; maxit, verbose, mintol, tol, savesteps, savedata, filekey, smart, weighted)

Smooths a finite element mesh using Laplacian smoothing (standard, weighted, smart).
"""
function laplacian_smooth!(mesh::Mesh; maxit::Int64=20, verbose::Bool=true, fixed=false,
    mintol::Float64=1e-2, tol::Float64=1e-4, facetol::Float64=1e-5,
    savesteps::Bool=false, savedata::Bool=false, bin::Float64=0.05,
    filekey::String="smooth", smart=false, weighted=false)

    verbose && printstyled("Mesh $(smart ? "smart-" : "")$(weighted ? "weighted-" : "")Laplacian smoothing:\n", bold=true, color=:cyan)

    ndim = mesh.ndim

    # find element patches
    points, patches = get_patches(mesh)  # key points and corresponding patches

    # list of points per patch (without key point)
    P = Array{Point,1}[]
    for (point, patch) in zip(points, patches)
        patch_points = get_points(patch)
        idx = findfirst( p -> hash(p) == hash(point), patch_points)
        splice!(patch_points, idx) # remove key point
        push!(P, patch_points)
    end

    # get a list of surface nodes (sNodes) that include a list of faces per point
    surf_cells = get_surface(mesh.cells)
    surf_points, surf_patches = get_patches(surf_cells)
    border_nodes = [ sNode(point, patch, nothing) for (point,patch) in zip(surf_points,surf_patches)]

    # find normals for border nodes
    for node in border_nodes
        node.normals = faces_normal(node.faces, facetol)
    end

    # arrays of flags
    in_border = falses(length(mesh.points))
    border_idxs = [ n.point.id for n in border_nodes ]
    in_border[border_idxs] .= true

    # map vector for border nodes
    map_pn = zeros(Int, length(mesh.points)) # map point-node
    for (i,node) in enumerate(border_nodes)
        map_pn[ node.point.id ] = i
    end

    savesteps && save(mesh, "$filekey-0.vtk", verbose=false)

    # stats
    Q    = Float64[ c.quality for c in mesh.cells]
    q    = mean(Q)
    qmin = minimum(Q)
    qmax = maximum(Q)
    dev  = stdm(Q, q)
    q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])

    stats = DTable()
    hists = DTable()
    push!(stats, OrderedDict(:qavg=>q, :qmin=>qmin, :qmax=>qmax, :dev=>dev))

    hist  = fit(Histogram, Q, 0.0:bin:1.0, closed=:right).weights
    push!(hists, OrderedDict(Symbol(r) => v for (r,v) in zip(0.0:bin:1-bin,hist)))

    #verbose && @printf("  it: %2d  qmin: %7.5f  qavg: %7.5f  dev: %7.5f", 0, qmin, q, dev)
    #verbose && @printf("  it: %2d  q-range: %5.3f…%5.3f  qavg: %5.3f  dev: %7.5f", 0, qmin, qmax, q, dev)
    #verbose && println("  hist: ", fit(Histogram, Q, 0.5:bin:1.0, closed=:right).weights)
    verbose && @printf("%4s  %5s  %5s  %5s  %5s  %5s  %5s  %7s  %9s  %10s\n", "it", "qmin", "q1", "q2", "q3", "qmax", "qavg", "sdev", "time", "histogram (0:$bin:1]")
    verbose && @printf("%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s", 0, qmin, q1, q2, q3, qmax, q, dev, "-")
    #verbose && println("  ", fit(Histogram, Q, 0.3:bin:1.0, closed=:right).weights)
    verbose && println("  ", str_histogram(hist))

    # find center point and update point position
    for i=1:maxit
        sw = StopWatch()

        for (point, patch, patch_points) in zip(points, patches, P)
            id = point.id
            X0 = [point.x, point.y, point.z][1:ndim]

            # skip surface nodes over non-flat locations
            if in_border[id]
                fixed && continue # skip if fixed boundary
                node = border_nodes[ map_pn[id] ]
                normals = node.normals
                nnorm   = length(normals)

                (nnorm == 1 && ndim==2) || (nnorm in (1,2) && ndim==3) || continue
            end

            # Weighed Laplacian smoothing
            if weighted
                Xcs = [ mean(getcoords(c),dims=1)[1:ndim] for c in patch ] # centroidal coordinates for each cell in patch
                Vs  = [ c.extent for c in patch ] # areas or volumes
                X = X0 + sum( Vs[j]*(Xcs[j]-X0) for j=1:length(patch) ) / sum(Vs)
            else # simplest Laplacian smoothing
                X  = mean(getcoords(patch_points),dims=1)[1:ndim]
            end

            # fix coordinates for surface points
            if in_border[id]
                ΔX = X - X0
                if nnorm==1
                    n1 = normals[1]
                    X = X0 + ΔX - dot(ΔX,n1)*n1 # projection to surface plane
                else
                    n3 = normalize(cross(normals[1], normals[2]))
                    X = X0 + dot(ΔX,n3)*n3 # projection to surface edge
                end
            end

            # update key point coordinates
            point.x = X[1]
            point.y = X[2]
            if ndim==3 point.z = X[3] end

            # Smart Laplacian smoothing
            if smart
                # get current minimum quality value in patch
                patch_qavg0 = mean( c.quality for c in patch )
                patch_qmin0 = minimum( c.quality for c in patch )

                # get patch new quality values
                patch_q = [ cell_quality(c) for c in patch ]
                patch_qavg = mean(patch_q)
                patch_qmin = minimum(patch_q)

                if patch_qmin < patch_qmin0
                #if patch_qavg < patch_qavg0
                    # restore point coordinates if no improvement
                    point.x = X0[1]
                    point.y = X0[2]
                    if ndim==3; point.z = X0[3] end
                else
                    # update quality values
                    for (c,q) in zip(patch, patch_q)
                        c.quality = q
                    end
                end
            end
        end

        if !smart
            for c in mesh.cells
                c.quality = cell_quality(c)
            end
        end

        # stats
        Q = Float64[ c.quality for c in mesh.cells]
        new_quality = mean(Q)
        new_qmin    = minimum(Q)

        mesh.cell_scalar_data["quality"] = Q
        savesteps && save(mesh, "$filekey-$i.vtk", verbose=false)

        #if any(Q .== 0.0)
            #error("smooth!: Smooth procedure got invalid element(s).")
        #end

        Δq    = abs(q - new_quality)
        Δqmin = new_qmin - qmin

        q    = new_quality
        qmin = new_qmin
        qmax = maximum(Q)
        dev  = stdm(Q, q)
        q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])

        push!(stats, OrderedDict(:qavg=>q, :qmin=>qmin, :qmax=>qmax, :dev=>dev))

        hist = fit(Histogram, Q, 0.0:bin:1.0, closed=:right).weights
        push!(hists, OrderedDict(Symbol(r) => v for (r,v) in zip(0.0:bin:1-bin,hist)))

        verbose && @printf("%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s", i, qmin, q1, q2, q3, qmax, q, dev, see(sw, format=:ms))
        #verbose && println("  ", fit(Histogram, Q, 0.3:bin:1.0, closed=:right).weights)
        verbose && println("  ", str_histogram(hist))

        #Δqmin<0.0 && break

        if Δq<tol && Δqmin<mintol && i>2
            break
        end
    end

    n_bad_cells = count(q -> q<=0, Q)
    n_bad_cells>0 && @warn "Invalid cells found: " n=n_bad_cells

    savedata && save(stats, "$filekey-stats.dat")
    savedata && save(hists, "$filekey-hists.dat")

    # update data at current mesh structure
    mesh.cell_scalar_data["quality"] = Q

    verbose && println("  done.")

    return nothing

end
