```@meta
CurrentModule = Amaru
DocTestSetup = quote
    using Amaru
end
```

```@setup 1
using Amaru
```

# Tutorial

## Mesh generation

A mesh is mainly composed by cells and nodes. in Amaru, a mesh is also composed by surface
cells and surface edges and are computed at the time of mesh generation.

### Nodes

A node is represented by a `Node` type.
It contains the node coordinates, and identification number and a tag
string that can be used to label a group of nodes.

```@example
using Amaru
node = Node(1.0, 2.0, 0.0, id=1, tag="vertex")
```

### Cells

A mesh cell is represented by a `Cell` type.
It is defined by a given shape, 
a list of nodes, a identification number an a tag.
The shape contains information related with the geometry of
a finite element.
For example, a `TRI3` shape represents the shape of a linear triangular 
finite element.

```@example
using Amaru;
node1 = Node(0.0, 0.0, id=1)
node2 = Node(1.0, 0.0, id=2)
node3 = Node(1.0, 1.0, id=3)
nodes = [ node1, node2, node3 ]

Cell(TRI3, nodes, id=1, tag="triangle")
```

There are many shapes defined in Amaru, for instance: 
`LIN2`, `LIN3`, `TRI3`, `TRI6`, `QUAD4`, `QUAD8`, `QUAD9` `TET4`, `TET10`, `PYR5`, `WED6`, `WED15`, `HEX8`, `HEX20`, `HEX27`, etc.

### Blocks

Blocks are entities used to aid the generation of structured meshes in 1D, 2D and 3D.
A `Block` structure represents a geometric segment, area (4 or 8 point quadrilateral) or volume (8 or 20 point hexahedron).
It is generated from a coordinates matrix and predefined number of divisions divisions (`nx`, `ny`, `nz`) in the ``x``, ``y`` and ``z`` directions of a local system.
The block below, represents a 2D `Block` with `nx=7` and `ny=5`.

```@example
using Amaru
block = Block([0 0; 4 0; 3 2; 1 2], nx=7, ny=5)
mplot(block, "block.png", markers=true)
```
![](./block.png)

Blocks can be combined to define a complex geometry.
To manipulate blocks, Amaru provides several operators as
`copy`, `move!`, `scale!`, `rotate`, `mirror`, `array`, `polar`, `extrude`, etc.

For example, let's start with a quadratic quadrilateral. Note that 
the vertex coordinates follow the same numbering as conventional
finite elements. Thus, the corner nodes are listed first; then, 
the middle nodes. All of them listed couterclockwise.

```@example 2
using Amaru

block1 = Block(
    [ 
        0.0      0.0    
        0.7643   0.7643 
        0.6667   1.0    
        0.0      1.0    
        0.38215  0.38215
        0.692    0.8724 
        0.3333   1.0    
        0.0      0.5 
    ],
    cellshape=QUAD8, nx=5, ny=5 
)
mplot(block1, "block1.png", markers=true)
```
![](./block1.png)

A mirror operation can be applyed to this block to get a second one.
Then, both blocks can be grouped in an array for further manipulation.
```@example 2
block2 = mirror(block1, base=[0, 0], axis=[0.7643, -0.7643])
blocks = [block1, block2]
mplot(blocks, "blocks.png", markers=true)
nothing # hide
```
![](./blocks.png)

To get the centre of the circle that is represented by the arc at the origin, 
we perform a `move!` operation.
```@example 2
move!(blocks, dx=-1.0, dy=-1.0)
mplot(blocks, "moved.png", markers=true, axis=true)
nothing # hide
```
![](./moved.png)

Now we can apply a `polar` operation to the current geometry to 
obtain a square region with a central hole.
The original geometry is replicated four times using a rotation axis that 
passes througth a base point.
```@example 2
hole = polar(blocks, base=[0, 0], axis=[0, 0, 1], angle=360, n=4)
mplot(hole, "hole.png", markers=true)
nothing # hide
```
![](./hole.png)

Then, the `extrude` operation is used to obtain a volume following a defined axis and length.
The number 4 here represents the number of divisions inteded along the new dimension.

```@example 2
solid = extrude(hole, axis=[0, 0, 1], length=1, n=4)
mplot(solid, "solid.png", markers=true, dist=13)
nothing # hide
```
![](./solid.png)

Note that, the `solid` variable contains a list of blocks as shown in the figure. This list is
used to generate the structured mesh.
```@example 2
mesh = Mesh(solid)
mplot(mesh, "mesh.png", dist=13)
nothing # hide
```
![](./mesh.png)

Finally, a `smooth` operation can be applied to improve the cells quality.
```@example 2
smooth_mesh = smooth!(mesh)
mplot(smooth_mesh, "smooth_mesh.png", dist=13)
nothing # hide
```
![](./smooth_mesh.png)



## Mesh generation examples

Below are presented some examples of structured mesh generation using Amaru.

Simple 2D mesh:
```@example
using Amaru

block = Block([0 0 ; 2 2], nx=8, ny=6)
mesh  = Mesh(block) 
mplot(mesh, "mesh-quad4.png", markers=true)
nothing # hide
```
![](./mesh-quad4.png)

Mesh with quadratic cells:
```@example
using Amaru

block = Block([0 0 ; 2 2], nx=8, ny=6, cellshape=QUAD8)
mesh  = Mesh(block) 
mplot(mesh, "mesh-quad8.png", markers=true)
nothing # hide
```
![](./mesh-quad8.png)

Simple 3D mesh:
```@example
using Amaru
block = Block([0 0 0; 2 4 3], nx=3, ny=6, nz=6)
mesh = Mesh(block)
mplot(mesh, "mesh-hex8.png", markers=true)
nothing # hide
```
![](./mesh-hex8.png)

2D mesh constructed from two quadratic blocks:
```@example
using Amaru
block1 = Block(
    [ 
        0.0      0.0    
        0.7643   0.7643 
        0.6667   1.0    
        0.0      1.0    
        0.38215  0.38215
        0.692    0.8724 
        0.3333   1.0    
        0.0      0.5 
    ],
    cellshape=QUAD8, nx=5, ny=5
)
block2 = mirror(block1, base=[0, 0], axis=[0.7643, -0.7643])
mesh = Mesh(block1, block2)
mplot(mesh, "mesh-2d.png", markers=true)
nothing # hide
```
![](./mesh-2d.png)


Mesh with cells growing in the ``x`` direction.
```@example x
using Amaru

R = 5.0
r = 1.0

block = BlockGrid(
    [ 0.0, r, R ], 
    [ 0.0, R ], 
    nx=[ 4, 8 ], 
    ny=[ 8 ], 
    rx=[ 1, 1.2] # elements growing rate in the x direction
)

mesh = Mesh(block)
mplot(mesh, "mesh-rate.png", markers=true)
nothing # hide
```
![](./mesh-rate.png)

3D mesh obtained revolving the last example:
```@example x
mesh = revolve(mesh, minangle=0, maxangle=90)
changeaxes!(mesh, "xzy")
rotate!(mesh, axis=[0,0,1], angle=90)
mplot(mesh, "mesh-rev.png", azim=-45) 
nothing # hide
```
![](./mesh-rev.png)


Wire mesh:
```@example y
using Amaru

block1 = Block([0 0 0; 1 0 1; 0.7 0 0.3 ], n=7) # curved wire (quadratic)
block2 = Block([1 0 1; 1 0 2], n=5) # straight wire
mesh = Mesh(block1, block2)

mplot(mesh, "mesh-arc.png", azim=-90, elev=0, markers=true)
nothing # hide
```
![](./mesh-arc.png)

3D mesh obtained revolving the last example:
```@example y
mesh = revolve(mesh, n=24, axis=[0,0,1]) # surface by revolution
mplot(mesh, "mesh-surf.png", elev=15)
nothing # hide
```
![](./mesh-surf.png)


## Finite element model

To performe a finite element analysis, a domain model object is required.
In Amaru, this object is called `Domain`. 
A mesh and a list of material models objects are required to construct a `Domain` object.
In Amaru, the `Domain` type is used for all analysis types, e.g. mechanical, thermal, seepage, etc.
For an specific analysis type, a consistent list of materials should be used.
For example, in a mechanical analysis, only material types related with a mechanical analysis can be used.
Also, based on the choosen materials, the `Domain` object will setup the corresponding finite elements and degrees of freedom.

### Material definitions

There are several material models implemented in Amaru, e.g. 
`ElasticSolid`, `DruckerPrager`, `VonMises`, `ElasticRod`, etc.
A material object is defined by its name and a set of corresponding parameters.
```@example
using Amaru

steel = ElasticSolid(E=2e8, nu=0.2) # E in kPa
```

### Model generation

The following example shows the walk through of creating a `Domain` object.
Note that the list of materials contains a pair with the domain region 
and the corresponding material.
```@example
using Amaru

# Mesh generation
block = Block([0 0 0; 1 0.1 0.1], nx=20, ny=2, nz=2, cellshape=HEX20)
mesh = Mesh(block)

# Analysis model
steel = ElasticSolid(E=2e8, nu=0.2) # E in kPa

materials = [
    :(x>=0) => steel,
]
model = Domain(mesh, materials)

nothing # hide
```


### Boundary conditions
The boundary conditions are given by objects that define the quantities that should be applied
to nodes, faces, edges or even elements (`NodeBC`, `FaceBC`, `EdgeBC` and `ElementBC`).

The example below shows a nodal boundary condition where all displacements are set to zero.
```@example
using Amaru

fixed_end = NodeBC(ux=0, uy=0, uz=0)
```

This other example shows a face boudary conditin where a traction is applied in the negative ``z`` direction.
```@example
using Amaru

load = FaceBC(tz=-10)
```
The keys to be set in a boundary condition are related to the analysis type, and ultimately to the material types in use.
For example, in a mechanical analysis, the keys for nodal boundary conditions are `ux`, `uy`, `uz`, `fx`, `fy` and `fz` where
`ux` represents displacement and `fz` concentrated force, both in the ``x`` direction.
Similarly, the keys for a face boudary condition are `ux`, `uy`, `uz`, `tx`, `ty` and `tz`, where `tx` is a traction (force per area)
in the ``x`` direction.


## Analysis example

This example shows a mechanical analysis in a steel beam.
The left end was clamped, thus all displacements at `x==0` are set to zero.
Also, a vertical traction is applied at the right end.
The `solve!` function performs the calculations according to the given boundary conditions and the domain properties.
The function also updates the state variables at the elements' integration points.
For post-processing in a visualization software (e.g. Paraview), the 
model can be saved in "vtk" and "vtu" formats. 
These output formats contain values related to nodes and elements, for example, nodal displacements, stresses, etc. corresponding to
a time-step.
To save the current state of a domain model, the output in "xml" format is also supported; thus, the model can be reused in future analyses.

```@example fem
using Amaru

# Mesh generation
block = Block([0 0 0; 1 0.1 0.1], nx=20, ny=2, nz=2, cellshape=HEX8)
mesh = Mesh(block)

# Analysis model
steel = ElasticSolid(E=200e6, nu=0.2) # E in kPa

materials = [
    :(x>=0) => steel,
]
model = Domain(mesh, materials)

bcs = [
    :(x==0) => NodeBC(ux=0, uy=0, uz=0)
    :(x==1) => FaceBC(tz=-1000)
]

solve!(model, bcs, verbosity=2)
save(model, "beam.vtu")
nothing # hide
```

```@example fem
mplot(model, "beam.png", warpscale=100, field=:uz, colorbarlabel=raw"u_z", colorbarscale=0.4, azim=-90, dist=7)
nothing # hide
```
![](./beam.png)

 
## Nonlinear analysis

```@example fem
using Amaru

# Mesh generation
block = Block([0 0 0; 1 0.1 0.1], nx=20, ny=2, nz=3, cellshape=HEX8)
mesh = Mesh(block)

# Analysis model
steel1 = ElasticSolid(E=2e8, nu=0.2) # E in kPa
steel2 = VonMises(E=2e8, nu=0.2, fy=500e3, H=1e4)

materials = [
    :(x>=0.5) => steel1,
    :(x<=0.5) => steel2
]
model = Domain(mesh, materials)

bcs = [
    :(x==0) => NodeBC(ux=0, uy=0, uz=0)
    :(x==1) => FaceBC(uz=-0.04)
]

solve!(model, bcs, autoinc=true)
nothing # hide
```

```@example fem
mplot(model, "beam-vm.png", warpscale=10, field=:sxx, fieldmult=1e-3, colorbarlabel=raw"$\sigma_{xx}$", colorbarscale=0.4, azim=-90, dist=7)
nothing # hide
```
![](./beam-vm.png)


### Loggers

```@example
using Amaru

nodelog = NodeLogger("nodelog.table")
iplog = IpLogger("nodelog.table")
facelog = FacesSumLogger("nodelog.table")
nodeslog = NodeGroupLogger("nodes.book")
nothing # hide
```


```@example loggers
using Amaru

# Mesh generation
block = Block([0 0 0; 1 0.1 0.1], nx=20, ny=2, nz=3, cellshape=HEX8)
mesh = Mesh(block)

# Analysis model
steel1 = ElasticSolid(E=2e8, nu=0.2) # E in kPa
steel2 = VonMises(E=2e8, nu=0.2, fy=500e3, H=1e4)

materials = [
    :(x>=0.5) => steel1,
    :(x<=0.5) => steel2
]
model = Domain(mesh, materials)

bcs = [
    :(x==0) => NodeBC(ux=0, uy=0, uz=0)
    :(x==1) => FaceBC(uz=-0.04)
]

loggers = [
    :(x==1) => FacesSumLogger("tip.table"),
    [0.05, 0.05, 0.1] => IpLogger("ip.table"),
    :(y==0.1 && z==0.1) => NodeGroupLogger("topgroup.book"),
    [0 0.05 0.1; 1 0.05 0.1] => SegmentLogger("top.book", n=20),
    [0 0.05 0.0; 1 0.05 0.0] => SegmentLogger("bottom.book", n=20),
]

setloggers!(model, loggers)

solve!(model, bcs, tol=0.01, nouts=2, autoinc=true)

nothing # hide
```

```@example
using Amaru

tip = DataTable("tip.table")
```

```@example
using Amaru

long = DataBook("topgroup.book").tables[end]
```


## Plotting

```@example plot
using Amaru

tip = DataTable("tip.table")
cplot(
    [
        (x=tip.uz, y=tip.fz, marker="o"),
    ],
    filename = "uz_vs_fz.png",
    xmult  = 1e3,
    xlabel = raw"$u_z$ [mm]",
    ylabel = raw"$f_z$ [kN]",
)

topgroup = DataBook("topgroup.book").tables[end]
cplot(
    [
        (x=topgroup.x, y=topgroup.uz, marker="o"),
    ],
    filename = "x_vs_uz.png",
    ymult  = 1e3,
    xlabel = raw"$x$ [m]",
    ylabel = raw"$u_z$ [mm]",
)

top = DataBook("top.book").tables[end]
bottom = DataBook("bottom.book").tables[end]
cplot(
    [
        (x=top.x, y=top.sxx, marker="o", label="top"),
        (x=bottom.x, y=bottom.sxx, marker="o", label="bottom"),
    ],
    filename = "x_vs_sxx.png",
    ymult  = 1e-3,
    xlabel = raw"$x$ [m]",
    ylabel = raw"$\sigma_{xx}$ [MPa]",
)

nothing # hide
```

![](./uz_vs_fz.png)

![](./x_vs_uz.png)

![](./x_vs_sxx.png)
