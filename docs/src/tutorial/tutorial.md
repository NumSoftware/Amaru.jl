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

A node is represented by an instance of the `Node` type.
It contains the node coordinates, a identification number and a tag
string that can be used to label a group of nodes.

```@example
using Amaru
node = Node(1.0, 2.0, 0.0, id=1, tag="vertex")
```

### Cells

A mesh cell is represented by an instace of the `Cell` type.
It is defined by a given shape object, 
a list of nodes, a identification number an a tag string.
The shape object contains information related with the geometry of
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
A `Block` type represents a geometric segment, area (4 or 8 point quadrilateral) or volume (8 or 20 point hexahedron).
An instance of the block type 
is generated from a coordinates matrix and predefined number of divisions divisions (`nx`, `ny`, `nz`) in the ``x``, ``y`` and ``z`` directions of a local system.
The block below represents a 2D `Block` with `nx=7` and `ny=5`.
The `mplot` function saves the block image to a file in svg format. Other formats as pdf and png can be used.

```@example
using Amaru
block = Block([0 0; 4 0; 3 2; 1 2], nx=7, ny=5)
mplot(block, "block.svg", markers=true)
nothing # hide
```
![](./block.svg)

Blocks can be combined to define a complex geometry.
To manipulate blocks, Amaru provides several operators as
`copy`, `move!`, `scale!`, `rotate`, `mirror`, `array`, `polar`, `extrude`, etc.

For example, let's generate a square mesh with a central hole.
We start with a quadratic quadrilateral block. Note that 
the vertex coordinates follow the same numbering as conventional
finite elements. Thus, the corner nodes are listed first, and then
the middle nodes. All the nodes are listed couterclockwise.

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
mplot(block1, "block1.svg", markers=true)
nothing # hide
```
![](./block1.svg)

A mirror operation can be applyed to this block to get a second one.
Then, both blocks can be grouped in an array for further manipulation.
```@example 2
block2 = mirror(block1, base=[0, 0], axis=[0.7643, -0.7643])
blocks = [block1, block2]
mplot(blocks, "blocks.svg", markers=true)
nothing # hide
```
![](./blocks.svg)

The center of the circle that contains the arc is located at the coordinates (1,1).
To place the center at the origin (0,0) we perform a `move!` operation.
```@example 2
move!(blocks, dx=-1.0, dy=-1.0)
mplot(blocks, "moved.svg", markers=true, axis=true)
nothing # hide
```
![](./moved.svg)

Next, we can apply a `polar` operation to the current geometry to 
obtain the square region with a central hole.
Note that the original geometry was replicated four times (`n=4`) with different rotation angles.
For this purpose a rotation axis and a base point were employed.
```@example 2
hole = polar(blocks, base=[0, 0], axis=[0, 0, 1], angle=360, n=4)
mplot(hole, "hole.svg", markers=true)
nothing # hide
```
![](./hole.svg)

Then, an `extrude` operation is used to obtain a volume following a defined axis and length.
The argument `n=4` represents the number of divisions intended along the new dimension.

```@example 2
solid = extrude(hole, axis=[0, 0, 1], length=1, n=4)
mplot(solid, "solid.svg", markers=true, dist=13)
nothing # hide
```
![](./solid.svg)

Note that, the `solid` variable is an array with four blocks, as shown in the figure. This array is
used to generate the structured mesh.
```@example 2
mesh = Mesh(solid)
mplot(mesh, "mesh.svg", dist=13)
nothing # hide
```
![](./mesh.svg)

Finally, a `smooth` operation can be applied to improve the cells quality.
```@example 2
smooth_mesh = smooth!(mesh)
mplot(smooth_mesh, "smooth_mesh.svg", dist=13)
nothing # hide
```
![](./smooth_mesh.svg)


## Mesh generation examples

Below are presented some examples of structured mesh generation using Amaru.

### Simple 2D mesh

```@example
using Amaru

block = Block([0 0 ; 2 2], nx=8, ny=6)
mesh  = Mesh(block) 
mplot(mesh, "mesh-quad4.svg", markers=true)
nothing # hide
```
![](./mesh-quad4.svg)

### Mesh with quadratic cells

```@example
using Amaru

block = Block([0 0 ; 2 2], nx=8, ny=6, cellshape=QUAD8)
mesh  = Mesh(block) 
mplot(mesh, "mesh-quad8.svg", markers=true)
nothing # hide
```
![](./mesh-quad8.svg)

###  Simple 3D mesh

```@example
using Amaru
block = Block([0 0 0; 2 4 3], nx=3, ny=6, nz=6)
mesh = Mesh(block)
mplot(mesh, "mesh-hex8.svg", markers=true)
nothing # hide
```
![](./mesh-hex8.svg)

###  2D mesh constructed from two quadratic blocks

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
mplot(mesh, "mesh-2d.svg", markers=true)
nothing # hide
```
![](./mesh-2d.svg)


###  Mesh with cells growing in the ``x`` direction

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
mplot(mesh, "mesh-rate.svg", markers=true)
nothing # hide
```
![](./mesh-rate.svg)

### 3D mesh obtained revolving the last example
```@example x
mesh = revolve(mesh, minangle=0, maxangle=90)
changeaxes!(mesh, "xzy")
rotate!(mesh, axis=[0,0,1], angle=90)
mplot(mesh, "mesh-rev.svg", azim=-45) 
nothing # hide
```
![](./mesh-rev.svg)


###  Wire mesh

```@example y
using Amaru

block1 = Block([0 0 0; 1 0 1; 0.7 0 0.3 ], n=7) # curved wire (quadratic)
block2 = Block([1 0 1; 1 0 2], n=5) # straight wire
mesh = Mesh(block1, block2)

mplot(mesh, "mesh-arc.svg", azim=-90, elev=0, markers=true)
nothing # hide
```
![](./mesh-arc.svg)

### 3D mesh obtained revolving the last example

```@example y
mesh = revolve(mesh, n=24, axis=[0,0,1]) # surface by revolution
mplot(mesh, "mesh-surf.svg", elev=15)
nothing # hide
```
![](./mesh-surf.svg)


## Finite element model

To performe a finite element analysis, an instance of the `Model` type is required.
This represents the domain model and contains nodes and finite elements.
The element instances are different from the cells in a `Mesh` object since
the elements are constructed according to the analysis problem.
To instantiate a `Model` object a mesh and a list of material models objects are needed.
The `Model` type is used for all analysis problems, e.g. mechanical, thermal, seepage, etc.
Thus, for a particular problem, a consistent list of materials should be used.
For example, in a mechanical analysis, all elements should be associated to mechanical material types.
Based on the choosen materials, the `Model` object will setup the corresponding finite elements and degrees of freedom.

### Material definitions

There are several material models implemented in Amaru and are associated to a corresponding type, e.g. 
`ElasticSolid`, `DruckerPrager`, `VonMises`, `ElasticRod`, etc.
A particular material instance is defined by calling the construcor with set of corresponding parameters.
```@example
using Amaru

steel = ElasticSolid(E=2e8, nu=0.2) # E in kPa
rock = DruckerPrager(E=2e8, nu=0.2, alpha=0.4, kappa=2500, H=1e4)
rebar = ElasticRod(E=2e8, A=1e-4)

nothing # hide
```

### Model generation

The following example shows the walk through of creating a `Model` object.
Note that the list of materials contains a pair relating a domain region 
to material instance.
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
model = Model(mesh, materials)

nothing # hide
```

### Boundary conditions

The boundary conditions are given by objects that define the quantities that should be applied
to nodes, faces, edges or even elements (`NodeBC`, `SurfaceBC`, `EdgeBC` and `ElementBC`).

The example below shows a nodal boundary condition where all displacements are set to zero.
```@example
using Amaru

fixed_end = NodeBC(ux=0, uy=0, uz=0)
```

This other example shows a face boudary condition where a traction is applied in the negative ``z`` direction.
```@example
using Amaru

load = SurfaceBC(tz=-10)
```
The keys (e.g. `ux`, `tz`, etc.) to be set in a boundary condition are related to the analysis problem, and ultimately to the material types in use.
For example, in a mechanical analysis, the keys for nodal boundary conditions are `ux`, `uy`, `uz`, `fx`, `fy` and `fz` where
`ux` represents displacement and `fx` concentrated force, both in the ``x`` direction.
Similarly, the keys for a face boudary condition are `ux`, `uy`, `uz`, `tx`, `ty` and `tz`, where `tx` is a traction (force per area)
in the ``x`` direction. Besides, the `tn` key can be used to apply traction normal to the face.


## Analysis example

This example shows a mechanical analysis in a steel beam.
The left end was clamped, thus all displacements at `x==0` are set to zero.
Also, a vertical traction is applied at the right end.
The `solve!` function performs the calculations of a `Model` instance
subjected to a set of boundary conditions.
The function also updates the state variables at the elements' integration points.

For post-processing in a visualization software (e.g. Paraview), the 
model can save the current per node and per element values (e.g. displacements, stresses, etc.) in "vtk" and "vtu" formats. 
To save the full state of a domain model, the output in "xml" format is also supported; thus, the model can be reused in future analyses.

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
model = Model(mesh, materials)

bcs = [
    :(x==0) => NodeBC(ux=0, uy=0, uz=0)
    :(x==1) => SurfaceBC(tz=-1000)
]

solve!(model, bcs, verbose=true)
save(model, "beam.vtu")

nothing # hide
```

The `mplot` function allow to save a plot of the domain in several formats (e.g. "pdf", "svg", "png", etc.).
In this case, the deformed state is plot using a mangnification scale of 100.
Also the ``u_y`` field is displayed and the corresponding label is placed next to the colorbar.
To adjunst the point of view, the `azim`, `elev` and `dist` parameters can be set. They represent, respectively, the azimut, elevation angle and distance to the rendered domain.

```@example fem
mplot(model, "beam.svg", warpscale=100, field=:uz, colorbarlabel=raw"$u_z$", colorbarscale=0.4, azim=-90, elev=30, dist=7)
nothing # hide
```
![](./beam.svg)

 
## Nonlinear analysis


The use of material models with nonlinear behavior lead to a nonlinear analyses.
For example, the `VonMises` material type represents a nonlinear material model.
For the analysis, the number of increments can be adjusted using the `nincs` parameter.
Also, the `autoinc` parameter can be set to true to le the `solve!` function to 
automatically compute and recompute the increment sizes.
The `tol` paramter specifies the maximum force discrepancy between internal and external forces.

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
model = Model(mesh, materials)

bcs = [
    :(x==0) => NodeBC(ux=0, uy=0, uz=0)
    :(x==1) => SurfaceBC(uz=-0.04)
]

solve!(model, bcs, tol=0.01, autoinc=true)
nothing # hide
```

```@example fem
mplot(model, "beam-vm.svg", warpscale=10, field=:sxx, fieldmult=1e-3, colorbarlabel=raw"$\sigma_{xx}$", colorbarscale=0.4, azim=-90, dist=7)
nothing # hide
```
![](./beam-vm.svg)


### Loggers

Loggers are objects designed to track the information of nodes, integration points, faces, edges, points and segments.
That information is updated at every increment and can be saved to a file.

```@example
using Amaru

nodelog = NodeLogger("nodelog.table")
iplog = IpLogger("nodelog.table")
facelog = FacesSumLogger("nodelog.table")
nodeslog = NodeGroupLogger("nodes.book")
nothing # hide
```

The loggers need to be associated to the entity they track by using a filter expression, 
coordinates or a tag string. 
Finally, the loggers have to be linked to a `Model` object by using the `setloggers!` function.

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
model = Model(mesh, materials)

bcs = [
    :(x==0) => NodeBC(ux=0, uy=0, uz=0)
    :(x==1) => SurfaceBC(uz=-0.04)
]

loggers = [
    :(x==1) => FacesSumLogger("tip.table"),
    [0.05, 0.05, 0.1] => IpLogger("ip.table"),
    :(y==0.1 && z==0.1) => NodeGroupLogger("topgroup.book"),
    [0 0.05 0.1; 1 0.05 0.1] => SegmentLogger("top.book", n=40),
    [0 0.05 0.0; 1 0.05 0.0] => SegmentLogger("bottom.book", n=40),
]

setloggers!(model, loggers)

solve!(model, bcs, tol=0.01, nouts=2, autoinc=true)

nothing # hide
```

The loggers are classified single loggers and group loggers.
Single loggers record a table of data. Group loggers record a book, which is a set of tables.
Once the analysis is finished, the loggers data can be recovered in another script 
calling the constructors of the types `DataTable` and `DataBook`.

```@example
using Amaru

tip = DataTable("tip.table")
```

```@example
using Amaru

long = DataBook("topgroup.book").tables[end]
```

## Plotting

Amaru provides a chart ploting function `cplot`
to ease the visualization of the data stored in `DataTable` and `DataBook` objects.

```@example plot
using Amaru

tip = DataTable("tip.table")
cplot(
    [
        (x=tip.uz, y=tip.fz, marker="o"),
    ],
    filename = "uz_vs_fz.svg",
    xmult  = 1e3,
    xlabel = raw"$u_z$ [mm]",
    ylabel = raw"$f_z$ [kN]",
)

topgroup = DataBook("topgroup.book").tables[end]
cplot(
    [
        (x=topgroup.x, y=topgroup.uz, marker="o"),
    ],
    filename = "x_vs_uz.svg",
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
    filename = "x_vs_sxx.svg",
    ymult  = 1e-3,
    xlabel = raw"$x$ [m]",
    ylabel = raw"$\sigma_{xx}$ [MPa]",
)

nothing # hide
```

![](./uz_vs_fz.svg)

![](./x_vs_uz.svg)

![](./x_vs_sxx.svg)
