# Example of a truss analysis

```@example
using Amaru

# 2D Truss
coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [[1, 2], [1, 5], [2, 3], [2, 6], [2, 5], [2, 4], [3, 6], [3, 5], [4, 5], [5, 6]]

# Generate a simple mesh based on coordinates and connectivities
mesh = Mesh(coord, conn)

# List of element types and material models
mats = [
        :all => MechBar => LinearElastic => (E=6.894757e7, A=0.043)
       ]

# A mechanical analysis context
ctx = MechContext(ndim=2)

# A finite element model object
model = FEModel(mesh, mats, ctx)

# A finite element analysis object
ana = MechAnalysis(model)

# List of boundary conditions
bcs = [
       :(x==0 && y==0) => NodeBC(ux=0, uy=0), # nodal restrain
       :(x==0 && y==9) => NodeBC(ux=0, uy=0),
       :(x==9 && y==0) => NodeBC(fy=-450.),  # nodal force
       :(x==18&& y==0) => NodeBC(fy=-450.),
       :(x>=0)         => BodyC(wy=-10.0)  # self weight
      ]

# Adds a load stage
addstage!(ana, bcs)

# Run the analysis
run!(ana)
```
