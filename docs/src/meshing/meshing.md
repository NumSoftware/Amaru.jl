# Mesh generation

## Node

### Node struct

```@docs
Node
```

### Node constructors

```@docs
Node()
Node(::Real,::Real,::Real)
Node(::AbstractArray)
```

### Node functions

```@docs
copy(::Node)
tag!(::Node, ::String)
tag!(::Array{Node,1}, ::String)
add_dof
```

## Cell

### Cell struct

```@docs
Cell
```

### Cell constructors

```@docs
Cell(::CellShape, ::Array{Node,1})
```

### Cell functions

```@docs
copy(::Cell)
tag!(::Cell, ::String)
tag!(::Array{Cell,1}, ::String)
```

## Blocks

### Block struct
```@docs
Block
```

### Block constructors
```@docs
Block(::Array{Real})
```


### Block functions
```@docs
copy(::Block)
tag!(::Block, ::String)
tag!(::Array{Block,1}, ::String)
array(::Block)
mirror(::Block)
polar(::Block)
rotate!(::Block)
scale!(::Block)
extrude(::Block)
```

## Mesh

### Mesh struct

```@docs
Mesh
```

### Mesh constructors

```@docs
Mesh(::Array{Real}, ::Array{Array{Int64,1},1}, ::Array{CellShape,1})
Mesh(::Block)
```

### Mesh functions