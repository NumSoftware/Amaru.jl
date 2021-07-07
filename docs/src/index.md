```@meta
DocTestSetup = quote
    using Amaru
end
```

# Introduction to Amaru

Amaru is a Finite Element library written in Julia language. The purpose of this library is to aid the research of new algorithms for the finite element method. Currently this library solves static and dynamic analysis in two and three-dimensions.

## Installation and basic usage

Install the package using the package manager:

```
] add Amaru
```

To use Amaru, type:

```
using Amaru
```

To test:

```
] test Amaru
```