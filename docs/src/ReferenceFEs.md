```@meta
CurrentModule = Gridap.ReferenceFEs
```

# Gridap.ReferenceFEs

The following are the following finite elements implemented in Gridap. The
conformities are the default ones, but every element also implement the
[L2Conformity](@ref).

A reference finite element can be manually constructed using
```@doc
ReferenceFE(::Polytope, ::ReferenceFEName, args...; kwargs...)
```

| Name                    | Gridap name     | PTFE name  | polytopes   | orders  | conf.     | ref.|
| :---------------------- | :-------------  | :--------- | :--------   | :------ | :-------- | :-- |
| Lagrangian              | `lagrangian`    | 𝓟⁽⁻⁾Λ⁰     | △           | ``≥1``  | `:H1`     |     |
|                         |                 | 𝓠⁻Λ⁰       | ``\square`` | ``≥1``  | `:H1`     |     |
|                         |                 |       |`WEDGE`,`PYRAMID` | ``≥1``  | `:H1`     |     |
| Serendipity             | `serendipity`   | 𝓢Λ⁰        | ``\square`` | ``≥1``  | `:H1`     |     |
| Bezier                  | `bezier`        | 𝓟⁻Λ⁰       | △           | ``≥1``  | `:H1`     |     |
|                         |                 | 𝓠⁻Λ⁰       | ``\square`` | ``≥1``  | `:H1`     |     |
| ModalC0                 | `modalC0`       | 𝓠⁻Λ⁰       | ``\square`` | ``≥1``  | `:H1`     |     |
|                                                                                                  |
| Nédélec (first kind)    | `nedelec`       | 𝓟⁻Λ¹       | △           | ``≥0``  | `:Hcurl`  |     |
|                         |                 | 𝓠⁻Λ¹       | ``\square`` | ``≥0``  | `:Hcurl`  |     |
| Nédélec (second kind)   | `nedelec`       | 𝓟Λ¹        | △           | `TODO`  | `:Hcurl`  |     |
|                                                                                                  |
| Raviart-Thomas          | `raviart_thomas`| 𝓟⁻Λᴰ⁻¹     | △           | ``≥0``  | `:Hdiv`   |     |
|                         |                 | 𝓠⁻Λᴰ⁻¹     | ``\square`` | ``≥0``  | `:Hdiv`   |     |
| Brezzi-Douglas-Marini   | `bdm`           | 𝓟Λᴰ⁻¹      | △           | ``≥1``  | `:Hdiv`   |     |
| Mardal-Tai-Winther      | `mtw`           |            | △           | ``3``   |           |     |
|                                                                                                  |
| Discontinuous Lagrangian| `lagrangian`    | 𝓟⁽⁻⁾Λᴰ     | △           | ``≥0``  | `:L2`     |     |
|                         |                 | 𝓠⁻Λᴰ       | ``\square`` | ``≥0``  | `:L2`     |     |
|                         |                 |       |`WEDGE`,`PYRAMID` | ``≥0``  | `:L2`     |     |
| Crouzeix-Raviart        |`couzeix_raviart`|            |  `TRI`      | ``1``   | `:L2`     |     |
|                                                                                                  |
| Arnold-Winther          | `arnoldwinther` |            | `TRI`       | ``2``   | `:Hdivdiv`|     |
| Hellan-Herrmann-Jhonson | `hhj`           |            | `TRI`       | `TODO`  | `:Hdivdiv`|     |

###### Legend
- Gridap name: the name ([ReferenceFEName](@ref) singleton) to use in Gridap
APIs, with a link to the docstring of the element constructor.
- PTFE name: name of the element family in the [Periodic Table of the Finite
Elements](https://www-users.cse.umn.edu/~arnold/femtable/index.html) [1]
- Polytopes:
    - △ simplices (`SEGMENT`, `TRI`  (triangle),      `TET` (tetrahedron))
    - ``\square`` n-cubes   (`SEGMENT`, `QUAD` (quadridateral), `HEX` (hexahedron)
- Implemented polynomial order parameter
- Reference: link to the [DefElement](https://defelement.org/) page of the
element, containing all the details defining the element and references.

#### Contents

```@contents
Pages = ["ReferenceFEs.md"]
Depth = 2:3
```

## Polytopes

### Abstract API

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/Polytopes.jl"]
```

### Extrusion Polytopes

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["ExtrusionPolytopes.jl"]
```

### General Polytopes

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["GeneralPolytopes.jl"]
```

## Quadratures

### Abstract API

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/Quadratures.jl"]
```

### Available Quadratures

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["TensorProductQuadratures.jl","DuffyQuadratures.jl","StrangQuadratures.jl","XiaoGimbutasQuadratures.jl"]
```

## ReferenceFEs

### Abstract API

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["ReferenceFEInterfaces.jl","Dofs.jl","LinearCombinationDofVectors.jl"]
```

### Nodal ReferenceFEs

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["LagrangianRefFEs.jl","LagrangianDofBases.jl","SerendipityRefFEs.jl","BezierRefFEs.jl","ModalC0RefFEs.jl"]
```

### Moment-Based ReferenceFEs

#### Framework

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["MomentBasedReferenceFEs.jl","Pullbacks.jl"]
```

#### Available Moment-Based ReferenceFEs

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["RaviartThomasRefFEs.jl","NedelecRefFEs.jl","BDMRefFEs.jl","CrouzeixRaviartRefFEs.jl"]
```

## References

[1] [D.N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014.](https://www-users.cse.umn.edu/~arnold/papers/periodic-table.pdf)
