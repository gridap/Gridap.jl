```@meta
CurrentModule = Gridap.ReferenceFEs
```

# Gridap.ReferenceFEs

#### Contents

```@contents
Pages = ["ReferenceFEs.md"]
Depth = 2:3
```


## Reference FE summary

A reference finite element is defined using the following signature
```@docs; canonical=false
ReferenceFE(::ReferenceFEName, args...; kwargs...)
```

The following table summarizes the elements implemented in Gridap (legend below).

| Name                                                                                    | Gridap name                                  | PTFE name  | Polytopes   | Order          | Conformity|
| :-------------------------------------------------------------------------------------- | :------------------------------------------- | :--------- | :---------- | :------------- | :-------- |
| [Lagrangian](https://defelement.org/elements/lagrange.html)                             | [`lagrangian`](@ref LagrangianRefFE)         | 𝓟⁽⁻⁾Λ⁰     | △           | ``{r≥1, r}``   | `:H1`     |
|                                                                                         |                                              | 𝓠⁻Λ⁰       | ``\square`` | ``{r≥1, r}``   | `:H1`     |
|                                                                                         |                                              |      |`WEDGE`, `PYRAMID` | ``{r≥1, r}``   | `:H1`     |
| [Serendipity](https://defelement.org/elements/serendipity.html)                         | [`serendipity`](@ref SerendipityRefFE)       | 𝓢Λ⁰        | ``\square`` | ``{r≥1, r}``   | `:H1`     |
| [Bezier](https://defelement.org/elements/bernstein.html)                                | [`bezier`](@ref BezierRefFE)                 | 𝓟⁻Λ⁰       | △           | ``{r≥1, r}``   | `:H1`     |
|                                                                                         |                                              | 𝓠⁻Λ⁰       | ``\square`` | ``{r≥1, r}``   | `:H1`     |
| [ModalC0](https://doi.org/10.48550/arXiv.2201.06632)                                    | [`modalC0`](@ref ModalC0RefFE)               | 𝓠⁻Λ⁰       | ``\square`` | ``{r≥1, r}``   | `:H1`     |
|                                                                                                                                                                                                |
| [Nédélec (first kind)](https://defelement.org/elements/nedelec1.html)                   | [`nedelec`](@ref NedelecRefFE)               | 𝓟⁻Λ¹       | `TRI`,`TET` | ``{r≥1, r+1}`` | `:Hcurl`  |
|                                                                                         |                                              | 𝓠⁻Λ¹       | `QUAD`,`HEX`| ``{r≥1, r+1}`` | `:Hcurl`  |
| [Nédélec (second kind)](https://defelement.org/elements/nedelec2.html)                  | `TODO`                                       | 𝓟Λ¹        | `TRI`,`TET` | ``{r≥1, r+2}`` | `:Hcurl`  |
|                                                                                                                                                                                                |
| [Raviart-Thomas](https://defelement.org/elements/raviart-thomas.html)                   | [`raviart_thomas`](@ref LagrangianRefFE)     | 𝓟⁻Λᴰ⁻¹     | `TRI`,`TET` | ``{r≥0, r+1}`` | `:Hdiv`   |
|                                                                                         |                                              | 𝓠⁻Λᴰ⁻¹     | `QUAD`,`HEX`| ``{r≥0, r+1}`` | `:Hdiv`   |
| [Brezzi-Douglas-Marini](https://defelement.org/elements/brezzi-douglas-marini.html)     | [`bdm`](@ref BDMRefFE)                       | 𝓟Λᴰ⁻¹      | `TRI`,`TET` | ``{r≥1, r}  `` | `:Hdiv`   |
| [Mardal-Tai-Winther](https://defelement.org/elements/mardal-tai-winther.html)           | `TODO` `mtw`                                 |            | `TRI`,`TET` | ``{r=1, D+1}`` | `:Hdiv`   |
|                                                                                                                                                                                                |
| [Crouzeix-Raviart](https://defelement.org/elements/crouzeix-raviart.html)               |[`couzeix_raviart`](@ref CrouzeixRaviartRefFE)|            |  `TRI`      | ``{r=1, r}``   | `:L2`     |
| [discontinuous Lagrangian](https://defelement.org/elements/discontinuous-lagrange.html) | [`lagrangian`](@ref LagrangianRefFE)         | ...Λᴰ      | as above    | ``{r≥0, r}``   | `:L2`     |
| Serendipity, Bezier, ModalC0                                                            | as above                                     |            |             | ``{r≥0, r}``   | `:L2`     |
|                                                                                                                                                                                                |
| [Arnold-Winther](https://defelement.org/elements/arnold-winther.html)                   | `TODO`  `arnoldwinther`                      |        | `TRI`       | ``{r=2, 4}``   | `:Hdiv`   |
| [Hellan-Herrmann-Jhonson](https://defelement.org/elements/hellan-herrmann-johnson.html) | `TODO`  `hhj`                                |        | `TRI`       | ``{TODO, r}``  | `:Hdiv`   |

###### Legend

- Name: usual name of the element and link to its
    [DefElement](https://defelement.org/) page, containing all the details
    defining the element and references.
- Gridap name: the name  to use in the Gridap APIs (it is a [`ReferenceFEName`](@ref ReferenceFEs.ReferenceFEName)
    singleton), with a link to the docstring of the element constructor.
- PTFE name: name of the element family in the [Periodic Table of the Finite
    Elements](https://www-users.cse.umn.edu/~arnold/femtable/index.html) [1]
- Polytopes:
    - △ simplices (`SEGMENT`, `TRI`  (triangle),      `TET` (tetrahedron))
    - ``\square`` n-cubes   (`SEGMENT`, `QUAD` (quadridateral), `HEX` (hexahedron)
- Order: Implemented range for the polynomial order parameter, and the actual
    maximum polynomial degree of the shape functions.
- Conformity: default [`Conformity`](@ref). All the elements also implement `:L2`
    conformity (discontinuous methods).

##### Additional information

The `lagrangian`, `modalC0` and `bezier` elements support anisotropic orders on
`QUAD` and `HEX`, leveraging the tensor product basis in each dimension.

Also, a Cartesian product finite-element space is available for `lagrangian`,
`serendipity` and `modalC0` elements, for all polytopes. It means that a tensor
type ([`<:MultiValue`](@ref Gridap.TensorValues)) can be given as value type
argument `T`, for example `VectorValue{3,Float64}` or
`SymTensorValue{2,Float64}`. The DoFs are duplicated for each independent
component of the tensor.

The `modalC0` element has the particularity that it's polytope and thus the
shape function support can be adapted to the physical element.

```@docs
ReferenceFEs
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
Pages   = ["ReferenceFEInterfaces.jl","Dofs.jl","LinearCombinationDofVectors.jl","Pullbacks.jl"]
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
Pages   = ["MomentBasedReferenceFEs.jl"]
```

#### Available Moment-Based ReferenceFEs

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["RaviartThomasRefFEs.jl","NedelecRefFEs.jl","BDMRefFEs.jl","CrouzeixRaviartRefFEs.jl"]
```

## References

[1] [D.N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014.](https://www-users.cse.umn.edu/~arnold/papers/periodic-table.pdf)
