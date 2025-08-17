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

| Name                                                                                    | Gridap name                                  | FEEC name | Polytopes   | Order           | Conformity|
| :-------------------------------------------------------------------------------------- | :------------------------------------------- | :-------- | :---------- | --------------: | :-------- |
| [Lagrangian](https://defelement.org/elements/lagrange.html)                             | [`lagrangian`](@ref LagrangianRefFE)         | 𝓟ᵣ⁽⁻⁾Λ⁰   | △           | ``{r=o≥1, o}``  | `:H1`     |
|                                                                                         |                                              | 𝓠ᵣ⁻Λ⁰     | ``\square`` | ``{r=o≥1, o}``  | `:H1`     |
|                                                                                         |                                              |     |`WEDGE`, `PYRAMID` | ``{o≥1, o}``    | `:H1`     |
| [Serendipity](https://defelement.org/elements/serendipity.html)                         | [`serendipity`](@ref SerendipityRefFE)       | 𝓢ᵣΛ⁰      | ``\square`` | ``{r=o≥1, o}``  | `:H1`     |
| [Bezier](https://defelement.org/elements/bernstein.html)                                | [`bezier`](@ref BezierRefFE)                 | 𝓟ᵣ⁻Λ⁰     | △           | ``{r=o≥1, o}``  | `:H1`     |
|                                                                                         |                                              | 𝓠ᵣ⁻Λ⁰     | ``\square`` | ``{r=o≥1, o}``  | `:H1`     |
| [ModalC0](https://doi.org/10.48550/arXiv.2201.06632)                                    | [`modalC0`](@ref ModalC0RefFE)               | 𝓠ᵣ⁻Λ⁰     | ``\square`` | ``{r=o≥1, o}``  | `:H1`     |
|                                                                                                                                                                                                |
| [Nédélec (first kind)](https://defelement.org/elements/nedelec1.html)                   | [`nedelec`](@ref NedelecRefFE)               | 𝓟ᵣ⁻Λ¹     | `TRI`,`TET` | ``{r=o+1≥1, r}``| `:Hcurl`  |
|                                                                                         |                                              | 𝓠ᵣ⁻Λ¹     | `QUAD`,`HEX`| ``{r=o+1≥1, r}``| `:Hcurl`  |
| [Nédélec (second kind)](https://defelement.org/elements/nedelec2.html)                  | `TODO`                                       | 𝓟ᵣΛ¹      | `TRI`,`TET` | ``{r=o≥1,   r}``| `:Hcurl`  |
|                                                                                                                                                                                                |
| [Raviart-Thomas](https://defelement.org/elements/raviart-thomas.html)                   | [`raviart_thomas`](@ref LagrangianRefFE)     | 𝓟ᵣ⁻Λᴰ⁻¹   | `TRI`,`TET` | ``{r=o+1≥1, r}``| `:Hdiv`   |
|                                                                                         |                                              | 𝓠ᵣ⁻Λᴰ⁻¹   | `QUAD`,`HEX`| ``{r=o+1≥1, r}``| `:Hdiv`   |
| [Brezzi-Douglas-Marini](https://defelement.org/elements/brezzi-douglas-marini.html)     | [`bdm`](@ref BDMRefFE)                       | 𝓟ᵣΛᴰ⁻¹    | `TRI`,`TET` | ``{r=o≥1,   r}``| `:Hdiv`   |
| [Mardal-Tai-Winther](https://defelement.org/elements/mardal-tai-winther.html)           | `TODO` `mtw`                                 |           | `TRI`,`TET` | ``{o=1, D+1}``  | `:Hdiv`   |
|                                                                                                                                                                                                |
| [Crouzeix-Raviart](https://defelement.org/elements/crouzeix-raviart.html)               |[`couzeix_raviart`](@ref CrouzeixRaviartRefFE)|           |  `TRI`      | ``{o=1, o}``    | `:L2`     |
| [discontinuous Lagrangian](https://defelement.org/elements/discontinuous-lagrange.html) | [`lagrangian`](@ref LagrangianRefFE)         | 𝓟ᵣ⁻Λᴰ     | △           | ``{r-1=o≥0, o}``| `:L2`     |
|                                                                                         |                                              | 𝓟ᵣΛᴰ      | △           | ``{r=o≥0, o}``  | `:L2`     |
|                                                                                         |                                              | 𝓠ᵣ⁻Λᴰ     | ``\square`` | ``{r-1=o≥0, o}``| `:L2`     |
|                                                                                         |                                              | 𝓢ᵣΛᴰ      | ``\square`` | ``{r=o≥0, o}``  | `:L2`     |
| [MINI bubble](@ref "Bubble reference FE")                                               | [`bubble`](@ref BubbleRefFE)                 |           |△,``\square``| ``{o=1, 2}``    | `:L2`     |
| Bezier, ModalC0                                                                         | as above                                     |           |             | ``{o≥0, o}``    | `:L2`     |
|                                                                                                                                                                                                |
| [Arnold-Winther](https://defelement.org/elements/arnold-winther.html)                   | `TODO`  `arnoldwinther`                      |           | `TRI`       | ``{o=2, 4}``    | `:Hdiv`   |
| [Hellan-Herrmann-Jhonson](https://defelement.org/elements/hellan-herrmann-johnson.html) | `TODO`  `hhj`                                |           | `TRI`       | ``{TODO, o}``   | `:Hdiv`   |

###### Legend

- Name: usual name of the element and link to its
    [DefElement](https://defelement.org/) page, containing all the details
    defining the element and references.
- Gridap name: the name  to use in the Gridap APIs (it is a [`ReferenceFEName`](@ref ReferenceFEs.ReferenceFEName)
    singleton), with a link to the docstring of the element constructor.
- FEEC name: name of the element family in the [Periodic Table of the Finite
    Elements](https://www-users.cse.umn.edu/~arnold/femtable/index.html) [1]
- Polytopes:
  - △ simplices (`SEGMENT`, `TRI`  (triangle),      `TET` (tetrahedron))
  - ``\square`` n-cubes   (`SEGMENT`, `QUAD` (quadridateral), `HEX` (hexahedron))
- Order: ( _definition of ``r`` and ``o``_; _``k``_ ) where
  - ``r`` is the FEEC polynomial degree parameter (if defined).
  - ``o`` is the `order` parameter of the non FEEC `ReferenceFE` constructor (using [`name::ReferenceFEName`](@ref ReferenceFEName)),
  - ``k`` is the maximum polynomial order of the shape functions in one direction (Lagrange superdegree), defined in function of `k` or `o`,
- Conformity: supported [`Conformity`](@ref). All the elements also implement `:L2` conformity.

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

###### Bubble reference FE

The [`BubbleRefFE`](@ref) currently implements bubble space for MINI element,
but the bubble can be fine-tuned. The MINI element on `TRI`angle is documented
[here](https://defelement.org/elements/mini.html).

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

| Name | Extrusion | Nodes | Edges | Faces |
| :--- | :-------- | :---- | :---- | :---- |
| `SEGMENT` | `(HEX,)` | ![](assets/polytopes/SEGMENT_0.svg) | |  |
| `TRI` | `(TET,TET)` | ![](assets/polytopes/TRI_0.svg) | ![](assets/polytopes/TRI_1.svg) | |
| `QUAD` | `(HEX,HEX)` | ![](assets/polytopes/QUAD_0.svg) | ![](assets/polytopes/QUAD_1.svg) | |
| `HEX` | `(HEX,HEX,HEX)` | ![](assets/polytopes/HEX_0.svg) | ![](assets/polytopes/HEX_1.svg) | ![](assets/polytopes/HEX_2.svg) |
| `TET` | `(TET,TET,TET)` | ![](assets/polytopes/TET_0.svg) | ![](assets/polytopes/TET_1.svg) | ![](assets/polytopes/TET_2.svg) |
| `WEDGE` | `(TET,TET,HEX)` | ![](assets/polytopes/WEDGE_0.svg) | ![](assets/polytopes/WEDGE_1.svg) | ![](assets/polytopes/WEDGE_2.svg) |
| `PYRAMID` | `(HEX,HEX,TET)` | ![](assets/polytopes/PYRAMID_0.svg) | ![](assets/polytopes/PYRAMID_1.svg) | ![](assets/polytopes/PYRAMID_2.svg) |

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
Pages   = ["LagrangianRefFEs.jl","LagrangianDofBases.jl","SerendipityRefFEs.jl","BezierRefFEs.jl","ModalC0RefFEs.jl","BubbleRefFEs.jl"]
```

### Moment-Based ReferenceFEs

#### Framework

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["MomentBasedReferenceFEs.jl","GeometricDecompitions.jl"]
```

##### Geometric decompositions

The kwarg `sh_is_pb=true` means that the shape functions are defined as the
basis polynomials of the pre-basis, the basis must verify a geometric
decomposition. `sh_is_pb` is only available for BDM, Raviart-Thomas and
Nédélec elements, and set to true by default on simplices. Otherwise, if
`sh_is_pb=false`, the shape functions are defined as the dual basis of the DoF
basis.

The geometric decomposition API consist in the methods
- [`has_geometric_decomposition(polybasis,p,conf)`](@ref),
- [`get_face_own_funs(polybasis,p,conf)`](@ref),
- [`get_facet_flux_sign_flip(polybasis,p,conf)`](@ref).
where `polybasis` is a polynomial basis, `p` a polytope and `conf` a conformity.

This API ensures that each polynomial `𝑝ᵢ` of the basis are all associated to a
face of `p`, such that the `conf`-trace of `𝑝ᵢ` (scalar trace, tangential trace,
normal trace) over another face `g` of `p` is zero whenever `g` does not
contain `f`, and that these polynomial can be glued together with conformity
`conf` in the physical space.

Currently, only [`BernsteinBasisOnSimplex`](@ref),
[`BarycentricPmΛBasis`](@ref) and [`BarycentricPΛBasis`](@ref) implement it,
see also Bernstein basis [Geometric decomposition](@ref "Geometric
decomposition"). It could be easily generalized to bases on n-cubes.

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["GeometricDecompositions.jl"]
```

#### Available Moment-Based ReferenceFEs

```@autodocs
Modules = [ReferenceFEs,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["RaviartThomasRefFEs.jl","NedelecRefFEs.jl","BDMRefFEs.jl","CrouzeixRaviartRefFEs.jl"]
```

## References

[1] [D.N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014.](https://www-users.cse.umn.edu/~arnold/papers/periodic-table.pdf)
