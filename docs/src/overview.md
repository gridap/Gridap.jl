# Gridap at a glance

Gridap provides a comprehensive suite of tools for grid-based approximation of partial differential equations (PDEs) in Julia. This page gives you a quick overview of the main capabilities and where to find detailed documentation for each feature.

## Finite Element Discretizations

Gridap supports a wide range of finite element discretizations for solving PDEs:

- **Linear and nonlinear PDE systems** - Handle both simple and complex mathematical problems
- **Transient problems** - Solve time-dependent equations
- **Multi-variable problems** - Couple different physics or solve systems of equations
- **Scalar and tensor variables** - Work with single or multi-component quantities
- **Conforming and nonconforming elements** - Choose from a variety of elements, with arbitrary polynomial order
- **Support for constraints** - Constraints through Lagrange multipliers or linear constraints

📖 **Learn more**: [Gridap.FESpaces](@ref), [Gridap.MultiField](@ref), [Gridap.ReferenceFEs](@ref), [Gridap.ODEs](@ref)

## Mesh Support

Work with flexible mesh representations:

- **Structured and unstructured meshes** - Regular grids or complex geometries
- **Polytopal meshes** - Triangular/tetrahedral or quadrilateral/hexahedral elements
- **Physical entity labeling** - Tag boundaries and regions for boundary conditions

📖 **Learn more**: [Gridap.Geometry](@ref)

## Expandable building blocks

We provide a modular architecture that allows you seamlessly build your own low-level features:

- **Tensor algebra** - Efficient tensor operations
- **Polynomials** - Comprehensive polynomial basis library
- **Fields** - Arbitrary operations on local fields
- **Quadratures** - Library of numerical quadratures for integration

📖 **Learn more**: [Gridap.Polynomials](@ref), [Gridap.TensorValues](@ref), [Gridap.Fields](@ref)

## Input/Output and Visualization

Get data in and results out:

- **VTK support** - Export your results in VTK format for visualization
- **JDL2** - Save and load Gridap objects

📖 **Learn more**: [Gridap.Io](@ref), [Gridap.Visualization](@ref)

## Extended Ecosystem

Gridap's capabilities are greatly extended through companion packages!

📖 **Learn more**: [Gridap ecosystem](@ref)

## Performance mode

Internally, Gridap has many checks to ensure errors are caught early. These are provided through the [`Gridap.Helpers.@check`](@ref) macro, which by default is equivalent to Julia's `@assert`. These checks can however be completely deactivated by through [`Gridap.Helpers.set_execution_mode`](@ref). Note this will recompile the code, so you have to restart Julia after changing the preference.

```@docs
Gridap
```
