# Gridap at a glance

## Meshes

## FESpaces

## FEOperators

## Performance mode

Internally, Gridap has many checks to ensure errors are caught early. These are provided through the [`Gridap.Helpers.@check`](@ref) macro, which by default is equivalent to Julia's `@assert`. These checks can however be completely deactivated by through [`Gridap.Helpers.set_execution_mode`](@ref). Note this will recompile the code, so you have to restart Julia after changing the preference.

```@docs
Gridap
```
