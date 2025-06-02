```@meta
CurrentModule = Gridap.ReferenceFEs
```

# Gridap.ReferenceFEs

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
