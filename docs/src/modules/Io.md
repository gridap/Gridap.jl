
```@meta
CurrentModule = Gridap.Io
```

# Gridap.Io

The `Gridap.Io` module provides a unified interface for the serialization and de-serialization of Gridap objects (such as models, fields, and finite element spaces) into various formats, including JSON, BSON, and JLD2.

Most serialization processes in Gridap follow a two-step approach:
1.  Converting the object into a plain Julia `Dict{Symbol,Any}` using `to_dict`.
2.  Serializing this dictionary into the target format (e.g., JSON or BSON).

De-serialization follows the reverse process.

## Dictionary Interface

The core of the serialization system is the `Dict` interface.

```@docs
to_dict
from_dict
check_dict
```

## JSON

Functions for converting objects to and from JSON (JavaScript Object Notation) strings and files.

```@docs
to_json
from_json
to_json_file
from_json_file
```

## BSON

Functions for binary JSON (BSON) serialization, which can be more efficient for large datasets.

```@docs
to_bson_file
from_bson_file
```

## JLD2

Functions for JLD2 (Julia Data Format) serialization, a HDF5-compatible format for Julia.

```@docs
to_jld2_file
from_jld2_file
```
