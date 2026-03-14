
"""
    from_jld2_file(path::AbstractString, dataset::AbstractString="data")

Loads an object from a JLD2 file.

The `dataset` parameter specifies the name/root location of the data
inside the file.
"""
function from_jld2_file(filename::AbstractString,dataset::AbstractString="data")
    load(filename)[dataset]
end

"""
    from_jld2_file(::Type{T}, path::AbstractString, dataset::AbstractString="data") where T

Loads an object of type `T` from a JLD2 file.

Throws an error if the loaded object is not of the expected type.
"""
function from_jld2_file(::Type{T},filename::AbstractString,dataset::AbstractString="data") where T
    obj = from_jld2_file(filename,dataset)
    if isa(obj,T)
        return obj
    else 
        error("Object from JLD2 file is of type ",typeof(obj)," instead of the expected ",T)
    end
end

"""
    to_jld2_file(object, path::AbstractString, dataset::AbstractString="data")

Serializes `object` into a JLD2 file.

The `dataset` parameter specifies the name/root location of the data
inside the file.
"""
function to_jld2_file(object,filename::AbstractString,dataset::AbstractString="data")
    save(DataFormat{:JLD2}, filename, dataset, object)
end
