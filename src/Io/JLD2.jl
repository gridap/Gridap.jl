
"""
    function from_jld2_file(filename::AbstractString,dataset::AbstractString="data")

Loads an object from a JLD2 file given its filename and, optionally, the dataset name.
The dataset specifies the name and the root location of the data inside the generated JLD2 file.
"""
function from_jld2_file(filename::AbstractString,dataset::AbstractString="data")
    load(filename)[dataset]
end

"""
    function from_jld2_file(::Type{T},filename::AbstractString,dataset::AbstractString="data") where T

Loads an object from a JLD2 file given its filename and, optionally, the dataset name.
The dataset specifies the name and the root location of the data inside the generated JLD2 file.
Checks if the returned object is of the expected `Type{T}`, if not return error.
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
    function to_jld2_file(object,filename::AbstractString,dataset::AbstractString="data")

Stores an object to a JLD2 file given its filename and, optionally, the dataset name.
The dataset specifies the name and the root location of the data inside the JLD2 file.
"""
function to_jld2_file(object,filename::AbstractString,dataset::AbstractString="data")
    save(DataFormat{:JLD2}, filename, dataset, object)
end
