
"""
    function from_jld2_file(filename::AbstractString,dataset::AbstractString="data")

Loads an object from a JLD2 file given its filename and, optionally, the dataset name.
"""
function from_jld2_file(filename::AbstractString,dataset::AbstractString="data")
    load(filename)[dataset]
end

"""
    function to_jld2_file(object,filename::AbstractString,dataset::AbstractString="data")

Stores an object to a JLD2 file given its filename and, optionally, the dataset name.
"""
function to_jld2_file(object,filename::AbstractString,dataset::AbstractString="data")
    save(DataFormat{:JLD2}, filename, dataset, object)
end
