function from_jld2_file(filename::AbstractString,dataset::AbstractString="data")
    load(filename)[dataset]
end

function to_jld2_file(object,filename::AbstractString,dataset::AbstractString="data")
    save(DataFormat{:JLD2}, filename, dataset, object)
end
