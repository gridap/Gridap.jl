
"""
    struct FaceLabeling
      d_to_dface_to_entity::Vector{Vector{Int32}}
      tag_to_entities::Vector{Vector{Int32}}
      tag_to_name::Vector{String}
    end
"""
struct FaceLabeling
  d_to_dface_to_entity::Vector{Vector{Int32}}
  tag_to_entities::Vector{Vector{Int32}}
  tag_to_name::Vector{String}
end

"""
    FaceLabeling(d_to_num_dfaces::Vector{Int})
"""
function FaceLabeling(d_to_num_dfaces::Vector{Int})
  d_to_dface_to_entity = [
    fill(Int32(UNSET),n) for n in d_to_num_dfaces]
  tag_to_entities = Vector{Int32}[]
  tag_to_name = String[]
  FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

"""
    num_dims(lab::FaceLabeling)
"""
function num_dims(lab::FaceLabeling)
  length(lab.d_to_dface_to_entity)-1
end

"""
    num_cell_dims(lab::FaceLabeling)
"""
function num_cell_dims(lab::FaceLabeling)
  num_dims(lab)
end

"""
    num_tags(lab::FaceLabeling)
"""
num_tags(lab::FaceLabeling) = length(lab.tag_to_name)

"""
    num_entities(lab::FaceLabeling)
"""
function num_entities(lab::FaceLabeling)
  maximum(map(maximum,lab.d_to_dface_to_entity))
end

"""
    num_faces(lab::FaceLabeling,d::Integer)
"""
function num_faces(lab::FaceLabeling,d::Integer)
  length(lab.d_to_dface_to_entity[d+1])
end

"""
    num_faces(lab::FaceLabeling)
"""
function num_faces(lab::FaceLabeling)
  sum(map(length,lab.d_to_dface_to_entity))
end

"""
    num_vertices(lab::FaceLabeling)
"""
num_vertices(lab::FaceLabeling) = num_faces(lab,0)

"""
    num_edges(lab::FaceLabeling)
"""
num_edges(lab::FaceLabeling) = num_faces(lab,1)

"""
    num_facets(lab::FaceLabeling)
"""
function num_facets(lab::FaceLabeling)
  D = num_dims(lab)
  if D > 0
    num_faces(lab,num_dims(lab)-1)
  else
    0
  end
end

"""
    num_cells(lab::FaceLabeling)
"""
num_cells(lab::FaceLabeling) = num_faces(lab,num_dims(lab))

"""
    get_face_entity(lab::FaceLabeling,d::Integer)
"""
function get_face_entity(lab::FaceLabeling,d::Integer)
  lab.d_to_dface_to_entity[d+1]
end

"""
    get_face_entity(lab::FaceLabeling)
"""
function get_face_entity(lab::FaceLabeling)
  # TODO return a view
  vcat(lab.d_to_dface_to_entity...)
end

"""
    get_tag_entities(lab::FaceLabeling,tag::Integer)
    get_tag_entities(lab::FaceLabeling,tag::String)
"""
function get_tag_entities(lab::FaceLabeling,tag::Integer)
  lab.tag_to_entities[tag]
end

function get_tag_entities(lab::FaceLabeling,tag::String)
  i = get_tag_from_name(lab,tag)
  get_tag_entities(lab,i)
end

"""
    get_tag_entities(lab::FaceLabeling)
"""
function get_tag_entities(lab::FaceLabeling)
  lab.tag_to_entities
end

"""
    get_tag_name(lab::FaceLabeling,tag::Integer)
"""
function get_tag_name(lab::FaceLabeling,tag::Integer)
  lab.tag_to_name[tag]
end

"""
    get_tag_name(lab::FaceLabeling)
"""
function get_tag_name(lab::FaceLabeling)
  lab.tag_to_name
end

"""
    get_tag_from_name(lab::FaceLabeling,name::String)
"""
function get_tag_from_name(lab::FaceLabeling,name::String)
  for tag in 1:num_tags(lab)
    if lab.tag_to_name[tag] == name
      return tag
    end
  end
  @unreachable
  0
end

"""
    get_tag_from_name(lab::FaceLabeling)
"""
function get_tag_from_name(lab::FaceLabeling)
  dict = Dict{String,Int}()
  for (tag, name) in enumerate(lab.tag_to_name)
    dict[name] = tag
  end
  dict
end

"""
    add_tag!(lab::FaceLabeling,name::String,entities::Vector{<:Integer})
"""
function add_tag!(lab::FaceLabeling,name::String,entities::Vector{<:Integer})
  @assert ! (name in lab.tag_to_name) "Tag name $(name) already present in this FaceLabeling object"
  push!(lab.tag_to_entities,entities)
  push!(lab.tag_to_name,name)
end

"""
    add_tag_from_tags!(lab::FaceLabeling, name::String, tags::Vector{Int})
    add_tag_from_tags!(lab::FaceLabeling, name::String, tags::Vector{String})
"""
function add_tag_from_tags!(lab::FaceLabeling, name::String, tags::Vector{Int})
  entities = Int32[]
  tag_to_entities = lab.tag_to_entities
  for tag in tags
     for entity in tag_to_entities[tag]
       push!(entities,entity)
     end
  end
  add_tag!(lab,name,sort(collect(Set(entities))))
end

function add_tag_from_tags!(
  labels::FaceLabeling, name::String, names::Vector{String})
  tags = [get_tag_from_name(labels,name) for name in names ]
  add_tag_from_tags!(labels,name,tags)
end

function Base.show(io::IO,lab::FaceLabeling)
  print(io,typeof(lab))
end

function Base.show(io::IO,::MIME"text/plain",labels::FaceLabeling)
  show(io,labels)
  print(io,":")
  for d = 0:num_dims(labels)
    print(io,"\n $d-faces: $(num_faces(labels,d))")
  end
  print(io,"\n tags: $(num_tags(labels))")
  print(io,"\n entities: $(num_entities(labels))")
end

