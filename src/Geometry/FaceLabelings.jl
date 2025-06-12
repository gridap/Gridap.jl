
"""
    struct FaceLabeling <: GridapType
      d_to_dface_to_entity::Vector{Vector{Int32}}
      tag_to_entities::Vector{Vector{Int32}}
      tag_to_name::Vector{String}
    end
"""
struct FaceLabeling <: GridapType
  d_to_dface_to_entity::Vector{Vector{Int32}}
  tag_to_entities::Vector{Vector{Int32}}
  tag_to_name::Vector{String}
end

"""
    FaceLabeling(d_to_num_dfaces::Vector{Int})
    FaceLabeling(topo::GridTopology)
    FaceLabeling(topo::GridTopology,cell_to_tag::Vector{<:Integer},tag_to_name::Vector{String})
"""
function FaceLabeling(d_to_num_dfaces::Vector{Int})
  d_to_dface_to_entity = [
    fill(Int32(UNSET),n) for n in d_to_num_dfaces]
  tag_to_entities = Vector{Int32}[]
  tag_to_name = String[]
  FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

function FaceLabeling(topo::GridTopology)
  D = num_cell_dims(topo)
  d_to_ndfaces = [ num_faces(topo,d) for d in 0:D ]
  labels = FaceLabeling(d_to_ndfaces)
  for d in 0:D
    dface_to_entity = get_face_entity(labels,d)
    if d != D
      dface_to_is_boundary = get_isboundary_face(topo,d)
      dface_to_entity .= dface_to_is_boundary
      dface_to_entity .+= 1 
    else
      dface_to_entity .= 1 
    end
  end
  add_tag!(labels,"interior",[1])
  add_tag!(labels,"boundary",[2])
  labels
end

function FaceLabeling(topo::GridTopology,cell_to_tag::Vector{<:Integer},tag_to_name::Vector{String})
  a = FaceLabeling(topo)
  b = face_labeling_from_cell_tags(topo,cell_to_tag,tag_to_name)
  merge!(a,b)
  return a
end

"""
"""
function FaceLabeling(labeling::FaceLabeling,D::Integer)
  if D == num_dims(labeling)
    return labeling
  else
    d_to_dface_to_entity = labeling.d_to_dface_to_entity[1:(D+1)]
    tag_to_entities = labeling.tag_to_entities
    tag_to_name = labeling.tag_to_name
    return FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
  end
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
"""
function get_cell_entity(lab::FaceLabeling)
  get_face_entity(lab,num_dims(lab))
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
  @unreachable "Unknown tag name $name"
  0
end

"""
    get_tags_from_names(lab::FaceLabeling,names::Vector{String})
"""
function get_tags_from_names(lab::FaceLabeling,names::Vector{String})
  f = (s) -> get_tag_from_name(lab,s)
  map(f, names)
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
    add_tag_from_tags!(lab::FaceLabeling, name::String, tag::Int)
    add_tag_from_tags!(lab::FaceLabeling, name::String, tag::String)

Adds a new tag `name`, given by the union of the tags `tags`.
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

function add_tag_from_tags!(labels::FaceLabeling, name::String, tag::Int)
  tags = [tag, ]
  add_tag_from_tags!(labels,name,tags)
end

function add_tag_from_tags!(labels::FaceLabeling, name::String, tag::String)
  tags = [tag, ]
  add_tag_from_tags!(labels,name,tags)
end

"""
    add_tag_from_tags_intersection!(labels::FaceLabeling, name::String, tags::Vector{Int})
    add_tag_from_tags_intersection!(labels::FaceLabeling, name::String, tags::Vector{String})

Adds a new tag `name`, given by the intersection of the tags `tags`.
"""
function add_tag_from_tags_intersection!(labels::FaceLabeling, name::String, tags::Vector{Int})
  f(entity_tags) = issubset(tags,entity_tags)
  add_tag_from_tag_filter!(labels,name,f)
end

function add_tag_from_tags_intersection!(labels::FaceLabeling, name::String, names::Vector{String})
  tags = [get_tag_from_name(labels,name) for name in names ]
  add_tag_from_tags_intersection!(labels,name,tags)
end

"""
    add_tag_from_tags_complementary!(lab::FaceLabeling, name::String, tags::Vector{Int})
    add_tag_from_tags_complementary!(lab::FaceLabeling, name::String, tags::Vector{String})
    add_tag_from_tags_complementary!(lab::FaceLabeling, name::String, tag::Int)
    add_tag_from_tags_complementary!(lab::FaceLabeling, name::String, tag::String)

Adds a new tag `name`, given by the complementary of the tags `tags`.
"""
function add_tag_from_tags_complementary!(labels::FaceLabeling, name::String, tags::Vector{Int})
  f(entity_tags) = isdisjoint(tags,entity_tags)
  add_tag_from_tag_filter!(labels,name,f)
end

function add_tag_from_tags_complementary!(labels::FaceLabeling, name::String, names::Vector{String})
  tags = [get_tag_from_name(labels,name) for name in names ]
  add_tag_from_tags_complementary!(labels,name,tags)
end

function add_tag_from_tags_complementary!(labels::FaceLabeling, name::String, tag::Int)
  tags = [tag,]
  add_tag_from_tags_complementary!(labels,name,tags)
end

function add_tag_from_tags_complementary!(labels::FaceLabeling, name::String, tag::String)
  tags = [tag,]
  add_tag_from_tags_complementary!(labels,name,tags)
end

"""
    add_tag_from_tags_setdiff!(lab::FaceLabeling, name::String, tags_include::Vector{Int}, tags_exclude::Vector{Int})
    add_tag_from_tags_setdiff!(lab::FaceLabeling, name::String, tags_include::Vector{String}, tags_exclude::Vector{String})

Adds a new tag `name`, given by the set difference between `tags_include` and `tags_exclude`.
"""
function add_tag_from_tags_setdiff!(
  labels::FaceLabeling, name::String, tags_include::Vector{Int}, tags_exclude::Vector{Int}
)
  f(entity_tags) = issubset(tags_include,entity_tags) && isdisjoint(tags_exclude,entity_tags)
  add_tag_from_tag_filter!(labels,name,f)
end

function add_tag_from_tags_setdiff!(
  labels::FaceLabeling, name::String, names_include::Vector{String}, names_exclude::Vector{String}
)
  tags_include = [get_tag_from_name(labels,name) for name in names_include ]
  tags_exclude = [get_tag_from_name(labels,name) for name in names_exclude ]
  add_tag_from_tags_setdiff!(labels,name,tags_include,tags_exclude)
end

"""
    add_tag_from_tag_filter!(lab::FaceLabeling, name::String, filter::Function)

Adds a new tag `name`, by including all entities selected by a filter function. 
The filter function must have signature 

    filter(entity_tags::Vector{Int}) -> Bool

where `entity_tags` are the tags of a particular geometrical entity.
"""
function add_tag_from_tag_filter!(
  labels::FaceLabeling, name::String, filter::Function
)
  tag_to_entities = labels.tag_to_entities
  n_entities = maximum(e -> maximum(e,init=0), tag_to_entities)
  entity_to_tags = [Int32[] for i in 1:n_entities]
  for (tag,entities) in enumerate(tag_to_entities)
    for e in entities
      push!(entity_to_tags[e],tag)
    end
  end

  entities = collect(Int32,findall(filter,entity_to_tags))
  add_tag!(labels,name,entities)
end

"""
    get_face_mask(labeling::FaceLabeling,tags::Vector{Int},d::Integer)
    get_face_mask(labeling::FaceLabeling,tags::Vector{String},d::Integer)
    get_face_mask(labeling::FaceLabeling,tag::Int,d::Integer)
    get_face_mask(labeling::FaceLabeling,tag::String,d::Integer)
"""
function get_face_mask(labeling::FaceLabeling,tags,d::Integer)
  _tags = _prepare_tags(labeling,tags)
  face_to_entity = get_face_entity(labeling,d)
  tag_to_entities = get_tag_entities(labeling)
  nfaces = num_faces(labeling,d)
  face_to_mask = fill(false,nfaces)
  for (face, entity) in enumerate(face_to_entity)
    for tag in _tags
      entities = tag_to_entities[tag]
      if entity in entities
        face_to_mask[face] = true
      end
    end
  end
  face_to_mask
end

"""
    get_face_tag(labeling::FaceLabeling,tags::Vector{Int},d::Integer)
    get_face_tag(labeling::FaceLabeling,tags::Vector{String},d::Integer)
    get_face_tag(labeling::FaceLabeling,tag::Int,d::Integer)
    get_face_tag(labeling::FaceLabeling,tag::String,d::Integer)
    get_face_tag(labeling::FaceLabeling,d::Integer)

The first of the given tags appearing in the face is taken.
If there is no tag on a face, this face will have a value equal to `UNSET`.
If not tag or tags are provided, all the tags in the model are considered
"""
function get_face_tag(labeling::FaceLabeling,tags,d::Integer)
  _tags = _prepare_tags(labeling,tags)
  face_to_entity = get_face_entity(labeling,d)
  tag_to_entities = get_tag_entities(labeling)
  nfaces = num_faces(labeling,d)
  face_to_tag = fill(Int8(UNSET),nfaces)
  for (face, entity) in enumerate(face_to_entity)
    for tag in _tags
      entities = tag_to_entities[tag]
      if entity in entities
        face_to_tag[face] = tag
        break
      end
    end
  end
  face_to_tag
end

function get_face_tag(labeling::FaceLabeling,d::Integer)
  tags = collect(1:num_tags(labeling))
  get_face_tag(labeling,tags,d)
end

"""
    get_face_tag_index(labeling::FaceLabeling,tags::Vector{Int},d::Integer)
    get_face_tag_index(labeling::FaceLabeling,tags::Vector{String},d::Integer)
    get_face_tag_index(labeling::FaceLabeling,tag::Int,d::Integer)
    get_face_tag_index(labeling::FaceLabeling,tag::String,d::Integer)

Like `get_face_tag` by provides the index into the array `tags` instead of the tag stored in `tags`.
"""
function get_face_tag_index(labeling::FaceLabeling,tags,d::Integer)
  _tags = _prepare_tags(labeling,tags)
  face_to_entity = get_face_entity(labeling,d)
  tag_to_entities = get_tag_entities(labeling)
  nfaces = num_faces(labeling,d)
  face_to_tag = fill(Int8(UNSET),nfaces)
  for (face, entity) in enumerate(face_to_entity)
    for (i,tag) in enumerate(_tags)
      entities = tag_to_entities[tag]
      if entity in entities
        face_to_tag[face] = i
        break
      end
    end
  end
  face_to_tag
end

"""
    restrict(labels::FaceLabeling,d_to_dface_to_parent_dface)
"""
function restrict(labels::FaceLabeling,d_to_dface_to_parent_dface)
  D = length(d_to_dface_to_parent_dface)-1

  d_to_dface_to_entity = [
    zeros(Int32,length(d_to_dface_to_parent_dface[d+1])) for d in 0:D
  ]

  for d in 0:D
    face_to_label = d_to_dface_to_entity[d+1]
    oface_to_label = labels.d_to_dface_to_entity[d+1]
    face_to_oface = d_to_dface_to_parent_dface[d+1]
    for face in eachindex(face_to_label)
      oface = face_to_oface[face]
      label = oface_to_label[oface]
      face_to_label[face] = label
    end
  end

  tag_to_entities = copy(labels.tag_to_entities)
  tag_to_name     = copy(labels.tag_to_name)
  return FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

"""
    Base.merge!(a::FaceLabeling,b::FaceLabeling...)
"""
function Base.merge!(a::FaceLabeling,b::FaceLabeling...)
  for b_ in b
    merge!(a,b_)
  end
  return a
end

function Base.merge!(a::FaceLabeling,b::FaceLabeling)
  msg = "merge! not implemented if FaceLabeling objects have common tags"
  @notimplementedif !isempty(intersect(a.tag_to_name,b.tag_to_name)) msg
  @assert length(a.d_to_dface_to_entity) == length(b.d_to_dface_to_entity)

  D = length(a.d_to_dface_to_entity) - 1
  offset = maximum(maximum,a.tag_to_entities;init=0)
  append!(a.tag_to_name, b.tag_to_name)
  append!(a.tag_to_entities, (e .+ offset for e in b.tag_to_entities))

  n_entities = offset + maximum(maximum,b.tag_to_entities;init=0)
  entity_to_tags = [Int32[] for i in 1:n_entities]
  for (tag,entities) in enumerate(a.tag_to_entities)
    for e in entities
      push!(entity_to_tags[e],tag)
    end
  end

  entities = Dict{Tuple,Int32}()
  for d in 0:D
    a_dface_to_entity = a.d_to_dface_to_entity[d+1]
    b_dface_to_entity = b.d_to_dface_to_entity[d+1]
    @assert length(a_dface_to_entity) == length(b_dface_to_entity)
    for dface in eachindex(a_dface_to_entity)
      a_unset = isequal(a_dface_to_entity[dface],UNSET)
      b_unset = isequal(b_dface_to_entity[dface],UNSET)
      if a_unset && !b_unset
        # If a is UNSET but not b, take b
        a_dface_to_entity[dface] = b_dface_to_entity[dface] + offset
      elseif !a_unset && !b_unset
        # If both are set, we consider the combined entity
        dface_entities = (a_dface_to_entity[dface],b_dface_to_entity[dface]+offset)

        if !haskey(entities,dface_entities)
          # If the combined entity is new, we add it
          n_entities += 1
          entities[dface_entities] = n_entities
          for tag in vcat((entity_to_tags[e] for e in dface_entities)...)
            push!(a.tag_to_entities[tag],n_entities)
          end
        end

        a_dface_to_entity[dface] = entities[dface_entities]
      end
      # If a is set but not b, or if both are unset, we do nothing
    end
  end

  return a
end

"""
    face_labeling_from_cell_tags(topo::GridTopology,cell_to_tag::Vector{<:Integer},tag_to_name::Vector{String})

Creates a FaceLabeling object from a GridTopology and a vector containing a tag for each cell.

You can use this function to add custom cell-defined tags to a pre-existing FaceLabeling, by 
calling `merge!` on the two FaceLabeling objects, i.e 

```julia
  labels = FaceLabeling(topo)
  additional_labels = face_labeling_from_cell_tags(topo,cell_to_tag,tag_to_name)
  merge!(labels,additional_labels)
```

"""
function face_labeling_from_cell_tags(topo::GridTopology, cell_to_tag, tag_to_name)
  D = num_cell_dims(topo)
  n_tags = length(tag_to_name)
  tag_to_entities = [Int32[] for tag in 1:n_tags]
  d_to_dface_to_entity = [fill(UNSET,num_faces(topo,d)) for d in 0:D]

  # Cell entities: 1:n_tags
  d_to_dface_to_entity[D+1] .= cell_to_tag
  for tag in 1:n_tags
    push!(tag_to_entities[tag],tag)
  end

  # Face entities:
  n_entities = n_tags
  entities = Dict{Set{Int},Int32}([Set{Int}(i) => i for i in 1:n_tags])
  for d in D-1:-1:0
    dface_to_cells = get_faces(topo,d,D)
    dface_to_entity = d_to_dface_to_entity[d+1]
    for (dface,cells) in enumerate(dface_to_cells)
      dface_tags = Set(cell_to_tag[cells])
      if haskey(entities,dface_tags)
        dface_entity = entities[dface_tags]
      else
        n_entities += 1
        dface_entity = n_entities
        entities[dface_tags] = dface_entity
        for tag in dface_tags
          push!(tag_to_entities[tag],dface_entity)
        end
      end
      dface_to_entity[dface] = dface_entity
    end
  end

  foreach(unique!,tag_to_entities)
  return FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

"""
    face_labeling_from_vertex_filter(topo::GridTopology, name::String, filter::Function)

Creates a FaceLabeling object from a GridTopology and a coordinate-based filter function. 
The filter function takes in vertex coordinates and returns a boolean values. A geometrical 
entity is tagged if all its vertices pass the filter.

You can use this function to add custom vertex-defined tags to a pre-existing FaceLabeling, by
calling `merge!` on the two FaceLabeling objects, i.e

```julia
  labels = FaceLabeling(topo)
  additional_labels = face_labeling_from_vertex_filter(topo,name,filter)
  merge!(labels,additional_labels)
```

"""
function face_labeling_from_vertex_filter(topo::GridTopology, name, filter)
  D = num_cell_dims(topo)
  d_to_dface_to_entity = [fill(UNSET,num_faces(topo,d)) for d in 0:D]

  vertex_mask = map(filter,get_vertex_coordinates(topo))

  for d in 0:D
    dface_to_vertices = get_faces(topo,d,0)
    dface_to_entity = d_to_dface_to_entity[d+1]
    for (dface,vertices) in enumerate(dface_to_vertices)
      if all(vertex_mask[vertices])
        dface_to_entity[dface] = 1
      end
    end
  end

  tag_to_entities = [Int32[1]]
  tag_to_name = [name]
  return FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

function _prepare_tags(labeling,tags::Vector{<:Integer})
  tags
end

function _prepare_tags(labeling,names::Vector{String})
  get_tags_from_names(labeling,names)
end

function _prepare_tags(labeling,tag::Integer)
  [tag,]
end

function _prepare_tags(labeling,name::String)
  get_tags_from_names(labeling,[name,])
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

function to_dict(labels::FaceLabeling)
  dict = Dict{Symbol,Any}()
  D = num_dims(labels)
  dict[:D] = D
  for d in 0:D
    k = Symbol("entities_$d")
    dict[k] = get_face_entity(labels,d)
  end
  dict[:tags] = labels.tag_to_entities
  dict[:names] = labels.tag_to_name
  dict
end

function from_dict(::Type{FaceLabeling},dict::Dict{Symbol,Any})
  D::Int = dict[:D]
  d_to_dface_to_entity = Vector{Vector{Int32}}(undef,D+1)
  for d in 0:D
    k = Symbol("entities_$d")
    dface_to_entity::Vector{Int32} = dict[k]
    d_to_dface_to_entity[d+1] = dface_to_entity
  end
  tag_to_entities::Vector{Vector{Int32}} = dict[:tags]
  tag_to_name::Vector{String} = dict[:names]
  FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end


