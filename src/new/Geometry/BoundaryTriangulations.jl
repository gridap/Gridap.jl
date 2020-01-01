
"""
    abstract type BoundaryTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}
"""
abstract type BoundaryTriangulation{Dc,Dp} <: Triangulation{Dc,Dp} end

"""
"""
function get_volume_triangulation(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_face_to_cell(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_face_to_cell_map(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_normal_vector(trian::BoundaryTriangulation)
  @abstractmethod
end

# Tester

function test_boundary_triangulation(trian::BoundaryTriangulation)
  test_triangulation(trian)
  @test isa(get_volume_triangulation(trian),Triangulation)
  @test isa(get_face_to_cell(trian),AbstractArray{<:Integer})
  @test isa(get_face_to_cell_map(trian),AbstractArray{<:Field})
  @test isa(get_normal_vector(trian),AbstractArray)
end

# Default API

function reindex(a::AbstractArray,trian::BoundaryTriangulation)
  face_to_cell = get_face_to_cell(trian)
  reindex(a,face_to_cell)
end

function restrict(f::AbstractArray, trian::BoundaryTriangulation)
  compose_field_arrays(reindex(f,trian), get_face_to_cell_map(trian))
end

# Constructors

"""
    BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
    BoundaryTriangulation(model::DiscreteModel)
"""
function BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
  GenericBoundaryTriangulation(model,face_to_mask)
end

function BoundaryTriangulation(model::DiscreteModel)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  face_to_mask = collect(Bool,get_isboundary_face(topo,D-1))
  GenericBoundaryTriangulation(model,face_to_mask)
end

"""
    BoundaryTriangulation(model::DiscreteModel,tags::Vector{Int})
    BoundaryTriangulation(model::DiscreteModel,tags::Vector{String})
    BoundaryTriangulation(model::DiscreteModel,tag::Int)
    BoundaryTriangulation(model::DiscreteModel,tag::String)
"""
function BoundaryTriangulation(model::DiscreteModel,tags::Vector{Int})
  labeling = get_face_labeling(model)
  D = num_cell_dims(model)
  face_to_mask = get_face_mask(labeling,tags,D-1)
  BoundaryTriangulation(model,face_to_mask)
end

function BoundaryTriangulation(model::DiscreteModel,names::Vector{String})
  labeling = get_face_labeling(model)
  tags = get_tags_from_names(labeling,names)
  BoundaryTriangulation(model,tags)
end

function BoundaryTriangulation(model::DiscreteModel,tag::Union{Int,String})
  tags = [tag,]
  BoundaryTriangulation(model,tags)
end



