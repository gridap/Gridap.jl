
"""
    abstract type BoundaryTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}
"""
abstract type BoundaryTriangulation{Dc,Dp} <: Triangulation{Dc,Dp} end

"""
    get_volume_triangulation(trian::BoundaryTriangulation)
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
function get_face_to_lface(trian::BoundaryTriangulation)
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

"""
    test_boundary_triangulation(trian::BoundaryTriangulation)
"""
function test_boundary_triangulation(trian::BoundaryTriangulation)
  test_triangulation(trian)
  @test isa(get_volume_triangulation(trian),Triangulation)
  @test isa(get_face_to_cell(trian),AbstractArray{<:Integer})
  @test isa(get_face_to_lface(trian),AbstractArray{<:Integer})
  @test isa(get_face_to_cell_map(trian),AbstractArray{<:Field})
  @test isa(get_normal_vector(trian),CellField)
end

# Default API

function restrict(f::AbstractArray, trian::BoundaryTriangulation)
  compose_field_arrays(reindex(f,trian), get_face_to_cell_map(trian))
end

function get_cell_id(trian::BoundaryTriangulation)
  get_face_to_cell(trian)
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
function BoundaryTriangulation(model::DiscreteModel,tags)
  labeling = get_face_labeling(model)
  BoundaryTriangulation(model,labeling,tags)
end

"""
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling,tags::Vector{Int})
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling,tags::Vector{String})
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling,tag::Int)
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling,tag::String)
"""
function BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling,tags)
  D = num_cell_dims(model)
  face_to_mask = get_face_mask(labeling,tags,D-1)
  BoundaryTriangulation(model,face_to_mask)
end

#function BoundaryTriangulation(model::DiscreteModel,names::Vector{String})
#  labeling = get_face_labeling(model)
#  tags = get_tags_from_names(labeling,names)
#  BoundaryTriangulation(model,tags)
#end
#
#function BoundaryTriangulation(model::DiscreteModel,tag::Union{Int,String})
#  tags = [tag,]
#  BoundaryTriangulation(model,tags)
#end
#
#function _convert_to_face_to_masks(labeling,tags)
#  D = num_cell_dims(model)
#  face_to_mask = get_face_mask(labeling,tags,D-1)
#end



