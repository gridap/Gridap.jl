
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

function compute_face_to_cell_map(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_face_to_cell_map(trian::BoundaryTriangulation)
  key = :get_face_to_cell_map
  memo = get_memo(trian)
  if ! haskey(memo,key)
    memo[key] = compute_face_to_cell_map(trian)
  end
  memo[key]
end

"""
"""
function get_normal_vector(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_face_to_face(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_cell_around(trian::BoundaryTriangulation)
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
  @test isa(get_face_to_cell_map(trian),CellField)
  @test isa(get_normal_vector(trian),CellField)
  @test get_face_to_cell_map(trian) === get_face_to_cell_map(trian)
end

# Default API

function compute_cell_map(trian::BoundaryTriangulation)
  vtrian = get_volume_triangulation(trian)
  ϕ = get_cell_map(vtrian)
  ϕ_Γ = _get_cell_map(trian)
  FaceMap(ϕ_Γ,ϕ,get_face_to_cell(trian),get_face_to_cell_map(trian))
end

function CellField(object,trian::BoundaryTriangulation)
  CellField(object,get_volume_triangulation(trian))
end

function get_cell_id(trian::BoundaryTriangulation)
  get_face_to_cell(trian)
end

# Constructors

"""
    BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
    BoundaryTriangulation(model::DiscreteModel)
"""
function BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool},icell_around)
  GenericBoundaryTriangulation(model,face_to_mask,icell_around)
end

function BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
  icell_around = 1
  BoundaryTriangulation(model,face_to_mask,icell_around)
end

function BoundaryTriangulation(model::DiscreteModel)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  face_to_mask = collect(Bool,get_isboundary_face(topo,D-1))
  BoundaryTriangulation(model,face_to_mask)
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

