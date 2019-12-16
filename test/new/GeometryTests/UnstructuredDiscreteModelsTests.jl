module UnstructuredDiscreteModelsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: ConformingTrianMock
using Gridap.Geometry: DiscreteModelMock

grid = ConformingTrianMock()

model = UnstructuredDiscreteModel(grid)
test_discrete_model(model)

m = DiscreteModelMock()

@test num_faces(model,0) == num_faces(m,0)
@test num_faces(model,1) == num_faces(m,1)
@test num_faces(model,2) == num_faces(m,2)
@test num_nodes(model) == num_nodes(m)
@test get_reffes(ReferenceFE{1},model) == get_reffes(ReferenceFE{1},m)
@test get_reffes(ReferenceFE{0},model) == get_reffes(ReferenceFE{0},m)
@test get_reffes(ReferenceFE{2},model) == get_reffes(ReferenceFE{2},m)
@test get_polytopes(Polytope{1},model) == get_polytopes(Polytope{1},m)
@test get_polytopes(Polytope{0},model) == get_polytopes(Polytope{0},m)
@test get_polytopes(Polytope{2},model) == get_polytopes(Polytope{2},m)
@test get_face_reffe_type(model,1) == get_face_reffe_type(m,1)
@test get_face_polytope_type(model,1) == get_face_polytope_type(m,1)
@test get_isboundary_face(model,0) == get_isboundary_face(m,0)
@test get_isboundary_face(model,1) == get_isboundary_face(m,1)
@test get_isboundary_face(model,2) == get_isboundary_face(m,2)
@test get_face_nodes(model,1) == get_face_nodes(m,1)
@test get_face_nodes(model,0) == get_face_nodes(m,0)
@test get_face_nodes(model,2) == get_face_nodes(m,2)
@test get_face_own_nodes(model,1) == get_face_own_nodes(m,1)
@test get_face_own_nodes(model,0) == get_face_own_nodes(m,0)
@test get_face_own_nodes(model,2) == get_face_own_nodes(m,2)
@test get_isboundary_node(model) == get_isboundary_node(m)
@test get_faces(model,2,0) == get_faces(m,2,0)
@test get_faces(model,2,1) == get_faces(m,2,1)
@test get_faces(model,1,2) == get_faces(m,1,2)
@test get_faces(model,2,1) == get_faces(m,2,1)
@test get_faces(model,0,1) == get_faces(m,0,1)
@test get_faces(model,1,0) == get_faces(m,1,0)
@test get_faces(model,1,1) == get_faces(m,1,1)
@test get_node_coordinates(model) == get_node_coordinates(m)
@test get_vertex_coordinates(model) == get_vertex_coordinates(m)           

labels = get_face_labeling(model)
@test num_entities(labels) == 0
@test num_tags(labels) == 0
@test get_face_labeling(model) === labels
@test get_face_entity(get_face_labeling(model),0) === get_face_entity(labels,0)

domain = (0,1,0,1,0,1)
partition = (2,2,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)

edge_to_isboundary = get_isboundary_face(model,1)
@test length(edge_to_isboundary) == num_edges(model)

edge_to_faces = get_faces(model,1,2)
@test length(edge_to_faces) == num_edges(model)
@test maximum(edge_to_faces.data) == num_facets(model)
face_to_edges = get_faces(model,2,1)
@test length(face_to_edges) == num_facets(model)
@test maximum(face_to_edges.data) == num_edges(model)

test_discrete_model(model)

domain = (0,1,0,1,0,1)
partition = (2,3,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)
test_discrete_model(model)
@test is_oriented(model) == true

domain = (0,1,0,1)
partition = (2,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)
order = 2

reffes = [ LagrangianRefFE(Float64,get_polytope(reffe),order) for reffe in get_reffes(model)]
model2 = replace_reffes(model,reffes)
test_discrete_model(model2)

m = DiscreteModelMock()
model = UnstructuredDiscreteModel(m)
test_discrete_model(model)

struct CellDofsNonOriented{T,P,C<:AbstractVector,F<:AbstractVector} <:AbstractVector{Vector{T}}
  cell_to_faces::Table{T,P}
  cell_to_lface_to_pindex::Table{T,P}
  cell_to_ctype::C
  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_to_num_dofs::Vector{Int}
  face_to_own_dofs::Table{T,P}
  face_to_ftype::F
  ftype_to_pindex_to_pdofs::Vector{Vector{Vector{Int}}}
end

Base.size(a::CellDofsNonOriented) = size(a.cell_to_faces)

Base.IndexStyle(::Type{<:CellDofsNonOriented}) = IndexStyle(Table)

function array_cache(a::CellDofsNonOriented)
  n_dofs = testitem(a.ctype_to_num_dofs)
  T = eltype(eltype(a))
  v = zeros(T,n_dofs)
  CachedArray(v)
end

function getindex!(cache,a::CellDofsNonOriented,cell::Integer)
  ctype = a.cell_to_ctype[cell]
  n_dofs = a.ctype_to_num_dofs[ctype]
  setsize!(cache,(n_dofs,))
  dofs = cache.array
  lface_to_own_ldofs = a.ctype_to_lface_to_own_ldofs[ctype]
  p = a.cell_to_faces.ptrs[cell]-1
  for (lface, own_ldofs) in enumerate(lface_to_own_ldofs)
    face = a.cell_to_faces.data[p+lface]
    ftype = a.face_to_ftype[face]
    pindex = a.cell_to_lface_to_pindex.data[p+lface]
    pdofs = ftype_to_pindex_to_pdofs[ftype][pindex]
    q = a.face_to_own_dofs.ptrs[face]-1
    for (i,ldof) in enumerate(own_ldofs)
      j = pdofs[i]
      dof = a.face_to_own_dofs.data[q+j]
      dofs[ldof] = dof
    end
  end
  dofs
end

function Base.getindex(a::CellDofsNonOriented,cell::Integer)
  cache = array_cache(a)
  getindex!(cache,a,cell)
end


end # module
