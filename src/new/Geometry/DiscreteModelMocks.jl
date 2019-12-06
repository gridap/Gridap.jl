
#
#  7 ---9--- 8 ---12-- 9
#
#  |         |         |
#  |         |         |
#  |         |         |
#  10   4    11   5    13
#  |         |         |
#  |         |         |
#  |         |         |
#
#  4 ---2--- 5 ---8--- 6
#
#  |         | \   3   |
#  |         |  \      |
#  |         |   \     |
#  3    1    4    6    7
#  |         |     \   |
#  |         |  2   \  |
#  |         |       \ |
#
#  1 ---1--- 2 ---5--- 3
#

struct DiscreteModelMock <: DiscreteModel{2,2} end

function get_faces(g::DiscreteModelMock,dimfrom::Integer,dimto::Integer)
  if dimfrom == 0
    if dimto == 0
      nvertex = 9
      return identity_table(Int,Int32,nvertex)
    elseif dimto == 1
      vertex_to_faces_data = [1,3,1,4,5,5,6,7,2,3,10,2,4,6,8,11,7,8,13,9,10,9,11,12,12,13]
      vertex_to_faces_ptrs = Int32[1,3,6,9,12,17,20,22,25,27]
      return Table(vertex_to_faces_data,vertex_to_faces_ptrs)
    elseif dimto == 2
      vertex_to_cells_data = [1,1,2,2,3,1,4,1,2,3,4,5,3,5,4,4,5,5]
      vertex_to_cells_ptrs = Int32[1,2,4,6,8,13,15,16,18,19]
      return Table(vertex_to_cells_data,vertex_to_cells_ptrs)
    else
      @unreachable
    end

  elseif dimfrom == 1
    if dimto == 0
      face_to_vertices_data = [1,2,4,5,1,4,2,5,2,3,3,5,3,6,6,5,7,8,4,7,5,8,8,9,6,9]
      face_to_vertices_ptrs = Int32[1,3,5,7,9,11,13,15,17,19,21,23,25,27]
      return Table(face_to_vertices_data, face_to_vertices_ptrs)
    elseif dimto == 1
      nedges = 13
      return identity_table(Int,Int32,nedges)
    elseif dimto == 2
      face_to_cells_data = [1,1,4,1,1,2,2,2,3,3,3,5,4,4,4,5,5,5]
      face_to_cells_ptrs = Int32[1,2,4,5,7,8,10,11,13,14,15,17,18,19]
      return Table(face_to_cells_data,face_to_cells_ptrs)
    else
      @unreachable
    end

  elseif dimfrom == 2
    if dimto == 0
      cell_to_vertices_data = [1,2,4,5,2,3,5,3,6,5,4,5,7,8,5,6,8,9]
      cell_to_vertices_ptrs = Int32[1,5,8,11,15,19]
      return Table(cell_to_vertices_data,cell_to_vertices_ptrs)
    elseif dimto == 1
      cell_to_faces_data = [1,2,3,4,5,4,6,7,6,8,2,9,10,11,8,12,11,13]
      cell_to_faces_ptrs = Int32[1,5,8,11,15,19]
      return Table(cell_to_faces_data, cell_to_faces_ptrs)
    elseif dimto == 2
      ncells = 5
      return identity_table(Int,Int32,ncells)
    else
      @unreachable
    end

  else
    @unreachable
  end
end

function get_vertex_node(g::DiscreteModelMock)
  collect(1:9)
end

function get_node_face_owner(g::DiscreteModelMock)
  collect(1:9)
end

function get_face_nodes(g::DiscreteModelMock,d::Integer)
  get_faces(g,d,0)
end

function get_isboundary_face(g::DiscreteModelMock,d::Integer)
  mask = fill(false, num_faces(g,d))
  if d == 0
    mask[[1,2,3,4,6,7,8,9]] .= true
  elseif d == 1
    mask[[1,5,3,7,9,10,12,13]] .= true
  elseif d == 2
    mask .= true
  else
    @unreachable
  end
  mask
end

function get_face_reffe_type(g::DiscreteModelMock,d::Integer)
  if d != 2
    return fill(Int8(1),num_faces(g,d))
  else
    quad = Int8(1)
    tri = Int8(2)
    [quad, tri, tri, quad, quad]
  end
end

function get_face_polytope_type(g::DiscreteModelMock,d::Integer)
  get_face_reffe_type(g,d)
end

function get_reffes(::Type{<:ReferenceFE{d}},g::DiscreteModelMock) where d
  [ LagrangianRefFE(Float64,p,1) for p in get_polytopes(Polytope{d},g) ]
end

function get_polytopes(::Type{<:Polytope{d}},g::DiscreteModelMock) where d
  if d == 0
    [VERTEX,]
  elseif d == 1
    [SEGMENT,]
  elseif d == 2
    [QUAD, TRI]
  else
    @unreachable
  end
end

function get_node_coordinates(g::DiscreteModelMock)
  Point{2,Float64}[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
end

function get_face_labeling(model::DiscreteModelMock)
  d_to_num_dfaces = [ num_faces(model,d) for d in 0:num_dims(model)]
  labels = FaceLabeling(d_to_num_dfaces)
  get_face_entity(labels,0) .= get_isboundary_face(model,0) .+ 1
  get_face_entity(labels,1) .= get_isboundary_face(model,1) .+ 1
  get_face_entity(labels,2) .= get_isboundary_face(model,2) .+ 1
  add_tag!(labels,"interior",[1,])
  add_tag!(labels,"boundary",[2,])
  add_tag_from_tags!(labels,"all",["interior","boundary"])
  labels
end

