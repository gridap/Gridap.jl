"""
"""
abstract type DiscreteModel{Dc,Dp} end

"""
"""
function num_faces(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
"""
function num_nodes(g::DiscreteModel)
  @abstractmethod
end

"""
"""
function get_faces(g::DiscreteModel,dimfrom::Integer,dimto::Integer)
  @abstractmethod
end

"""
"""
function get_vertex_node(g::DiscreteModel)
  @abstractmethod
end

"""
"""
function get_node_vertex(g::DiscreteModel)
  @abstractmethod
end

"""
"""
function get_face_nodes(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
"""
function get_faces_around_node(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
"""
function get_isboundary_face(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
"""
function get_isboundary_node(g::DiscreteModel)
  @abstractmethod
end

"""

Index to the vector get_reffes(g,d)
"""
function get_face_reffe_type(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
"""
function get_face_polytope_type(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
"""
function get_reffes(::Type{<:ReferenceFE{d}},g::DiscreteModel) where d
  @abstractmethod
end

"""
Index to the vector get_polytopes(g,d)
"""
function get_polytopes(::Type{<:Polytope{d}},g::DiscreteModel) where d
  @abstractmethod
end

"""
"""
function get_vertex_coordinates(g::DiscreteModel)
  @abstractmethod
end

"""
"""
function get_node_coordinates(g::DiscreteModel)
  @abstractmethod
end

# Default API

"""
"""
num_dims(g::DiscreteModel{Dc}) where Dc = Dc

num_dims(::Type{<:DiscreteModel{Dc}}) where Dc = Dc

"""
"""
num_cell_dims(g::DiscreteModel{Dc}) where Dc = Dc

num_cell_dims(::Type{<:DiscreteModel{Dc}}) where Dc = Dc

"""
"""
num_point_dims(g::DiscreteModel{Dc,Dp}) where {Dc,Dp} = Dp

num_point_dims(::Type{<:DiscreteModel{Dc,Dp}}) where {Dc,Dp} = Dp

#get_dimranges(g)
#get_offsets(g)
#get_offset(g,d)

#get_dface_to_faces(g,d)
#get_faces_around_dface(g,d)
#
#
#num_faces(g)
#num_cells(g)
#num_facets(g)
#num_edges(g)
#num_vertices(g)

"""
"""
function test_discrete_model(model::DiscreteModel{Dc,Dp}) where {Dc,Dp}

end

