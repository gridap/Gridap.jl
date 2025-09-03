
"""
"""
struct TreeNode{T}
  node::T
  showid::Bool
end

"""
    similar_tree_node(n::TreeNode,node)
"""
similar_tree_node(n::TreeNode,node) = TreeNode(node,n.showid)

AbstractTrees.children(a::TreeNode) = get_children(a,a.node)

AbstractTrees.printnode(io::IO,a::TreeNode) = print_node(io,a,a.node)

get_children(n::TreeNode,a) = ()

function print_node(io,n::TreeNode,a)
  if length(get_children(n,a)) == 0
    print(io,typeof(a))
    if n.showid
      print(io," (")
      show(io,objectid(a))
      print(io,")")
    end
  else
    print(io,nameof(typeof(a)))
    if n.showid
      print(io," (")
      show(io,objectid(a))
      print(io,")")
    end
  end
end

"""
    print_op_tree(a,args...;kwargs...)
    print_op_tree(io::IO,a,args...;showid=false,kwargs...)

Print the operation tree of a lazy operation/array, in `stdout` if `io` is not given.
"""
function print_op_tree(a,args...;kwargs...)
  print_op_tree(stdout,a,args...;kwargs...)
end

function print_op_tree(io::IO,a,args...;showid=false,kwargs...)
  tree = TreeNode(a,showid)
  AbstractTrees.print_tree(io,tree,args...;kwargs...)
end

# Specialization

function get_children(n::TreeNode, a::Fill)
   (similar_tree_node(n,a.value),)
end

function get_children(n::TreeNode, a::CompressedArray)
  (similar_tree_node(n,a.values),similar_tree_node(n,a.ptrs))
end

function get_children(n::TreeNode, a::LazyArray)
  (similar_tree_node(n,a.maps),map(i->similar_tree_node(n,i),a.args)...)
end

function get_children(n::TreeNode, a::Reindex)
  (similar_tree_node(n,a.values),)
end

function get_children(n::TreeNode, a::PosNegReindex)
  (similar_tree_node(n,a.values_pos),similar_tree_node(n,a.values_neg))
end

function get_children(n::TreeNode, a::PosNegPartition)
  (similar_tree_node(n,a.i_to_iposneg),similar_tree_node(n,a.ipos_to_i),similar_tree_node(n,a.ineg_to_i))
end

function get_children(n::TreeNode, a::AppendedArray)
  (similar_tree_node(n,a.a),similar_tree_node(n,a.b))
end

function get_children(n::TreeNode, a::Table)
  (similar_tree_node(n,a.data),similar_tree_node(n,a.ptrs))
end

function get_children(n::TreeNode, a::VectorWithEntryRemoved)
  (similar_tree_node(n,a.a),)
end

function get_children(n::TreeNode, a::VectorWithEntryInserted)
  (similar_tree_node(n,a.a),)
end

