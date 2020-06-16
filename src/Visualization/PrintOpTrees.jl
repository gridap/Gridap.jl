
struct TreeNode{T}
  node::T
  showid::Bool
end

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

function print_op_tree(a,args...;kwargs...)
  print_op_tree(stdin,a,args...;kwargs...)
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

function get_children(n::TreeNode, a::Arrays.AppliedArray)
  (similar_tree_node(n,a.g),map(i->similar_tree_node(n,i),a.f)...)
end

function get_children(n::TreeNode, a::Arrays.SubVector)
  (similar_tree_node(n,a.vector),)
end

function get_children(n::TreeNode, a::Arrays.Reindexed)
  (similar_tree_node(n,a.i_to_v),similar_tree_node(n,a.j_to_i))
end

function get_children(n::TreeNode, a::Arrays.AppendedArray)
  (similar_tree_node(n,a.a),similar_tree_node(n,a.b))
end

function get_children(n::TreeNode, a::Arrays.ArrayPair)
  (similar_tree_node(n,a.a),similar_tree_node(n,a.b))
end

function get_children(n::TreeNode, a::Arrays.LocalToGlobalArray)
  (similar_tree_node(n,a.lid_to_gid),similar_tree_node(n,a.gid_to_val))
end

function get_children(n::TreeNode, a::Arrays.LocalToGlobalPosNegArray)
 (
  similar_tree_node(n,a.lid_to_gid),
  similar_tree_node(n,a.gid_to_val_pos),
  similar_tree_node(n,a.gid_to_val_neg)
 )
end

function get_children(n::TreeNode, a::Arrays.Table)
  (similar_tree_node(n,a.data),similar_tree_node(n,a.ptrs))
end

function get_children(n::TreeNode, a::FESpaces.ExtendedVector)
 (
  similar_tree_node(n,a.void_to_val),
  similar_tree_node(n,a.cell_to_val)
 )
end
