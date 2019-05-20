function globaldofs(reffe::RefFE{D,T},
  gridgr::FullGridGraph,
  labels::FaceLabels,
  dirt::NTuple{N,Int}) where {D,T,N}
  dim_eqclass = Int[]
  c=1
  c_n = -1
  for vef_dim in 0:D
    vefcells= connections(gridgr,vef_dim,D)
    cellvefs= connections(gridgr,D,vef_dim)
    vef_labels = labels_on_dim(labels,vef_dim)
    num_vefs = length(vefcells)
    nfdofs=Array{Array{Int64},1}(undef,num_vefs)
    nfdofs_l = Int[]
    nfdofs_g = zeros(Int, num_vefs+1)
    nfdofs_g[1] = 1
    for (ignf,nf) in enumerate(vefcells)
      owner_cell = nf[1]
      lid_vef_dim = findfirst(i->i==ignf,cellvefs[owner_cell])
      lid_vef = reffe.polytope.nf_dim[end][vef_dim+1][1]+lid_vef_dim-1
      # @santiagobadia : Better a method for nfs of a particular type...
      num_nf_dofs = length(reffe.nfacedofs[lid_vef])
      if ( _is_fixed(vef_labels[ignf],dirt,labels) )
        nfdofs_l = Int[nfdofs_l..., c_n:-1:c_n-num_nf_dofs+1... ]
        c_n -= num_nf_dofs
      else
        nfdofs_l = Int[nfdofs_l..., c:c+num_nf_dofs-1... ]
        c += num_nf_dofs
      end
      nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
    end
    dim_eqclass = [ dim_eqclass..., CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g) ]
  end
  return [ dim_eqclass , c-1, -c_n-1 ]
end
