function StructHexMesh(nparts,nprocs=-1,iproc=0)
  D = length(nparts)
  if (iproc == 0)
    iproc = tuple(zeros(Int64,D)...)
  end
  extrusion = Polytopes.PointInt{D}(ones(Int64,D))
  polytope = Polytope(extrusion)
  numnfs=length(polytope.nfaces)
  nvefs=length(polytope.nfaces)-1
  aux = 2*ones(Int64,D)
  offst = [ prod(aux[1:i-1]) for i=1:D]
  nfspol = hcat([polytope.nfaces[i].extrusion for i=1:nvefs+1]...)
  nftype = [offst'*polytope.nfaces[i].extrusion for i=1:nvefs+1].+1
  spdimv=2*ones(Int64,D)
  exttypes = Array{Array{Int64,1}}(undef,spdimv...)
  cartesianindexmatrix!(exttypes)
  exttypes = hcat(reshape(exttypes, length(exttypes))...)'
  dimtypes = [sum(exttypes[i,:]) for i=1:size(exttypes,1)]
  extli = hcat([(nparts+flip.(exttypes[i,:])) for i=1:size(exttypes,1)]...)
  offnf= hcat([[prod(extli[1:i-1,j]) for i=1:D] for j=1:size(exttypes,1)]...)
  extli = [ prod(extli[1:D,j]) for j=1:size(exttypes,1)]
  perm = sortperm(dimtypes)
  invperm = sortperm(perm)
  offst = accumulate(+, extli[perm])
  offst[2:end]=offst[1:end-1]; offst[1]=0
  offst = offst[invperm]
  offstnf=offst[nftype]
  cells = Array{Array{Int64,1}}(undef,nparts...)
  cartesianindexmatrixoffset!(cells,iproc)
  cells=reshape(cells,length(cells))
  numcells=length(cells)
  cellvefs=gidscellxtype(cells,polytope,nparts)
  for ivef=1:nvefs+1
    cellvefs[ivef,:].+=offstnf[ivef]
  end
  cellvefs = vec(cellvefs)
  cellvefs_arr = CellVectorFromDataAndStride(cellvefs, numnfs)
  vefcells = [ Int64[] for i=1:offstnf[nvefs+1]+numcells]
  for (cell,cvefs) in enumerate(cellvefs_arr)
    for iface=1:nvefs+1
      push!(vefcells[cvefs[iface]],cell)
    end
  end
  # @santiagobadia : Why don't we create a List such that we add entries to the
  # list and it finally creates the data an ptrs arrays? I think it is the same
  # as the coo to csr but without matrix entries... Is there a "symbolic" version
  # in Julia?
  nvefsg = offstnf[nvefs+1]+numcells
  vefcells_l = []
  vefcells_p = zeros(Int, nvefsg+1)
  for (cell,cvefs) in enumerate(cellvefs_arr)
    for (lnf,gnf) in enumerate(cvefs)
      vefcells_p[gnf+1] += 1
    end
  end
  vefcells_l = zeros(Int, sum(vefcells_p))
  vefcells_p[1] = 1
  vefcells_p
  for p in 2:length(vefcells_p)
    vefcells_p[p] += vefcells_p[p-1]
  end
  aux = copy(vefcells_p)
  for (cell,cvefs) in enumerate(cellvefs_arr)
    for (lnf,gnf) in enumerate(cvefs)
      vefcells_l[aux[gnf]] = cell
      aux[gnf] +=1
    end
  end
  vefcells_l = vec(vefcells_l)
  vefcells_p = vec(vefcells_p)
  vefcells_arr = CellVectorFromDataAndPtrs(vefcells_l, vefcells_p)
  # coordinates = zeros(Float64,2,2)
  cellvefs_arr
  vefcells_arr
  # StructHexMesh{D}(cellvefs_arr,vefcells_arr,coordinates,polytope)
  StructHexMesh{D}(polytope,cellvefs_arr,vefcells_arr)
end

# function StructHexMesh(nparts,nprocs=-1,iproc=0)
#   D = length(nparts)
#   if (iproc == 0)
#     iproc = tuple(zeros(Int64,D)...)
#   end
#   extrusion = Polytopes.PointInt{D}(ones(Int64,D))
#   polytope = Polytope(extrusion)
#   numnfs=length(polytope.nfaces)
#   nvefs=length(polytope.nfaces)-1
#   aux = 2*ones(Int64,D)
#   offst = [ prod(aux[1:i-1]) for i=1:D]
#   # offst = linfs.offset
#   nfspol = hcat([polytope.nfaces[i].extrusion for i=1:nvefs+1]...)
#   # nftypep1 = lex2int(linfs,nfspol)
#   nftype = [offst'*polytope.nfaces[i].extrusion for i=1:nvefs+1].+1
#   spdimv=2*ones(Int64,D)
#   exttypes = Array{Array{Int64,1}}(undef,spdimv...)
#   cartesianindexmatrix!(exttypes)
#   exttypes = hcat(reshape(exttypes, length(exttypes))...)'
#   dimtypes = [sum(exttypes[i,:]) for i=1:size(exttypes,1)]
#   extli = hcat([(nparts+flip.(exttypes[i,:])) for i=1:size(exttypes,1)]...)
#   offnf= hcat([[prod(extli[1:i-1,j]) for i=1:D] for j=1:size(exttypes,1)]...)
#   extli = [ prod(extli[1:D,j]) for j=1:size(exttypes,1)]
#   perm = sortperm(dimtypes)
#   invperm = sortperm(perm)
#   offst = accumulate(+, extli[perm])
#   offst[2:end]=offst[1:end-1]; offst[1]=0
#   offst = offst[invperm]
#   offstnf=offst[nftype]
#   cells = Array{Array{Int64,1}}(undef,nparts...)
#   cartesianindexmatrixoffset!(cells,iproc)
#   cells=reshape(cells,length(cells))
#   numcells=length(cells)
#   cellvefs=gidscellxtype(cells,polytope,nparts)
#   for ivef=1:nvefs+1
#     cellvefs[ivef,:].+=offstnf[ivef]
#   end
#   cellvefs = [cellvefs[:,icell] for icell=1:size(cellvefs,2)]
#   vefcells = [ Int64[] for i=1:offstnf[nvefs+1]+numcells]
#   for cell=1:numcells
#     for iface=1:nvefs+1
#       push!(vefcells[cellvefs[cell][iface]],cell)
#     end
#   end
#   # for cell=1:length(cells)
#   #   cellvefs[:,cell]+=offstnf
#   # end
#   coordinates = zeros(Float64,2,2)
#   StructHexMesh{D}(cellvefs,vefcells,coordinates,polytope)
# end

function lex2int(liset::LexIndexSet,points::Array{Int64,1})
  return liset.offset'*points.+1
end

function gidscellxtype(cells,polytope,nparts)
  licells = LexIndexSet(nparts)
  numnfs=length(polytope.nfaces)
  gid1=zeros(Int64,numnfs)
  gid=zeros(Int64,numnfs,length(cells))
  offst = [nparts+flip.(polytope.nfaces[i].extrusion) for i=1:numnfs]
  for inface=1:numnfs
    offst[inface] = [ prod(offst[inface][1:i-1]) for i=1:length(offst[inface])]
    nface=polytope.nfaces[inface]
    anc=nface.anchor
    gid1[inface] = offst[inface]'*anc.+1
  end
  cellvefs = Array{Int64,2}(undef,numnfs,length(cells))
  for i=1:numnfs; cellvefs[i,:].=gid1[i]; end
  cellvefs+=hcat(offst...)'*hcat(cells...)
  return cellvefs
end

function flip(a)
  if (a==0)
    a=1
  elseif a==1
    a=0
  end
end

@generated function cartesianindexmatrixoffset!(A::Array{Array{Int64,M},N},offst) where {M,N}
  quote
    offst=offst.-1
    @nloops $N i A begin (@nref $N A i) = [map(+,(@ntuple $N i),offst)...] end
  end
end

@generated function cartesianindexmatrix!(A::Array{Array{T,M},N}) where {T,M,N}
  quote
    @nloops $N i A begin (@nref $N A i) = [((@ntuple $N i).-1)...] end
  end
end
