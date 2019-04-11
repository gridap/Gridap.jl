
struct VisualizationGrid{D,Z,G<:Grid{D,Z},C<:IndexCellValue{Int},P<:CellPoints{Z}} <: Grid{D,Z}
  grid::G
  coarsecells::C
  samplingpoints::P
end

points(vg::VisualizationGrid) = points(vg.grid)

cells(vg::VisualizationGrid) = cells(vg.grid)

celltypes(vg::VisualizationGrid) = celltypes(vg.grid)

function writevtk(vg::VisualizationGrid,filebase;celldata=Dict(),cellfields=Dict())
  cdata = prepare_cdata(celldata,vg.coarsecells)
  pdata = prepare_pdata(cellfields,vg.samplingpoints)
  writevtk(vg.grid,filebase,celldata=cdata,pointdata=pdata)
end

function prepare_cdata(celldata,fine_to_coarse)
  cdata = Dict()
  for (k,v) in celldata
    acoarse = collect(v)
    afine = allocate_afine(acoarse,length(fine_to_coarse))
    fill_afine!(afine,acoarse,fine_to_coarse)
    cdata[k] = afine
  end
  k2 = "cellid"
  @assert ! haskey(cdata,k2)
  cdata[k2] = fine_to_coarse
  cdata
end

allocate_afine(acoarse::Array{T},l) where T = Array{T,1}(undef,(l,))

function fill_afine!(afine,acoarse,fine_to_coarse)
  for (i,coarse) in enumerate(fine_to_coarse)
    afine[i] = acoarse[coarse]
  end
end

function prepare_pdata(cellfields,samplingpoints)
  pdata = Dict()
  for (k,v) in cellfields
    pdata[k] = collect(flatten(evaluate(v,samplingpoints)))
  end
  pdata
end

