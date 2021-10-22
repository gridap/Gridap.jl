module TMP

using WriteVTK
using WriteVTK: VTKFile, VTKDataType, AbstractVTKDataset, XMLDocument
using WriteVTK: XMLElement, DatasetFile, AbstractFieldData, VTKUnstructuredGrid
using WriteVTK: create_root, xml_name, set_attribute, IS_LITTLE_ENDIAN
using WriteVTK: new_child, file_extension, add_extension, root
using WriteVTK: attribute, find_element, child_elements, datatype_str
using WriteVTK: save_file
import WriteVTK: vtk_save

struct ParallelDatasetFile <: VTKFile
    pvtkargs::Dict
    xdoc::XMLDocument
    dataset::DatasetFile
    path::AbstractString
    function ParallelDatasetFile(
        pvtkargs::Dict,
        xdoc::XMLDocument,
        dataset::DatasetFile,
        path::AbstractString)
        new(pvtkargs,xdoc,dataset,path)
    end
end

function _default_pvtkargs(pvtkargs)
    mandatory = [:part,:nparts]
    optional = [:ismain,:ghost_level]
    allkeys = vcat(mandatory,optional)

    d = Dict{Symbol,Any}(pvtkargs)
    for k in keys(d)
        if !(k in allkeys)
            throw(ArgumentError("$k is not a valid key in pvtkargs."))
        end
    end
    for k in mandatory
        if !(haskey(d,k))
            throw(ArgumentError("$k is a mandatory key in pvtkargs."))
        end
    end
    if ! haskey(d,:ismain)
      d[:ismain] = (d[:part] == 1)
    end
    if ! haskey(d,:ghost_level)
        d[:ghost_level] = 0
    end
    d
end

function vtk_save(vtk::ParallelDatasetFile)
  if isopen(vtk)
    _pvtk_grid_body(vtk)
    if vtk.pvtkargs[:ismain]
       save_file(vtk.xdoc, vtk.path)
    end
    vtk_save(vtk.dataset)
    close(vtk)
  end
  return [vtk.path] :: Vector{String}
end

function new_pdata_array(xParent,
                         type::Type{<:VTKDataType},
                         name::AbstractString,
                         Nc=nothing)
  xDA = new_child(xParent, "PDataArray")
  set_attribute(xDA, "type", datatype_str(type))
  set_attribute(xDA, "Name", name)
  if Nc != nothing
    set_attribute(xDA, "NumberOfComponents", Nc)
  end
end

function get_extension(filename::AbstractString)
  path, ext = splitext(filename)
  ext
end

function get_path(filename::AbstractString)
  path, ext = splitext(filename)
  path
end

function parallel_file_extension(g::AbstractVTKDataset)
  ext=file_extension(g)
  replace(ext,"."=>".p")
end

function parallel_xml_name(g::AbstractVTKDataset)
  "P"*xml_name(g)
end

function parallel_xml_write_header(pvtx::XMLDocument,dtype::AbstractVTKDataset)
  xroot = create_root(pvtx, "VTKFile")
  set_attribute(xroot, "type", parallel_xml_name(dtype))
  set_attribute(xroot, "version", "1.0")
  if IS_LITTLE_ENDIAN
    set_attribute(xroot, "byte_order", "LittleEndian")
  else
    set_attribute(xroot, "byte_order", "BigEndian")
  end
  xroot
end

function add_pieces!(grid_xml_node,
                     prefix::AbstractString,
                     extension::AbstractString,
                     num_pieces::Int)
  for i=1:num_pieces
    piece=new_child(grid_xml_node,"Piece")
    fn = _serial_filename(i,num_pieces,prefix,extension)
    set_attribute(piece,"Source",fn)
  end
  grid_xml_node
end

function add_ppoints!(grid_xml_node;
                      type::Type{<:VTKDataType}=Float64)
  ppoints=new_child(grid_xml_node,"PPoints")
  new_pdata_array(ppoints,type,"Points",3)
end

function add_ppoint_data!(pdataset::ParallelDatasetFile,
                          name::AbstractString;
                          type::Type{<:VTKDataType}=Float64,
                          Nc::Integer=1)
  xroot=root(pdataset.xdoc)
  dtype = xml_name_to_VTK_Dataset(pdataset.dataset.grid_type)
  grid_xml_node=find_element(xroot,parallel_xml_name(dtype))
  @assert grid_xml_node != nothing
  grid_xml_node_ppoint=find_element(grid_xml_node,"PPointData")
  if (grid_xml_node_ppoint==nothing)
    grid_xml_node_ppoint=new_child(grid_xml_node,"PPointData")
  end
  new_pdata_array(grid_xml_node_ppoint,type,name,Nc)
end

function add_pcell_data!(pdataset::ParallelDatasetFile,
                         name::AbstractString;
                         type::Type{<:VTKDataType}=Float64,
                         Nc=nothing)
  xroot=root(pdataset.xdoc)
  dtype = xml_name_to_VTK_Dataset(pdataset.dataset.grid_type)
  grid_xml_node=find_element(xroot,parallel_xml_name(dtype))
  @assert grid_xml_node != nothing
  grid_xml_node_pcell=find_element(grid_xml_node,"PCellData")
  if (grid_xml_node_pcell==nothing)
    grid_xml_node_pcell=new_child(grid_xml_node,"PCellData")
  end
  new_pdata_array(grid_xml_node_pcell,type,name,Nc)
end

function Base.setindex!(pvtk::ParallelDatasetFile,
                        data,
                        name::AbstractString,
                        loc::AbstractFieldData)
  pvtk.dataset[name,loc]=data
end

function Base.setindex!(vtk::ParallelDatasetFile, data, name::AbstractString)
  vtk.dataset[name]=data
end

function get_dataset_extension(dataset::DatasetFile)
  path, ext = splitext(dataset.path)
  @assert ext != ""
  ext
end

function get_dataset_path(dataset::DatasetFile)
  path, ext = splitext(dataset.path)
  @assert ext != ""
  path
end

function xml_name_to_VTK_Dataset(xml_name::AbstractString)
    if (xml_name=="ImageData")
      VTKImageData()
    elseif (xml_name=="RectilinearGrid")
      VTKRectilinearGrid()
    elseif (xml_name=="PolyData")
      VTKPolyData()
    elseif (xml_name=="StructuredGrid")
      VTKStructuredGrid()
    elseif (xml_name=="UnstructuredGrid")
      VTKUnstructuredGrid()
    else
      @assert false
    end
end

function num_children(element)
  length(collect(child_elements(element)))
end

function get_child_attribute(element,child,attr)
  children=collect(child_elements(element))
  attribute(children[child],attr)
end

function get_dataset_xml_type(dataset::DatasetFile)
  r=root(dataset.xdoc)
  attribute(r,"type")
end

function get_dataset_xml_grid_element(dataset::DatasetFile,
                                      element::AbstractString)
   r=root(dataset.xdoc)
   uns=find_element(r,get_dataset_xml_type(dataset))
   piece=find_element(uns,"Piece")
   find_element(piece,element)
end


const string_to_VTKDataType = Dict("Int8"=>Int8,
                                   "UInt8"=>UInt8,
                                   "Int16"=>Int16,
                                   "UInt16"=>UInt16,
                                   "Int32"=>Int32,
                                   "UInt32"=>UInt32,
                                   "Int64"=>Int64,
                                   "UInt64"=>UInt64,
                                   "Float32"=>Float32,
                                   "Float64"=>Float64,
                                   "String"=>String)

"""
    pvtk_grid(args...;pvtkargs,kwargs...)

Return a handler representing a parallel vtk file, which can be
eventually written to file with `vtk_save`.

Positional and keyword arguments in `args` and `kwargs`
are passed to `vtk_grid` verbatim (except file names that are augmented with the 
corresponding part id).

The extra keyword argument `pvtkargs` contains
extra options (as a `Dict{Symbol,Any}` or a `Vector{Pair{Symbol,Any}}`)
that only apply for parallel vtk file formats.

Mandatory keys in `pvtkargs` are:

- `:part` current (1-based) part id
- `:nparts` total number of parts

Default key/value pairs in `pvtkargs` are
- `:ismain=>part==1` True if the current part id `part` is the main (the only one that will write the .pvtk file).
- `:ghost_level=>0` Ghost level.
"""
function pvtk_grid(filename::AbstractString,args...;pvtkargs,kwargs...)
    _pvtkargs = _default_pvtkargs(pvtkargs)
    part = _pvtkargs[:part]
    nparts = _pvtkargs[:nparts]
    bname = basename(filename)
    path = mkpath(filename)
    prefix = joinpath(path,bname)
    extension = ""
    fn = _serial_filename(part,nparts,prefix,extension)
    vtk = vtk_grid(fn,args...;kwargs...)
    _pvtk_grid_header(_pvtkargs,vtk,path)
end

function _serial_filename(part,nparts,prefix,extension)
  p = lpad(part,ceil(Int,log10(nparts)),'0')
  fn = prefix*"_$p"*extension
end

function _pvtk_grid_header(
    pvtkargs::Dict,
    dataset::DatasetFile,
    filename::AbstractString)

    pvtx  = XMLDocument()
    dtype = xml_name_to_VTK_Dataset(dataset.grid_type)
    xroot = parallel_xml_write_header(pvtx,dtype)

    # Generate parallel grid node
    grid_xml_node=new_child(xroot,"P"*dataset.grid_type)
    set_attribute(grid_xml_node, "GhostLevel", string(pvtkargs[:ghost_level]))

    prefix=joinpath(filename,basename(filename))
    extension=file_extension(dtype)
    add_pieces!(grid_xml_node,prefix,extension,pvtkargs[:nparts])

    pextension = parallel_file_extension(dtype)
    pfilename = add_extension(filename,pextension)
    pdataset=ParallelDatasetFile(pvtkargs,pvtx,dataset,pfilename)

    # Generate PPoints
    points=get_dataset_xml_grid_element(dataset,"Points")
    type=get_child_attribute(points,1,"type")
    add_ppoints!(grid_xml_node,type=string_to_VTKDataType[type])

    pdataset
  end

function _pvtk_grid_body(pdataset::ParallelDatasetFile)
    dataset=pdataset.dataset
    # Generate PPointData
    pointdata=get_dataset_xml_grid_element(dataset,"PointData")
    if pointdata != nothing
      if (num_children(pointdata)>0)
        for child=1:num_children(pointdata)
          name=get_child_attribute(pointdata,child,"Name")
          type=get_child_attribute(pointdata,child,"type")
          Nc=get_child_attribute(pointdata,child,"NumberOfComponents")
          add_ppoint_data!(pdataset,
                           name;
                           type=string_to_VTKDataType[type],
                           Nc=parse(Int64,Nc))
        end
      end
    end

    # Generate PCellData
    celldata=get_dataset_xml_grid_element(dataset,"CellData")
    if celldata!=nothing
      if (num_children(celldata)>0)
        for child=1:num_children(celldata)
          name=get_child_attribute(celldata,child,"Name")
          type=get_child_attribute(celldata,child,"type")
          Nc=get_child_attribute(celldata,child,"NumberOfComponents")
          add_pcell_data!(pdataset,
                          name;
                          type=string_to_VTKDataType[type],
                          Nc=(Nc==nothing ? nothing : parse(Int64,Nc)))
        end
      end
    end
    pdataset
  end

end #module
