
"""
"""
abstract type SingleFieldFESpace <: FESpace end

"""
"""
function num_dirichlet_dofs(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function get_cell_fe(fs::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function get_cell_dofs(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function zero_dirichlet_values(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function gather_free_and_dirichlet_values(f::SingleFieldFESpace,cell_vals)
  @abstractmethod
end

"""
"""
function scatter_free_and_dirichlet_values(f::SingleFieldFESpace,free_values,dirichlet_values)
  @abstractmethod
end

"""
"""
function test_single_field_fe_space(f::SingleFieldFESpace,cellmat,cellvec,pred=(==))
  test_fe_space(f,cellmat,cellvec)
  dirichlet_values = zero_dirichlet_values(f)
  @test length(dirichlet_values) == num_dirichlet_dofs(f)
  cell_fe = get_cell_fe(f)
  @test isa(cell_fe, SingleFieldCellFE)
  fe_basis = get_cell_fe_basis(f)
  @test isa(fe_basis,SingleFieldCellBasis)
  cell_dofs = get_cell_dofs(f)
  free_values = zero_free_values
  cell_vals = scatter_free_and_dirichlet_values(f,free_values,dirichlet_values)
  fv, dv = gather_free_and_dirichlet_values(f,cell_vals)
  @test pred(fv,free_values)
  @test pred(dv,dirichlet_values)
  fe_function = FEFunction(f,free_values,dirichlet_values)
  @test isa(fe_function, SingleFieldFEFunction)
end

"""
"""
function get_cell_map(fs::SingleFieldFESpace)
  cell_fe = get_cell_fe(fs)
  get_cell_map(cell_fe)
end

function get_cell_fe_basis(f::SingleFieldFESpace)
  cell_fe = get_cell_fe(f)
  cell_basis = get_cell_shapefuns(cell_fe)
  cell_map = get_cell_map(cell_fe)
  SingleFieldCellBasis(cell_basis,cell_map)
end

"""
    FEFunction(
      fs::SingleFieldFESpace, free_values::AbstractVector, dirichlet_values::AbstractVector)

The resulting FEFunction will be in the space if and only if `dirichlet_values`
are the ones provided by `get_dirichlet_values(fs)`
"""
function FEFunction(
  fs::SingleFieldFESpace, free_values::AbstractVector, dirichlet_values::AbstractVector)
  cell_vals = scatter_free_and_dirichlet_values(fs,free_values,diri_values)
  cell_fe = get_cell_fe(s)
  cell_basis = get_cell_shapefuns(cell_fe)
  cell_field = lincomb(cell_basis,cell_vals)
  SingleFieldFEFunction(cell_field,free_values,fs)
end

function FEFunction(fe::SingleFieldFESpace, free_values)
  diri_values = get_dirichlet_values(fe)
  FEFunction(fe,free_values,diri_values)
end

"""
"""
function get_dirichlet_values(f::SingleFieldFESpace)
  zero_dirichlet_values(f)
end

"""
"""
function gather_dirichlet_values(f::SingleFieldFESpace,cell_vals)
  _, dirichlet_values = gather_free_and_dirichlet_values(f,cell_vals)
  dirichlet_values
end

"""
"""
function gather_free_values(f::SingleFieldFESpace,cell_vals)
  free_values, _ = gather_free_and_dirichlet_values(f,cell_vals)
  free_values
end

"""
cell_field defined in the reference space with derivatives in the physical one
"""
function compute_free_and_dirichlet_values(f::SingleFieldFESpace, cell_field)
  cell_vals = _compute_cell_vals(f, cell_field)
  gather_free_and_dirichlet_values(f,cell_vals)
end

"""
"""
function compute_dirichlet_values(f::SingleFieldFESpace,cell_field)
  cell_vals = _compute_cell_vals(f, cell_field)
  gather_dirichlet_values(f,cell_vals)
end

"""
"""
function compute_free_values(f::SingleFieldFESpace,cell_field)
  cell_vals = _compute_cell_vals(f, cell_field)
  gather_free_values(f,cell_vals)
end

function _compute_cell_vals(f,cell_field)
  cell_fe = get_cell_fe(f)
  cell_dof_basis = get_cell_dof_basis(cell_fe)
  cell_dofs = get_cell_dofs(f)
  cell_vals = evaluate_dof_array(cell_dof_basis,cell_field)
  cell_vals
end

"""
The resulting FE function is in the space (in particular it fulfills Dirichlet BCs
even in the case that the given cell field does not fulfill them)
"""
function interpolate(fs::SingleFieldFESpace,object)
  cell_map = get_cell_map(fs)
  cell_field = _convert_to_interpolable(object,cell_map)
  free_values = compute_free_values(fs,cell_field)
  FEFunction(fs,free_values)
end

"""
"""
function compute_dirichlet_values_for_tags(f::SingleFieldFESpace,tag_to_object)
  dirichlet_dof_to_tag = get_dirichlet_dof_tag(f)
  cell_map = get_cell_map(f)
  dirichlet_values = zero_dirichlet_values(f)
  for (tag, object) in enumerate(_convert_to_collectable(tag_to_object))
    cell_field = _convert_to_interpolable(object,cell_map)
    dv = compute_dirichlet_values(f,cell_field)
    _fill_dirichlet_values_for_tag!(dirichlet_values,dv,tag,dirichlet_dof_to_tag)
  end
  dirichlet_values
end

function _fill_dirichlet_values_for_tag!(dirichlet_values,dv,tag,dirichlet_dof_to_tag)
  for (dof, v) in enumerate(dv)
    if dirichlet_dof_to_tag[dof] == tag
      dirichlet_values[dof] == v
    end
  end
end

function _convert_to_interpolable(object,cell_map)
  object
end

function _convert_to_interpolable(object::Function,cell_map)
  compose(object,cell_map)
end

function _convert_to_interpolable(object::Number,cell_map)
  Fill(object,length(cell_map))
end

function _convert_to_collectable(object)
  object
end

function _convert_to_collectable(object::Function)
  [object,]
end

# Ancillary types: SingleFieldFEFunction

"""
"""
struct SingleFieldFEFunction{A,B,C} <: FEFunction
  cell_field::A
  free_values::B
  fe_space::C
  @doc """
  """
  function SingleFieldFEFunction(
    cell_field::AbstractArray{<:Field},
    free_values::AbstractVector,
    fe_space::SingleFieldFESpace)

    A = typeof(cell_field)
    B = typeof(free_values)
    C = typeof(fe_space)
    new{A,B,C}(cell_field,free_values,fe_space)
  end
end

get_array(f::SingleFieldFEFunction) = f.cell_field

get_free_values(f::SingleFieldFEFunction) = f.free_values

get_fe_space(f::SingleFieldFEFunction) = f.fe_space

# Ancillary types: SingleFieldFEFunction

struct SingleFieldCellBasis{A,B} <: GridapType
  cell_basis::A
  cell_map::B
  @doc """
  """
  function SingleFieldCellBasis(
    cell_basis::AbstractArray{<:Field},
    cell_map::AbstractArray{<:Field})

    A = typeof(cell_basis)
    B = typeof(cell_map)
    new{A,B}(cell_basis,cell_map)
  end
end

function get_array(cb::SingleFieldCellBasis)
  cb.cell_basis
end

function get_cell_map(cb::SingleFieldCellBasis)
  cb.cell_map
end

# Ancillary types: SingleFieldCellFE

"""
"""
abstract type SingleFieldCellFE <: GridapType end

"""
"""
function get_cell_shapefuns(cf::SingleFieldCellFE)
  @abstractmethod
end

"""
"""
function get_cell_map(cf::SingleFieldCellFE)
  @abstractmethod
end

"""
"""
function get_reffes(cf::SingleFieldCellFE)
  @abstractmethod
end

"""
"""
function get_cell_type(cf::SingleFieldCellFE)
  @abstractmethod
end

"""
"""
function get_cell_dof_basis(cell_fe::SingleFieldCellFE)
  @abstractmethod
end

"""
"""
function test_single_field_cell_fe(f::SingleFieldCellFE)
  _ = get_cell_shapefuns(f)
  _ = get_reffes(f)
  _ = get_cell_map(f)
  _ = get_cell_type(f)
  _ = get_cell_dof_basis(f)
end

# Concrete implementations of SingleFieldCellFE

"""
"""
struct RefFEsWithMap{A,B,C,D} <: SingleFieldCellFE
  ctype_to_reffe::A
  cell_to_ctype::B
  cell_map::C
  cell_shapefuns::D
  @doc """
  """
  function RefFEsWithMap(
    ctype_to_reffe::Vector{<:ReferenceFE},
    cell_to_ctype::AbstractArray{<:Integer},
    cell_map::AbstractArray{<:Field})

    values =  map(get_shapefuns,ctype_to_reffe)
    ptrs = cell_to_ctype
    refshapefuns = CompressedArray(values,ptrs)
    cell_shapefuns = attachmap(refshapefuns,cell_map)

    A = typeof(ctype_to_reffe)
    B = typeof(cell_to_ctype)
    C = typeof(cell_map)
    D = typeof(cell_shapefuns)

    new{A,B,C,D}(ctype_to_reffe,cell_to_ctype,cell_map,cell_shapefuns)

  end
end

function get_cell_map(cf::RefFEsWithMap)
  cf.cell_map
end

function get_cell_shapefuns(cf::RefFEsWithMap)
  cf.cell_shapefuns
end

function get_reffes(cf::RefFEsWithMap)
  cf.ctype_to_reffe
end

function get_cell_type(cf::RefFEsWithMap)
  cf.cell_to_ctype
end

function get_cell_dof_basis(cf::RefFEsWithMap)
  reffes = cf.ctype_to_reffe
  values = map(get_dof_basis,reffes)
  ptrs = cf.cell_to_ctype
  #TODO in general we would need to attach a cell_map
  CompressedArray(values,ptrs)
end


