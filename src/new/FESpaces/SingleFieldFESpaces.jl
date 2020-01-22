
"""
"""
abstract type SingleFieldFESpace <: FESpace end

"""
"""
function get_cell_dofs(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function get_cell_dof_basis(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function num_dirichlet_dofs(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function zero_dirichlet_values(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function num_dirichlet_tags(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function get_dirichlet_dof_tag(f::SingleFieldFESpace)
  @abstractmethod
end

"""
"""
function scatter_free_and_dirichlet_values(f::SingleFieldFESpace,free_values,dirichlet_values)
  @abstractmethod
end

"""
"""
function gather_free_and_dirichlet_values(f::SingleFieldFESpace,cell_vals)
  @abstractmethod
end

"""
"""
function test_single_field_fe_space(f::SingleFieldFESpace,cellmat,cellvec,cellidsrows,cellidscols,pred=(==))
  fe_basis = get_cell_basis(f)
  @test isa(fe_basis,CellBasis)
  test_fe_space(f,cellmat,cellvec,cellidsrows,cellidscols)
  cell_dofs = get_cell_dofs(f)
  dirichlet_values = zero_dirichlet_values(f)
  @test length(dirichlet_values) == num_dirichlet_dofs(f)
  free_values = zero_free_values(f)
  cell_vals = scatter_free_and_dirichlet_values(f,free_values,dirichlet_values)
  fv, dv = gather_free_and_dirichlet_values(f,cell_vals)
  @test pred(fv,free_values)
  @test pred(dv,dirichlet_values)
  fe_function = FEFunction(f,free_values,dirichlet_values)
  @test isa(fe_function, SingleFieldFEFunction)
  test_fe_function(fe_function)
  if length(get_dirichlet_dof_tag(f)) == 0
    @test num_dirichlet_tags(f) == 0
  else
    @test maximum(get_dirichlet_dof_tag(f)) == num_dirichlet_tags(f)
  end
  cell_dof_basis = get_cell_dof_basis(f)
end

"""
"""
function get_cell_map(fs::SingleFieldFESpace)
  fe_basis = get_cell_basis(fs)
  get_cell_map(fe_basis)
end

"""
"""
function get_cell_shapefuns(fs::SingleFieldFESpace)
  fe_basis = get_cell_basis(fs)
  get_array(fe_basis)
end

"""
    FEFunction(
      fs::SingleFieldFESpace, free_values::AbstractVector, dirichlet_values::AbstractVector)

The resulting FEFunction will be in the space if and only if `dirichlet_values`
are the ones provided by `get_dirichlet_values(fs)`
"""
function FEFunction(
  fs::SingleFieldFESpace, free_values::AbstractVector, dirichlet_values::AbstractVector)
  cell_vals = scatter_free_and_dirichlet_values(fs,free_values,dirichlet_values)
  cell_shapefuns = get_cell_shapefuns(fs)
  cell_field = lincomb(cell_shapefuns,cell_vals)
  SingleFieldFEFunction(cell_field,free_values,dirichlet_values,fs)
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
function compute_free_and_dirichlet_values(f::SingleFieldFESpace, cell_field::CellField)
  cell_vals = _compute_cell_vals(f, cell_field)
  gather_free_and_dirichlet_values(f,cell_vals)
end

"""
"""
function compute_dirichlet_values(f::SingleFieldFESpace,cell_field::CellField)
  cell_vals = _compute_cell_vals(f, cell_field)
  gather_dirichlet_values(f,cell_vals)
end

"""
"""
function compute_free_values(f::SingleFieldFESpace,cell_field::CellField)
  cell_vals = _compute_cell_vals(f, cell_field)
  gather_free_values(f,cell_vals)
end

function _compute_cell_vals(f,cell_field)
  cell_dof_basis = get_cell_dof_basis(f)
  cell_dofs = get_cell_dofs(f)
  cell_vals = evaluate_dof_array(cell_dof_basis,get_array(cell_field))
  cell_vals
end

"""
The resulting FE function is in the space (in particular it fulfills Dirichlet BCs
even in the case that the given cell field does not fulfill them)
"""
function interpolate(fs::SingleFieldFESpace,object)
  cell_map = get_cell_map(fs)
  cell_field = convert_to_cell_field(object,cell_map)
  free_values = compute_free_values(fs,cell_field)
  FEFunction(fs,free_values)
end

"""
like interpolate, but also compute new degrees of freedom for the dirichlet component.
The resulting FEFunction does not necessary belongs to the underlying space
"""
function interpolate_everywhere(fs::SingleFieldFESpace,object)
  cell_map = get_cell_map(fs)
  cell_field = convert_to_cell_field(object,cell_map)
  free_values, dirichlet_values = compute_free_and_dirichlet_values(fs,cell_field)
  FEFunction(fs,free_values, dirichlet_values)
end

"""
"""
function interpolate_dirichlet(fs::SingleFieldFESpace,object)
  cell_map = get_cell_map(fs)
  cell_field = convert_to_cell_field(object,cell_map)
  dirichlet_values = compute_dirichlet_values(fs,cell_field)
  free_values = zero_free_values(fs)
  FEFunction(fs,free_values, dirichlet_values)
end

"""
"""
function compute_dirichlet_values_for_tags(f::SingleFieldFESpace,tag_to_object)
  dirichlet_dof_to_tag = get_dirichlet_dof_tag(f)
  cell_map = get_cell_map(f)
  dirichlet_values = zero_dirichlet_values(f)
  _tag_to_object = _convert_to_collectable(tag_to_object,num_dirichlet_tags(f))
  for (tag, object) in enumerate(_tag_to_object)
    cell_field = convert_to_cell_field(object,cell_map)
    dv = compute_dirichlet_values(f,cell_field)
    _fill_dirichlet_values_for_tag!(dirichlet_values,dv,tag,dirichlet_dof_to_tag)
  end
  dirichlet_values
end

function _fill_dirichlet_values_for_tag!(dirichlet_values,dv,tag,dirichlet_dof_to_tag)
  for (dof, v) in enumerate(dv)
    if dirichlet_dof_to_tag[dof] == tag
      dirichlet_values[dof] = v
    end
  end
end

function _convert_to_collectable(object,ntags)
  @assert ntags == length(object) "Incorrect number of dirichlet tags provided"
  object
end

function _convert_to_collectable(object::Function,ntags)
  _convert_to_collectable(fill(object,ntags),ntags)
end

function _convert_to_collectable(object::Integer,ntags)
  _convert_to_collectable(fill(object,ntags),ntags)
end

