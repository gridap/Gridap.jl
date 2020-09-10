
"""
"""
abstract type SingleFieldFESpace <: FESpace end

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
function gather_free_and_dirichlet_values!(free_values, dirichlet_values,fs::SingleFieldFESpace,cell_vals)
  @abstractmethod
end

"""
"""
function test_single_field_fe_space(f::SingleFieldFESpace,pred=(==))
  fe_basis = get_cell_basis(f)
  @test isa(fe_basis,CellField)
  @test is_basis(fe_basis)
  test_fe_space(f)
  cell_dofs = get_cell_dofs(f)
  dirichlet_values = zero_dirichlet_values(f)
  @test length(dirichlet_values) == num_dirichlet_dofs(f)
  free_values = zero_free_values(f)
  cell_vals = scatter_free_and_dirichlet_values(f,free_values,dirichlet_values)
  fv, dv = gather_free_and_dirichlet_values(f,cell_vals)
  @test pred(fv,free_values)
  @test pred(dv,dirichlet_values)
  gather_free_and_dirichlet_values!(fv,dv,f,cell_vals)
  @test pred(fv,free_values)
  @test pred(dv,dirichlet_values)
  fv, dv = gather_free_and_dirichlet_values!(fv,dv,f,cell_vals)
  @test pred(fv,free_values)
  @test pred(dv,dirichlet_values)
  fe_function = FEFunction(f,free_values,dirichlet_values)
  @test isa(fe_function, SingleFieldFEFunction)
  test_fe_function(fe_function)
  ddof_to_tag = get_dirichlet_dof_tag(f)
  @test length(ddof_to_tag) == num_dirichlet_dofs(f)
  if length(get_dirichlet_dof_tag(f)) != 0
    @test maximum(get_dirichlet_dof_tag(f)) <= num_dirichlet_tags(f)
  end
  cell_dof_basis = get_cell_dof_basis(f)
  @test isa(CellField(f,get_cell_dofs(f)),CellField)
end

function test_single_field_fe_space(f,matvecdata,matdata,vecdata,pred=(==))
  test_single_field_fe_space(f,pred)
  test_fe_space(f,matvecdata,matdata,vecdata)
end

#"""
#"""
#function get_cell_map(fs::SingleFieldFESpace)
#  fe_basis = get_cell_basis(fs)
#  get_cell_map(fe_basis)
#end

function CellData.get_cell_axes(f::SingleFieldFESpace)
  get_cell_axes(get_cell_basis(f))
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
  cell_field = CellField(fs,cell_vals)
  SingleFieldFEFunction(cell_field,cell_vals,free_values,dirichlet_values,fs)
end

function CellData.CellField(fs::SingleFieldFESpace,cell_vals)
  _default_cell_field(fs,cell_vals)
end

function _default_cell_field(fs,cell_vals)
  cell_basis = get_cell_basis(fs)
  cell_field = lincomb(cell_basis,cell_vals)
  cell_field
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
function gather_free_and_dirichlet_values(f::SingleFieldFESpace,cell_vals)
  free_values = zero_free_values(f)
  dirichlet_values = zero_dirichlet_values(f)
  gather_free_and_dirichlet_values!(free_values,dirichlet_values,f,cell_vals)
end

"""
"""
function gather_dirichlet_values(f::SingleFieldFESpace,cell_vals)
  dirichlet_values = zero_dirichlet_values(f)
  gather_dirichlet_values!(dirichlet_values,f,cell_vals)
  dirichlet_values
end

"""
"""
function gather_free_values(f::SingleFieldFESpace,cell_vals)
  free_values = zero_free_values(f)
  gather_free_values!(free_values,f,cell_vals)
  free_values
end

"""
"""
function gather_dirichlet_values!(dirichlet_values,f::SingleFieldFESpace,cell_vals)
  free_values = zero_free_values(f)
  gather_free_and_dirichlet_values!(free_values,dirichlet_values,f,cell_vals)
  dirichlet_values
end

"""
"""
function gather_free_values!(free_values,f::SingleFieldFESpace,cell_vals)
    dirichlet_values = zero_dirichlet_values(f)
    gather_free_and_dirichlet_values!(free_values,dirichlet_values,f,cell_vals)
    free_values
end

@deprecate(
  interpolate(fs::SingleFieldFESpace, object),
  interpolate(object, fs::SingleFieldFESpace)
)

@deprecate(
  interpolate!(free_values,fs::SingleFieldFESpace, object),
  interpolate!(object, free_values,fs::SingleFieldFESpace)
)

@deprecate(
  interpolate_everywhere(fs::SingleFieldFESpace, object),
  interpolate_everywhere(object, fs::SingleFieldFESpace)
)

@deprecate(
  interpolate_everywhere!(free_values,dirichlet_values,fs::SingleFieldFESpace, object),
  interpolate_everywhere!(object, free_values,dirichlet_values,fs::SingleFieldFESpace)
)

@deprecate(
  interpolate_dirichlet(fs::SingleFieldFESpace, object),
  interpolate_dirichlet(object, fs::SingleFieldFESpace)
)

@deprecate(
  interpolate_dirichlet!(free_values,dirichlet_values,fs::SingleFieldFESpace, object),
  interpolate_dirichlet!(object, free_values,dirichlet_values,fs::SingleFieldFESpace)
)

"""
The resulting FE function is in the space (in particular it fulfills Dirichlet BCs
even in the case that the given cell field does not fulfill them)
"""
function interpolate(object, fs::SingleFieldFESpace)
  free_values = zero_free_values(fs)
  interpolate!(object, free_values,fs)
end

"""
"""
function interpolate!(object, free_values,fs::SingleFieldFESpace)
    cell_vals = _cell_vals(fs,object)
    gather_free_values!(free_values,fs,cell_vals)
    FEFunction(fs,free_values)
end

function _cell_vals(fs::SingleFieldFESpace,object)
  cdb = get_cell_dof_basis(fs)
  cf = convert_to_cell_field(object,length(cdb))
  cell_vals = evaluate(cdb,cf)
end

"""
like interpolate, but also compute new degrees of freedom for the dirichlet component.
The resulting FEFunction does not necessary belongs to the underlying space
"""
function interpolate_everywhere(object, fs::SingleFieldFESpace)
  free_values = zero_free_values(fs)
  dirichlet_values = zero_dirichlet_values(fs)
  interpolate_everywhere!(object, free_values,dirichlet_values,fs)
end

"""
"""
function interpolate_everywhere!(object, free_values,dirichlet_values,fs::SingleFieldFESpace)
  cell_vals = _cell_vals(fs,object)
  gather_free_and_dirichlet_values!(free_values,dirichlet_values,fs,cell_vals)
  FEFunction(fs,free_values,dirichlet_values)
end

"""
"""
function interpolate_dirichlet(object, fs::SingleFieldFESpace)
  free_values = zero_free_values(fs)
  dirichlet_values = zero_dirichlet_values(fs)
  interpolate_dirichlet!(object, free_values,dirichlet_values,fs)
end

"""
"""
function interpolate_dirichlet!(object, free_values,dirichlet_values,fs::SingleFieldFESpace)
  cell_vals = _cell_vals(fs,object)
  gather_dirichlet_values!(dirichlet_values,fs,cell_vals)
  fill!(free_values,zero(eltype(free_values)))
  FEFunction(fs,free_values, dirichlet_values)
end

"""
"""
function compute_dirichlet_values_for_tags(f::SingleFieldFESpace,tag_to_object)
  dirichlet_values = zero_dirichlet_values(f)
  compute_dirichlet_values_for_tags!(dirichlet_values,f,tag_to_object)
end

function compute_dirichlet_values_for_tags!(dirichlet_values,f::SingleFieldFESpace,tag_to_object)
  dirichlet_dof_to_tag = get_dirichlet_dof_tag(f)
  _tag_to_object = _convert_to_collectable(tag_to_object,num_dirichlet_tags(f))
  for (tag, object) in enumerate(_tag_to_object)
    cell_vals = _cell_vals(f,object)
    dv = gather_dirichlet_values(f,cell_vals)
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


function _convert_to_collectable(object::Number,ntags)
  _convert_to_collectable(fill(object,ntags),ntags)
end
