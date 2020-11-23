
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
function zero_dirichlet_values(f::SingleFieldFESpace)
  V = get_vector_type(f)
  allocate_vector(V,num_dirichlet_dofs(f))
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
  fe_basis = get_cell_shapefuns(f)
  @test isa(fe_basis,CellField)
  test_fe_space(f)
  cell_dofs = get_cell_dof_ids(f)
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
  @test isa(cell_dof_basis,CellDof)
end

function test_single_field_fe_space(f,matvecdata,matdata,vecdata,pred=(==))
  test_single_field_fe_space(f,pred)
  test_fe_space(f,matvecdata,matdata,vecdata)
end

# Some default API

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
function gather_dirichlet_values!(dirichlet_values,f::SingleFieldFESpace,cell_vals)
  free_values = zero_free_values(f)
  gather_free_and_dirichlet_values!(free_values,dirichlet_values,f,cell_vals)
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
function gather_free_values!(free_values,f::SingleFieldFESpace,cell_vals)
    dirichlet_values = zero_dirichlet_values(f)
    gather_free_and_dirichlet_values!(free_values,dirichlet_values,f,cell_vals)
    free_values
end

function CellField(fs::SingleFieldFESpace,cell_vals)
  v = get_cell_shapefuns(fs)
  cell_basis = get_cell_data(v)
  cell_field = lazy_map(linear_combination,cell_vals,cell_basis)
  GenericCellField(cell_field,get_triangulation(v),DomainStyle(v))
end

struct SingleFieldFEFunction{T<:CellField} <: FEFunction
  cell_field::T
  cell_dof_values::AbstractArray{<:AbstractVector{<:Number}}
  free_values::AbstractVector{<:Number}
  dirichlet_values::AbstractVector{<:Number}
  fe_space::SingleFieldFESpace
end

get_cell_data(f::SingleFieldFEFunction) = get_cell_data(f.cell_field)
get_triangulation(f::SingleFieldFEFunction) = get_triangulation(f.cell_field)
DomainStyle(::Type{SingleFieldFEFunction{T}}) where T = DomainStyle(T)

get_free_values(f::SingleFieldFEFunction) = f.free_values
get_cell_dof_values(f::SingleFieldFEFunction) = f.cell_dof_values
get_fe_space(f::SingleFieldFEFunction) = f.fe_space
Base.real(f::SingleFieldFEFunction) = FEFunction(f.fe_space,real(f.free_values),real(f.dirichlet_values))
Base.imag(f::SingleFieldFEFunction) = FEFunction(f.fe_space,imag(f.free_values),imag(f.dirichlet_values))

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

function FEFunction(fe::SingleFieldFESpace, free_values)
  diri_values = get_dirichlet_values(fe)
  FEFunction(fe,free_values,diri_values)
end

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
  s = get_cell_dof_basis(fs)
  trian = get_triangulation(s)
  f = CellField(object,trian,DomainStyle(s))
  cell_vals = s(f)
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
  dirichlet_values_scratch = zero_dirichlet_values(f)
  compute_dirichlet_values_for_tags!(dirichlet_values,dirichlet_values_scratch,f,tag_to_object)
end

function compute_dirichlet_values_for_tags!(
  dirichlet_values,
  dirichlet_values_scratch,
  f::SingleFieldFESpace,tag_to_object)

  dirichlet_dof_to_tag = get_dirichlet_dof_tag(f)
  _tag_to_object = _convert_to_collectable(tag_to_object,num_dirichlet_tags(f))
  for (tag, object) in enumerate(_tag_to_object)
    cell_vals = _cell_vals(f,object)
    fill!(dirichlet_values_scratch,zero(eltype(dirichlet_values_scratch)))
    gather_dirichlet_values!(dirichlet_values_scratch,f,cell_vals)
    _fill_dirichlet_values_for_tag!(dirichlet_values,dirichlet_values_scratch,tag,dirichlet_dof_to_tag)
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
  _convert_to_collectable(Fill(object,ntags),ntags)
end

function _convert_to_collectable(object::Number,ntags)
  _convert_to_collectable(Fill(object,ntags),ntags)
end
