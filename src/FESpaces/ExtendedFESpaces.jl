
"""
"""
struct ExtendedFESpace{S<:SingleFieldFESpace} <: SingleFieldFESpace
  space::S
  model::RestrictedDiscreteModel
  partition::PosNegPartition
  function ExtendedFESpace(space::SingleFieldFESpace,model::RestrictedDiscreteModel)
    model_portion = model.model
    @check get_triangulation(model_portion) === get_triangulation(space)
    @notimplementedif ConstraintStyle(space) == Constrained()
    full_cell_to_cell = Geometry.get_cell_to_parent_cell(model)
    partition = PosNegPartition(full_cell_to_cell,num_cells(get_parent_model(model)))
    new{typeof(space)}(space,model,partition)
  end
end

ConstraintStyle(::Type{ExtendedFESpace{S}}) where S = ConstraintStyle(S)

function get_cell_dof_ids(f::ExtendedFESpace)
  nfull, nvoid = Arrays.pos_and_neg_length(f.partition)
  fullcell_to_ids = get_cell_dof_ids(f.space)
  @check length(fullcell_to_ids) == nfull
  T = eltype(fullcell_to_ids)
  voidcell_to_ids = Fill(similar(T,0),nvoid)
  lazy_map( PosNegReindex(fullcell_to_ids,voidcell_to_ids), f.partition)
end

function get_cell_shapefuns(f::ExtendedFESpace)
  nfull, nvoid = Arrays.pos_and_neg_length(f.partition)
  dv = get_cell_shapefuns(f.space)
  data = get_cell_data(dv)
  fullcell_shapefuns = lazy_map(VoidBasisMap(false),data)
  @check length(fullcell_shapefuns) == nfull
  voidcell_shapefuns = Fill(VoidBasis(testitem(data),true),nvoid)
  a = lazy_map( PosNegReindex(fullcell_shapefuns,voidcell_shapefuns),f.partition)
  FEBasis(a,get_triangulation(f),TestBasis(),DomainStyle(dv))
end

function get_cell_dof_basis(f::ExtendedFESpace)
  nfull, nvoid = Arrays.pos_and_neg_length(f.partition)
  s = get_cell_dof_basis(f.space)
  data = get_cell_data(s)
  fullcell_dof_basis = lazy_map(VoidBasisMap(false),data)
  @check length(fullcell_dof_basis) == nfull
  voidcell_dof_basis = Fill(VoidBasis(testitem(data),true),nvoid)
  a = lazy_map( PosNegReindex(fullcell_dof_basis,voidcell_dof_basis),f.partition)
  CellDof(a,get_triangulation(f),DomainStyle(s))
end

function scatter_free_and_dirichlet_values(f::ExtendedFESpace,fv,dv)
  nfull, nvoid = Arrays.pos_and_neg_length(f.partition)
  fullcell_val = scatter_free_and_dirichlet_values(f.space,fv,dv)
  @check length(fullcell_val) == nfull
  T = eltype(fullcell_val)
  voidcell_val = Fill(similar(T,0),nvoid)
  lazy_map( PosNegReindex(fullcell_val,voidcell_val), f.partition)
end

function gather_free_and_dirichlet_values!(fv,dv,f::ExtendedFESpace,cv)
  _cv = lazy_map(Reindex(cv),f.partition.ipos_to_i)
  gather_free_and_dirichlet_values!(fv,dv,f.space,_cv)
end

# Delegated functions

function num_free_dofs(f::ExtendedFESpace)
  num_free_dofs(f.space)
end

function get_vector_type(f::ExtendedFESpace)
  get_vector_type(f.space)
end

function num_dirichlet_dofs(f::ExtendedFESpace)
  num_dirichlet_dofs(f.space)
end

function num_dirichlet_tags(f::ExtendedFESpace)
  num_dirichlet_tags(f.space)
end

function get_dirichlet_dof_tag(f::ExtendedFESpace)
  get_dirichlet_dof_tag(f.space)
end

function get_triangulation(f::ExtendedFESpace)
  get_triangulation(get_parent_model(f.model))
end

# Helpers

struct VoidBasisMap <: Map
  isvoid::Bool
end

@inline Arrays.evaluate!(cache,k::VoidBasisMap,b) = VoidBasis(b,k.isvoid)

struct VoidBasis{T,A} <: AbstractVector{T}
  basis::A
  isvoid::Bool
  function VoidBasis(basis::AbstractVector{T},isvoid::Bool) where T
    new{T,typeof(basis)}(basis,isvoid)
  end
end

function Base.size(a::VoidBasis)
  if a.isvoid
    (0,)
  else
    size(a.basis)
  end
end

Base.IndexStyle(::Type{<:VoidBasis}) = IndexLinear()

function Base.getindex(a::VoidBasis,i::Integer)
  if a.isvoid
    @unreachable "Unable to access 0-length array"
  else
    a.basis[i]
  end
end

Arrays.testitem(a::VoidBasis) = testitem(a.basis)

function Fields.return_cache(a::VoidBasis,x::Point) 
  cb = return_cache(a.basis,x)
  bx = return_value(a.basis,x)
  r = similar(bx,(0,))
  cb,r
end

function Fields.return_cache(a::VoidBasis,x::Field) 
  cb = return_cache(a.basis,x)
  bx = return_value(a.basis,x)
  r = similar(bx,(0,))
  cb,r
end

function Fields.return_cache(a::VoidBasis,x::AbstractVector{<:Point}) 
  cb = return_cache(a.basis,x)
  bx = return_value(a.basis,x)
  r = similar(bx,(length(x),0))
  cb,r
end

function Fields.return_cache(a::VoidBasis,v::AbstractVector{<:Field}) 
  cb = return_cache(a.basis,v)
  bx = return_value(a.basis,v)
  r = similar(bx,(0,length(v)))
  cb,r
end

for T in (:Point,:Field,:(AbstractVector{<:Point}),:(AbstractVector{<:Field}))
  @eval begin

    @inline function Fields.evaluate!(cache,a::VoidBasis,x::$T)
      cb, r = cache
      if a.isvoid
        r
      else
        evaluate!(cb,a.basis,x)
      end
    end

  end
end

function Fields.evaluate!(cache,k::Broadcasting{typeof(∇)},a::VoidBasis)
  VoidBasis(k(a.basis),a.isvoid)
end

function Fields.evaluate!(cache,k::Broadcasting{typeof(∇∇)},a::VoidBasis)
  VoidBasis(k(a.basis),a.isvoid)
end

