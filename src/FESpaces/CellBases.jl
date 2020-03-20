# CellBasis

"""
"""
abstract type CellBasis <: CellFieldLike end

"""
"""
FECellBasisStyle(::Type{<:CellBasis}) = Val{true}()

"""
"""
function TrialStyle(::Type{<:CellBasis})
  Val{false}()
end

TrialStyle(cb) = TrialStyle(typeof(cb))

"""
"""
is_trial(cb) = is_trial(typeof(cb))

"""
"""
is_test(cb) = ! is_trial(cb)

function is_trial(::Type{T}) where T
  get_val_parameter(TrialStyle(T))
end

function is_test(::Type{T}) where T
  ! is_trial(T)
end

"""
"""
function test_cell_basis(cf::CellBasis,args...;kwargs...)
  @test is_a_fe_cell_basis(cf)
  test_cell_field_like(cf,args...;kwargs...)
end

# Define how the metadata is preserved

function change_ref_style(cf::CellBasis)
  ref_sty = RefStyle(cf)
  bool = !get_val_parameter(ref_sty)
  new_sty = Val{bool}()
  trial_style = TrialStyle(cf)
  ar = get_array(cf)
  cm = get_cell_map(cf)
  GenericCellBasis(trial_style,ar,cm,new_sty)
end

function similar_object(cf::CellBasis,array::AbstractArray)
  cm = get_cell_map(cf)
  trial_style = TrialStyle(cf)
  GenericCellBasis(trial_style,array,cm,RefStyle(cf))
end

function similar_object(a::CellBasis,b::CellField,v::AbstractArray)
  similar_object(a,v)
end

function similar_object(a::CellField,b::CellBasis,v::AbstractArray)
  similar_object(b,v)
end

function similar_object(a::CellBasis,b::CellBasis,v::AbstractArray)
  _similar_cell_basis(v,a,b,TrialStyle(a),TrialStyle(b))
end

function _similar_cell_basis(v,a,b,a_trial::Val{T},b_trial::Val{T}) where T
  similar_object(a,v)
end

function _similar_cell_basis(v,a,b,a_trial::Val{false},b_trial::Val{true})
  _similar_cell_basis_test_trial(v,a,b)
end

function _similar_cell_basis(v,a,b,a_trial::Val{true},b_trial::Val{false})
  _similar_cell_basis_test_trial(v,b,a)
end

function _similar_cell_basis_test_trial(v,a,b)
  cm = get_cell_map(a)
  @assert is_in_ref_space(a) == is_in_ref_space(b)
  GenericCellMatrixField(v,cm,RefStyle(a))
end

# Operations

function gradient(cf::CellBasis)
  a = get_array(cf)
  g = field_array_gradient(a)
  similar_object(cf,g)
end

function grad2curl(cf::CellBasis)
  a = get_array(cf)
  g = grad2curl(UnimplementedField,a)
  similar_object(cf,g)
end

function operate(op,cf::CellBasis)
  a = get_array(cf)
  b = field_array_operation(UnimplementedField,op,a)
  similar_object(cf,b)
end

function operate(op,cf1::CellBasis,cf2::CellField)
  @assert length(cf1) == length(cf2)
  a1 = get_array(cf1)
  a2 = get_array(cf2)
  b = field_array_operation(UnimplementedField,op,a1,a2)
  similar_object(cf1,cf2,b)
end

function operate(op,cf1::CellField,cf2::CellBasis)
  @assert length(cf1) == length(cf2)
  a1 = get_array(cf1)
  a2 = get_array(cf2)
  b = field_array_operation(UnimplementedField,op,a1,a2)
  similar_object(cf1,cf2,b)
end

function operate(op,cf1::CellBasis,object)
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm,RefStyle(cf1))
  operate(op,cf1,cf2)
end

function operate(op,object,cf2::CellBasis)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm,RefStyle(cf2))
  operate(op,cf1,cf2)
end

function operate(op,a::CellBasis,b::CellBasis)
  _operate_cell_basis(op,a,b,TrialStyle(a),TrialStyle(b))
end

function _operate_cell_basis(op,a,b,atrial::Val{T},btrial::Val{T}) where T
  aa = get_array(a)
  ab = get_array(b)
  k = bcast(op)
  c = apply_to_field_array(UnimplementedField,k,aa,ab)
  similar_object(a,b,c)
end

function _operate_cell_basis(op,a,b,atrial::Val{false},btrial::Val{true})
  _operate_cell_basis_test_trial(op,a,b)
end

function _operate_cell_basis(op,a,b,atrial::Val{true},btrial::Val{false})
  _operate_cell_basis_test_trial(op,b,a)
end

function _operate_cell_basis_test_trial(op,cf1,cf2)
  @assert length(cf1) == length(cf2)
  a1 = get_array(cf1)
  a2 = get_array(cf2)
  b = field_array_operation(UnimplementedField,op,a1,a2)
  similar_object(cf1,cf2,b)
end

# Operations with extra arguments

function operate(op,a::CellField,b,objects...)
  arrays = map(get_array,(a,b,objects...))
  v = apply_to_field_array(UnimplementedField,bcast(op),arrays...)
  similar_object(a,v)
end

function operate(op,a::CellBasis,b::CellField,objects...)
  arrays = map(get_array,(a,b,objects...))
  v = apply_to_field_array(UnimplementedField,bcast(op),arrays...)
  similar_object(a,v)
end

function operate(op,a::CellField,b::CellBasis,objects...)
  arrays = map(get_array,(a,b,objects...))
  v = apply_to_field_array(UnimplementedField,bcast(op),arrays...)
  similar_object(b,v)
end

function operate(op,a::CellBasis,b::CellBasis,objects...)
  arrays = map(get_array,(a,b,objects...))
  v = apply_to_field_array(UnimplementedField,bcast(op),arrays...)
  r = similar_object(a,b,v)
  function fun(x)
    if isa(x,CellBasis) && is_trial(x)
      return true
    else
      return false
    end
  end
  @notimplementedif isa(r,CellMatrixField) && any( map(fun,objects) )
  r
end

# Concrete CellBases

"""
"""
struct GenericCellBasis{T,R} <: CellBasis
  trial_style::Val{T}
  array
  cell_map
  ref_trait::Val{R}
end

function GenericCellBasis(trial_style::Val{T},array::AbstractArray,cell_map::AbstractArray) where T
  ref_trait = Val{true}()
  GenericCellBasis(trial_style,array,cell_map,ref_trait)
end

get_array(a::GenericCellBasis) = a.array

get_cell_map(a::GenericCellBasis) = a.cell_map

function TrialStyle(::Type{<:GenericCellBasis{T}}) where T
  Val{T}()
end

function RefStyle(::Type{<:GenericCellBasis{T,R}}) where {T,R}
  Val{R}()
end

"""
"""
abstract type CellMatrixField <: CellFieldLike end

"""
"""
function test_cell_matrix_field(cb::CellMatrixField,args...; kwargs...)
  test_cell_field_like(cb,args...;kwargs...)
end

# Define how the metadata is preserved

function change_ref_style(cf::CellMatrixField)
  ref_sty = RefStyle(cf)
  bool = !get_val_parameter(ref_sty)
  new_sty = Val{bool}()
  ar = get_array(cf)
  cm = get_cell_map(cf)
  GenericCellMatrixField(ar,cm,new_sty)
end

function similar_object(cf::CellMatrixField,a::AbstractArray)
  cm = get_cell_map(cf)
  GenericCellMatrixField(a,cm,RefStyle(cf))
end

function similar_object(a::CellMatrixField,b::CellField,v::AbstractArray)
  similar_object(a,v)
end

function similar_object(a::CellField,b::CellMatrixField,v::AbstractArray)
  similar_object(b,v)
end

function similar_object(a::CellMatrixField,b::CellMatrixField,v::AbstractArray)
  similar_object(a,v)
end

# Define operations

function operate(op,a::CellMatrixField)
  aa = get_array(a)
  k = bcast(op)
  c = apply_to_field_array(UnimplementedField,k,aa)
  similar_object(a,c)
end

function operate(op,a::CellMatrixField,b::CellMatrixField)
  _operate_cell_matrix_field(op,a,b)
end

function operate(op,a::CellMatrixField,b::CellField)
  _operate_cell_matrix_field(op,a,b)
end

function operate(op,a::CellField,b::CellMatrixField)
  _operate_cell_matrix_field(op,a,b)
end

function _operate_cell_matrix_field(op,a,b)
  aa = get_array(a)
  ab = get_array(b)
  k = bcast(op)
  c = apply_to_field_array(UnimplementedField,k,aa,ab)
  similar_object(a,b,c)
end

function operate(op,a::CellMatrixField,b)
  cm = get_cell_map(a)
  _b = convert_to_cell_field(b,cm,RefStyle(a))
  operate(op,a,_b)
end

function operate(op,a,b::CellMatrixField)
  cm = get_cell_map(b)
  _a = convert_to_cell_field(a,cm,RefStyle(b))
  operate(op,_a,b)
end


# Concrete CellMatrixField

"""
"""
struct GenericCellMatrixField{A,B,R} <: CellMatrixField
  array::A
  cell_map::B
  ref_style::Val{R}
end

RefStyle(::Type{GenericCellMatrixField{A,B,R}}) where {A,B,R} = Val{R}()

get_array(a::GenericCellMatrixField) = a.array

get_cell_map(a::GenericCellMatrixField) = a.cell_map

# Restrictions

function restrict(cf::CellBasis,trian::Triangulation)
  _cf = to_ref_space(cf)
  a = get_array(_cf)
  r = restrict(a,trian)
  trial_style = TrialStyle(cf)
  _restrict_cell_basis(trial_style,r,trian)
end

function _restrict_cell_basis(trial_style,r::SkeletonPair,trian)
  cm = get_cell_map(trian)
  la = r.left
  ra = r.right
  l = GenericCellBasis(trial_style,la,cm)
  r = GenericCellBasis(trial_style,ra,cm)
  SkeletonCellBasis(trial_style,l,r)
end

function _restrict_cell_basis(trial_style,r::AbstractArray,trian)
  cm = get_cell_map(trian)
  GenericCellBasis(trial_style,r,cm)
end

# Integration of CellBasis and CellMatrixField

function integrate(cell_basis::CellBasis,trian::Triangulation,quad::CellQuadrature)
  cell_map = get_cell_map(trian)
  q = get_coordinates(quad)
  w = get_weights(quad)
  j = gradient(cell_map)
  @assert length(cell_basis) == length(cell_map) "Are you using the right triangulation to integrate?"
  @assert length(cell_basis) == length(w) "Are you using the right quadrature to integrate?"
  integrate(get_array(cell_basis),q,w,j)
end

function integrate(cell_basis::CellMatrixField,trian::Triangulation,quad::CellQuadrature)
  cell_map = get_cell_map(trian)
  q = get_coordinates(quad)
  w = get_weights(quad)
  j = gradient(cell_map)
  @assert length(cell_basis) == length(cell_map) "Are you using the right triangulation to integrate?"
  @assert length(cell_basis) == length(w) "Are you using the right quadrature to integrate?"
  integrate(get_array(cell_basis),q,w,j)
end

# Skeleton-related stuff

struct SkeletonCellBasis{T} <: GridapType
  trial_style::Val{T}
  left::CellBasis
  right::CellBasis
end

function Base.getproperty(x::SkeletonCellBasis, sym::Symbol)
  if sym == :inward
    x.left
  elseif sym == :outward
    x.right
  else
    getfield(x, sym)
  end
end

TrialStyle(::Type{<:SkeletonCellBasis{T}}) where T = Val{T}()

function get_cell_map(a::SkeletonCellBasis)
  get_cell_map(a.left)
end

function gradient(cf::SkeletonCellBasis)
  left = gradient(cf.left)
  right = gradient(cf.right)
  SkeletonCellBasis(cf.trial_style,left,right)
end

function grad2curl(cf::SkeletonCellBasis)
  left = grad2curl(cf.left)
  right = grad2curl(cf.right)
  SkeletonCellBasis(cf.trial_style,left,right)
end

function operate(op,cf::SkeletonCellBasis)
  left = operate(op,cf.left)
  right = operate(op,cf.right)
  SkeletonCellBasis(cf.trial_style,left,right)
end

function operate(op,cf1::SkeletonCellBasis,cf2::CellField)
  left = operate(op,cf1.left,cf2)
  right = operate(op,cf1.right,cf2)
  SkeletonCellBasis(cf1.trial_style,left,right)
end

function operate(op,cf1::CellField,cf2::SkeletonCellBasis)
  left = operate(op,cf1,cf2.left)
  right = operate(op,cf1,cf2.right)
  SkeletonCellBasis(cf2.trial_style,left,right)
end

function operate(op,cf1::SkeletonCellBasis,object)
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm,RefStyle(cf1.left))
  operate(op,cf1,cf2)
end

function operate(op,object,cf2::SkeletonCellBasis)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm,RefStyle(cf2.left))
  operate(op,cf1,cf2)
end

function jump(a::SkeletonCellBasis)
  ReducedSkeletonCellBasis(a.trial_style,a.left,-a.right)
end

function mean(a::SkeletonCellBasis)
  ReducedSkeletonCellBasis(a.trial_style,0.5*a.left,0.5*a.right)
end

# Result of reducing a SkeletonCellBasis

struct ReducedSkeletonCellBasis{T} <: GridapType
  trial_style::Val{T}
  left
  right
end

TrialStyle(::Type{<:ReducedSkeletonCellBasis{T}}) where T = Val{T}()

function get_cell_map(a::ReducedSkeletonCellBasis)
  get_cell_map(a.left)
end

function operate(op,cf::ReducedSkeletonCellBasis)
  left = operate(op,cf.left)
  right = operate(op,cf.right)
  ReducedSkeletonCellBasis(cf.trial_style,left,right)
end

function operate(op,cf1::ReducedSkeletonCellBasis,cf2::CellField)
  left = operate(op,cf1.left,cf2)
  right = operate(op,cf1.right,cf2)
  ReducedSkeletonCellBasis(cf1.trial_style,left,right)
end

function operate(op,cf1::CellField,cf2::ReducedSkeletonCellBasis)
  left = operate(op,cf1,cf2.left)
  right = operate(op,cf1,cf2.right)
  ReducedSkeletonCellBasis(cf2.trial_style,left,right)
end

function operate(op,cf1::ReducedSkeletonCellBasis,object)
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm,RefStyle(cf1.left))
  operate(op,cf1,cf2)
end

function operate(op,object,cf2::ReducedSkeletonCellBasis)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm,RefStyle(cf2.left))
  operate(op,cf1,cf2)
end

function operate(op,a::ReducedSkeletonCellBasis,b::ReducedSkeletonCellBasis)
  _operate_reduced_skeleton_cell_basis(op,a,b,a.trial_style,b.trial_style)
end

function  _operate_reduced_skeleton_cell_basis(
  op,a,b,a_trial::Val{T},b_trial::Val{T}) where T
  left = operate(op,a.left,b.left)
  right = operate(op,a.right,b.right)
  trial_style = Val{T}()
  ReducedSkeletonCellBasis(trial_style,left,right)
end

function  _operate_reduced_skeleton_cell_basis(
  op,a,b,a_trial::Val{false},b_trial::Val{true})
  _operate_skeleton_test_trial(op,a,b)
end

function  _operate_reduced_skeleton_cell_basis(
  op,a,b,a_trial::Val{true},b_trial::Val{false})
  _operate_skeleton_test_trial(op,b,a)
end

function _operate_skeleton_test_trial(op,a,b)
  ll = operate(op,a.left,b.left)
  lr = operate(op,a.left,b.right)
  rl = operate(op,a.right,b.left)
  rr = operate(op,a.right,b.right)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

struct SkeletonCellMatrixField <: GridapType
  ll
  lr
  rl
  rr
end

get_cell_map(a::SkeletonCellMatrixField) = get_cell_map(a.ll)

function operate(op,cf::SkeletonCellMatrixField)
  ll = operate(op,cf.ll)
  lr = operate(op,cf.lr)
  rl = operate(op,cf.rl)
  rr = operate(op,cf.rr)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

function operate(op,a::SkeletonCellMatrixField,b::SkeletonCellMatrixField)
  ll = operate(op,a.ll,b.ll)
  lr = operate(op,a.lr,b.lr)
  rl = operate(op,a.rl,b.rl)
  rr = operate(op,a.rr,b.rr)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

function operate(op,a::SkeletonCellMatrixField,b::CellField)
  ll = operate(op,a.ll,b)
  lr = operate(op,a.lr,b)
  rl = operate(op,a.rl,b)
  rr = operate(op,a.rr,b)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

function operate(op,a::CellField,b::SkeletonCellMatrixField)
  ll = operate(op,a,b.ll)
  lr = operate(op,a,b.lr)
  rl = operate(op,a,b.rl)
  rr = operate(op,a,b.rr)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

function operate(op,object,cf2::SkeletonCellMatrixField)
  cm = get_cell_map(cf2)
  cf1 = convert_to_cell_field(object,cm,RefStyle(cf2.ll))
  operate(op,cf1,cf2)
end

function operate(op,cf1::SkeletonCellMatrixField,object)
  cm = get_cell_map(cf1)
  cf2 = convert_to_cell_field(object,cm,RefStyle(cf1.ll))
  operate(op,cf1,cf2)
end

# Integration of Skeleton quantities

struct SkeletonMatrix{A} <: GridapType
  ll::A
  lr::A
  rl::A
  rr::A
end

struct SkeletonVector{A} <: GridapType
  left::A
  right::A
end

struct SkeletonCellMatrix <: GridapType
  ll
  lr
  rl
  rr
end

struct SkeletonCellVector <: GridapType
  left
  right
end

function integrate(
  a::SkeletonCellMatrixField,trian::Triangulation,quad::CellQuadrature)
  ll = integrate(a.ll,trian,quad)
  lr = integrate(a.lr,trian,quad)
  rl = integrate(a.rl,trian,quad)
  rr = integrate(a.rr,trian,quad)
  SkeletonCellMatrix(ll,lr,rl,rr)
end

function integrate(
  a::ReducedSkeletonCellBasis,trian::Triangulation,quad::CellQuadrature)
  left = integrate(a.left,trian,quad)
  right = integrate(a.right,trian,quad)
  SkeletonCellVector(left,right)
end

# Dirichlet related

function compute_dirichlet_cell_vector(cellmat,cellvals)
  k = DirichletVecKernel()
  apply(k,cellmat,cellvals)
end

"""
"""
function attach_dirichlet_bcs(cellmatvec,cellvals)
  k = DirichletMatVecKernel()
  apply(k,cellmatvec,cellvals)
end

struct DirichletVecKernel <: Kernel end

function kernel_cache(k::DirichletVecKernel,mat::AbstractMatrix,vals)
  vec = mat*vals
  CachedArray(vec)
end

@inline function apply_kernel!(cache,k::DirichletVecKernel,mat::AbstractMatrix,vals)
  n = size(mat,1)
  setsize!(cache,(n,))
  vec = cache.array
  mul!(vec,mat,vals)
  @inbounds for i in eachindex(vec)
    vec[i] = -vec[i]
  end
  vec
end

struct DirichletMatVecKernel <: Kernel
  k::DirichletVecKernel
  function DirichletMatVecKernel()
    k = DirichletVecKernel()
    new(k)
  end
end

function kernel_cache(k::DirichletMatVecKernel,mat,vals)
  kernel_cache(k.k,mat,vals)
end

function kernel_cache(k::DirichletMatVecKernel,matvec::Tuple,vals)
  mat, = matvec
  kernel_cache(k.k,mat,vals)
end

@inline function apply_kernel!(cache,k::DirichletMatVecKernel,mat,vals)
  vec = apply_kernel!(cache,k.k,mat,vals)
  (mat,vec)
end

@inline function apply_kernel!(cache,k::DirichletMatVecKernel,matvec::Tuple,vals)
  mat, vec = matvec
  vecd = apply_kernel!(cache,k.k,mat,vals)
  add_to_array!(vecd,vec)
  (mat, vecd)
end
