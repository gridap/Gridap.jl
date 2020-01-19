# Abstract types

"""
"""
abstract type CellBasis <: CellField end

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

is_test(cb) = ! is_trial(cb)

function is_trial(::Type{T}) where T
  get_val_parameter(TrialStyle(T))
end

function is_test(::Type{T}) where T
  ! is_trial(T)
end

"""
"""
function test_cell_basis(cb::CellBasis,args...; kwargs...)
  @test is_a_fe_cell_basis(cb)
  test_cell_field(cb,args...;kwargs...)
end

"""
"""
abstract type CellMatrixField <: CellField end

"""
"""
function test_cell_matrix_field(cb::CellMatrixField,args...; kwargs...)
  test_cell_field(cb,args...;kwargs...)
end

# Define how the metadata is preserved

function similar_cell_field(cf::CellBasis,a::AbstractArray)
  cm = get_cell_map(cf)
  trial_style = TrialStyle(cf)
  GenericCellBasis(trial_style,a,cm)
end

function similar_cell_field(cf::CellMatrixField,a::AbstractArray)
  cm = get_cell_map(cf)
  GenericCellMatrixField(a,cm)
end

function similar_cell_field(a::CellBasis,b::CellField,v::AbstractArray)
  similar_cell_field(a,v)
end

function similar_cell_field(a::CellField,b::CellBasis,v::AbstractArray)
  similar_cell_field(b,v)
end

function similar_cell_field(::CellBasis,::CellMatrixField,::AbstractArray)
  @notimplemented
end

function similar_cell_field(::CellMatrixField,::CellBasis,a::AbstractArray)
  @notimplemented
end

function similar_cell_field(a::CellMatrixField,b::CellField,v::AbstractArray)
  similar_cell_field(a,v)
end

function similar_cell_field(a::CellField,b::CellMatrixField,v::AbstractArray)
  similar_cell_field(b,v)
end

function similar_cell_field(a::CellMatrixField,b::CellMatrixField,v::AbstractArray)
  similar_cell_field(a,v)
end

function similar_cell_field(a::CellBasis,b::CellBasis,v::AbstractArray)
  _similar_cell_basis(v,a,b,TrialStyle(a),TrialStyle(b))
end

function _similar_cell_basis(v,a,b,a_trial::Val{T},b_trial::Val{T}) where T
  similar_cell_field(a,v)
end

function _similar_cell_basis(v,a,b,a_trial::Val{false},b_trial::Val{true})
  _similar_cell_basis_test_trial(v,a,b)
end

function _similar_cell_basis(v,a,b,a_trial::Val{true},b_trial::Val{false})
  _similar_cell_basis_test_trial(v,b,a)
end

function _similar_cell_basis_test_trial(v,a,b)
  cm = get_cell_map(a)
  GenericCellMatrixField(v,cm)
end

# Define operations

function operate_cell_field(op,a::CellBasis,b::CellMatrixField)
  @notimplemented
end

function operate_cell_field(op,a::CellMatrixField,b::CellBasis)
  @notimplemented
end

function operate_cell_field(op,a::CellBasis,b::CellBasis)
  _operate_cell_basis(op,a,b,TrialStyle(a),TrialStyle(b))
end

function _operate_cell_basis(op,a,b,atrial::Val{T},btrial::Val{T}) where T
  aa = get_array(a)
  ab = get_array(b)
  k = bcast(op)
  c = apply_to_field_array(UnimplementedField,k,aa,ab)
  similar_cell_field(a,b,c)
end

function _operate_cell_basis(op,a,b,atrial::Val{false},btrial::Val{true})
  _operate_cell_basis_test_trial(op,a,b)
end

function _operate_cell_basis(op,a,b,atrial::Val{true},btrial::Val{false})
  _operate_cell_basis_test_trial(op,b,a)
end

function _operate_cell_basis_test_trial(op,a,b)
  operate_cell_field_default(op,a,b)
end

function operate_cell_field(op,a::CellMatrixField,b::CellMatrixField)
  _operate_cell_matrix_field(op,a,b)
end

function operate_cell_field(op,a::CellMatrixField,b::CellField)
  _operate_cell_matrix_field(op,a,b)
end

function operate_cell_field(op,a::CellField,b::CellMatrixField)
  _operate_cell_matrix_field(op,a,b)
end

function _operate_cell_matrix_field(op,a,b)
  aa = get_array(a)
  ab = get_array(b)
  k = bcast(op)
  c = apply_to_field_array(UnimplementedField,k,aa,ab)
  similar_cell_field(a,b,c)
end

# Concrete types

"""
"""
struct GenericCellBasis{T} <: CellBasis
  trial_style::Val{T}
  array
  cell_map
end

"""
"""
function GenericCellBasis(array::AbstractArray,cell_map::AbstractArray)
  trial_style = Val{false}()
  GenericCellBasis(trial_style,array,cell_map)
end

get_array(a::GenericCellBasis) = a.array

get_cell_map(a::GenericCellBasis) = a.cell_map

function TrialStyle(::Type{<:GenericCellBasis{T}}) where T
  Val{T}()
end

"""
"""
struct GenericCellMatrixField{A,B} <: CellMatrixField
  array::A
  cell_map::B
end

get_array(a::GenericCellMatrixField) = a.array

get_cell_map(a::GenericCellMatrixField) = a.cell_map

# Skeleton-related stuff

function restrict(cf::CellBasis,trian::Triangulation)
  a = get_array(cf)
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

struct SkeletonCellBasis{T,A,B} <: GridapType
  trial_style::Val{T}
  left::A
  right::B
end

TrialStyle(::Type{<:SkeletonCellBasis{T}}) where T = Val{T}()

function jump(a::SkeletonCellBasis)
  ReducedSkeletonCellBasis(a.trial_style,a.left,-a.right)
end

function mean(a::SkeletonCellBasis)
  ReducedSkeletonCellBasis(a.trial_style,0.5*a.left,0.5*a.right)
end

struct ReducedSkeletonCellBasis{T,A,B} <: CellBasis
  trial_style::Val{T}
  left::A
  right::B
end

TrialStyle(::Type{<:ReducedSkeletonCellBasis{T}}) where T = Val{T}()

function get_cell_map(a::ReducedSkeletonCellBasis)
  get_cell_map(a.left)
end

function get_array(a::ReducedSkeletonCellBasis)
  @unreachable
end

function similar_cell_field(cf::ReducedSkeletonCellBasis,a::AbstractArray)
  @unreachable
end

function similar_cell_field(a::ReducedSkeletonCellBasis,b::CellField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::CellField,b::ReducedSkeletonCellBasis,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::ReducedSkeletonCellBasis,b::CellBasis,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::CellBasis,b::ReducedSkeletonCellBasis,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::ReducedSkeletonCellBasis,b::CellMatrixField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::CellMatrixField,b::ReducedSkeletonCellBasis,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::ReducedSkeletonCellBasis,b::ReducedSkeletonCellBasis,v::AbstractArray)
  @unreachable
end

function operate_cell_field(op,a::ReducedSkeletonCellBasis)
  left = operate_cell_field(op,a.left)
  right = operate_cell_field(op,a.right)
  ReducedSkeletonCellBasis(a.trial_style,left,right)
end

function operate_cell_field(op,a::ReducedSkeletonCellBasis,b::CellField)
  left = operate_cell_field(op,a.left,b)
  right = operate_cell_field(op,a.right,b)
  ReducedSkeletonCellBasis(a.trial_style,left,right)
end

function operate_cell_field(op,a::CellField,b::ReducedSkeletonCellBasis)
  left = operate_cell_field(op,a,b.left)
  right = operate_cell_field(op,a,b.right)
  ReducedSkeletonCellBasis(b.trial_style,left,right)
end

function operate_cell_field(op,a::ReducedSkeletonCellBasis,b::CellBasis)
  @unreachable
end

function operate_cell_field(op,a::CellBasis,b::ReducedSkeletonCellBasis)
  @unreachable
end

function operate_cell_field(op,a::ReducedSkeletonCellBasis,b::CellMatrixField)
  @unreachable
end

function operate_cell_field(op,a::CellMatrixField,b::ReducedSkeletonCellBasis)
  @unreachable
end

function operate_cell_field(op,a::ReducedSkeletonCellBasis,b::ReducedSkeletonCellBasis)
  _operate_reduced_skeleton_cell_basis(op,a,b,a.trial_style,b.trial_style)
end

function  _operate_reduced_skeleton_cell_basis(
  op,a,b,a_trial::Val{T},b_trial::Val{T}) where T
  @notimplemented
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
  ll = operate_cell_field(op,a.left,b.left)
  lr = operate_cell_field(op,a.left,b.right)
  rl = operate_cell_field(op,a.right,b.left)
  rr = operate_cell_field(op,a.right,b.right)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

struct SkeletonCellMatrixField <: CellMatrixField
  ll
  lr
  rl
  rr
end

get_cell_map(a::SkeletonCellMatrixField) = get_cell_map(a.ll)

function get_array(a::SkeletonCellMatrixField)
  @notimplemented
end

function similar_cell_field(cf::SkeletonCellMatrixField,a::AbstractArray)
  @unreachable
end

function similar_cell_field(a::SkeletonCellMatrixField,b::CellField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::CellField,b::SkeletonCellMatrixField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::SkeletonCellMatrixField,b::CellBasis,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::CellBasis,b::SkeletonCellMatrixField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::SkeletonCellMatrixField,b::CellMatrixField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::CellMatrixField,b::SkeletonCellMatrixField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::SkeletonCellMatrixField,b::ReducedSkeletonCellBasis,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::ReducedSkeletonCellBasis,b::SkeletonCellMatrixField,v::AbstractArray)
  @unreachable
end

function similar_cell_field(a::SkeletonCellMatrixField,b::SkeletonCellMatrixField,v::AbstractArray)
  @unreachable
end

function operate_cell_field(op,cf::SkeletonCellMatrixField)
  ll = operate_cell_field(op,cf.ll)
  lr = operate_cell_field(op,cf.lr)
  rl = operate_cell_field(op,cf.rl)
  rr = operate_cell_field(op,cf.rr)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

function operate_cell_field(op,a::SkeletonCellMatrixField,b::CellField)
  ll = operate_cell_field(op,a.ll,b)
  lr = operate_cell_field(op,a.lr,b)
  rl = operate_cell_field(op,a.rl,b)
  rr = operate_cell_field(op,a.rr,b)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

function operate_cell_field(op,a::CellField,b::SkeletonCellMatrixField)
  ll = operate_cell_field(op,a,b.ll)
  lr = operate_cell_field(op,a,b.lr)
  rl = operate_cell_field(op,a,b.rl)
  rr = operate_cell_field(op,a,b.rr)
  SkeletonCellMatrixField(ll,lr,rl,rr)
end

function operate_cell_field(op,a::SkeletonCellMatrixField,b::CellBasis)
  @unreachable
end

function operate_cell_field(op,a::CellBasis,b::SkeletonCellMatrixField)
  @unreachable
end

function operate_cell_field(op,a::SkeletonCellMatrixField,b::CellMatrixField)
  @unreachable
end

function operate_cell_field(op,a::CellMatrixField,b::SkeletonCellMatrixField)
  @unreachable
end

function operate_cell_field(op,a::SkeletonCellMatrixField,b::ReducedSkeletonCellBasis)
  @unreachable
end

function operate_cell_field(op,a::ReducedSkeletonCellBasis,b::SkeletonCellMatrixField)
  @unreachable
end

function operate_cell_field(op,a::SkeletonCellMatrixField,b::SkeletonCellMatrixField)
  ll = operate_cell_field(op,a.ll,b.ll)
  lr = operate_cell_field(op,a.lr,b.lr)
  rl = operate_cell_field(op,a.rl,b.rl)
  rr = operate_cell_field(op,a.rr,b.rr)
  SkeletonCellMatrixField(ll,lr,rl,rr)
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


