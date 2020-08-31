module KK

using Gridap
using Gridap.CellData
using Gridap.Arrays
using Gridap.Helpers
using Gridap.Geometry
using Gridap.Fields
using Gridap.FESpaces

function get_new_cell_map(trian::Triangulation)
  array = get_cell_map(trian)
  GenericCellField(array,array)
end

function get_new_normal_vector(trian::GenericBoundaryTriangulation)
  k = Geometry.NormalVectorValued()
  refn = Geometry.ReferenceNormal(trian)
  cell_map = restrict(get_cell_map(trian.cell_trian),trian)
  J = gradient(cell_map)
  a = apply(k,J,refn)
  GenericCellField(a,get_cell_map(trian))
end

function get_new_normal_vector(trian::SkeletonTriangulation)
  get_new_normal_vector(trian.left)
end

function get_new_cell_map(btrian::BoundaryTriangulation)
  trian = get_volume_triangulation(btrian)
  cell_map = get_cell_map(trian)
  face_map = get_cell_map(btrian)
  face_to_cell = get_face_to_cell(btrian)
  refface_to_refcell_map = get_face_to_cell_map(btrian)
  cell_field = GenericCellField(face_map,face_map)
  FaceMap(
    cell_field,
    cell_map,
    face_to_cell,
    refface_to_refcell_map)
end

# TODO not sure if we want to dispatch on SkeletonTriangulation which is a concrete type
# Perhaps we can create an abstract type like for the boundary?
function get_new_cell_map(strian::SkeletonTriangulation)
  left = get_new_cell_map(strian.left)
  _right = get_new_cell_map(strian.right)
  right = change_face_map(_right,left.face_map)
  SkeletonFaceMap(left,right)
end

struct FaceMap{T<:CellField} <: CellField
  face_map::T
  cell_map::AbstractArray
  face_to_cell::AbstractArray
  refface_to_refcell_map::AbstractArray
end

function change_face_map(a::FaceMap,face_map::CellField)
  FaceMap(
    face_map,
    a.cell_map,
    a.face_to_cell,
    a.refface_to_refcell_map)
end

Arrays.get_array(a::FaceMap) = get_array(a.face_map)
CellData.get_cell_map(a::FaceMap) = get_cell_map(a.face_map)
CellData.get_cell_axes(a::FaceMap) = get_cell_axes(a.face_map)
CellData.get_memo(a::FaceMap) = get_memo(a.face_map)
CellData.RefStyle(::Type{FaceMap{T}}) where T = RefStyle(T)
CellData.MetaSizeStyle(::Type{FaceMap{T}}) where T = MetaSizeStyle(T)

struct SkeletonFaceMap{L<:CellField} <:CellField
  left::L
  right::FaceMap
end

Arrays.get_array(a::SkeletonFaceMap) = get_array(a.left)
CellData.get_cell_map(a::SkeletonFaceMap) = get_cell_map(a.left)
CellData.get_cell_axes(a::SkeletonFaceMap) = get_cell_axes(a.left)
CellData.get_memo(a::SkeletonFaceMap) = get_memo(a.left)
CellData.RefStyle(::Type{SkeletonFaceMap{T}}) where T = RefStyle(T)
CellData.MetaSizeStyle(::Type{SkeletonFaceMap{T}}) where T = MetaSizeStyle(T)

#TODO equalty of cell maps is the crucial ingredient, refine this part (overload == with efficient code)

function Base.:∘(f::CellField,ϕ::CellField)
  if get_cell_map(f) === get_array(ϕ) || get_cell_map(f) == get_array(ϕ)
    to_ref_space(f)
  else
    @notimplemented
  end
end

function Base.:∘(object,ϕ::CellField)
  convert_to_cell_field(object,get_array(ϕ))
end

function Base.:∘(f::CellField,ϕ::FaceMap)
  if get_cell_map(f) === get_array(ϕ) || get_cell_map(f) == get_array(ϕ) || typeof(get_cell_map(f)) == typeof(get_array(ϕ)) # TODO the last part is an inadmissible hack
    to_ref_space(f)
  elseif get_cell_map(f) === ϕ.cell_map || get_cell_map(f) == ϕ.cell_map
    _to_ref_face_space(f,ϕ.face_to_cell,ϕ.refface_to_refcell_map,get_array(ϕ.face_map))
  else
    @notimplemented
  end
end

function _to_ref_face_space(f,face_to_cell,refface_to_refcell_map,face_map)
  # This function is not efficient for f in the physical space
  cell_to_f = to_ref_space(f)
  face_to_f = reindex(cell_to_f,face_to_cell)
  array = compose_field_arrays(get_array(face_to_f),refface_to_refcell_map)
  axs = reindex(get_cell_axes(f),face_to_cell)
  similar_object(f,array,face_map,axs,MetaSizeStyle(f))
end

function Base.:∘(f::CellField,ϕ::SkeletonFaceMap)
  f∘ϕ.left
end

struct CellFieldFromOperation{F} <: GridapType
  op::F
  args
end

function Base.:∘(f::CellFieldFromOperation,ϕ::CellField)
  f.op(map(i->i∘ϕ,f.args)...)
end

function Base.:∘(f::CellFieldFromOperation{typeof(get_inward)},ϕ::SkeletonFaceMap)
  scf = merge_cell_fields_at_skeleton(f.args[1]∘ϕ.left,f.args[1]∘ϕ.right)
  get_inward(scf)
end

function Base.:∘(f::CellFieldFromOperation{typeof(get_outward)},ϕ::SkeletonFaceMap)
  scf = merge_cell_fields_at_skeleton(f.args[1]∘ϕ.left,f.args[1]∘ϕ.right)
  get_outward(scf)
end

for op in (:get_inward, :get_outward, :+ ,:-)
  @eval begin

    function Fields.gradient(f::CellFieldFromOperation{typeof($op)})
      a = gradient(f.args[1])
      CellFieldFromOperation($op,(a,))
    end

  end
end

## Just a temporary wrapper to change the behavior of operate()
# Do not include in Gridap
struct TmpCellField{T<:CellField} <: CellField
  cf::T
end

for op in (:+ ,:-)
  @eval begin

    function Fields.gradient(f::CellFieldFromOperation{typeof($op)})
      a = gradient(f.args[1])
      b = gradient(f.args[2])
      CellFieldFromOperation($op,(a,b))
    end

    function Helpers.operate(::typeof($op),a::Union{TmpCellField,CellFieldFromOperation}...)
      CellFieldFromOperation($op,a)
    end

    function Helpers.operate(::typeof($op),a::Union{TmpCellField,CellFieldFromOperation},b)
      CellFieldFromOperation($op,(a,b))
    end

    function Helpers.operate(::typeof($op),a,b::Union{TmpCellField,CellFieldFromOperation})
      CellFieldFromOperation($op,(a,b))
    end

  end
end

Arrays.get_array(a::TmpCellField) = get_array(a.cf)
CellData.get_cell_map(a::TmpCellField) = get_cell_map(a.cf)
CellData.get_cell_axes(a::TmpCellField) = get_cell_axes(a.cf)
CellData.get_memo(a::TmpCellField) = get_memo(a.cf)
CellData.RefStyle(::Type{TmpCellField{T}}) where T = RefStyle(T)
CellData.MetaSizeStyle(::Type{TmpCellField{T}}) where T = MetaSizeStyle(T)

function Helpers.operate(op,a::Union{TmpCellField,CellFieldFromOperation}...)
  f(args...) = operate(op,args...)
  CellFieldFromOperation(f,a)
end

function Helpers.operate(op,a::Union{TmpCellField,CellFieldFromOperation},b)
  f(a,b) = operate(op,a,b)
  CellFieldFromOperation(f,(a,b))
end

function Helpers.operate(op,a,b::Union{TmpCellField,CellFieldFromOperation})
  f(a,b) = operate(op,a,b)
  CellFieldFromOperation(f,(a,b))
end

function CellData.get_inward(a::Union{TmpCellField,CellFieldFromOperation})
  CellFieldFromOperation(get_inward,(a,))
end

function CellData.get_outward(a::Union{TmpCellField,CellFieldFromOperation})
  CellFieldFromOperation(get_outward,(a,))
end

function Base.getproperty(x::Union{TmpCellField,CellFieldFromOperation}, sym::Symbol)
  if sym in (:inward,:⁺)
    get_inward(x)
  elseif sym in (:outward,:⁻)
    get_outward(x)
  else
    getfield(x, sym)
  end
end

function Fields.gradient(a::TmpCellField)
  TmpCellField(gradient(a.cf))
end

function Base.:∘(f::TmpCellField,ϕ::CellField)
  f.cf∘ϕ
end

function Base.:∘(f::TmpCellField,ϕ::FaceMap)
  f.cf∘ϕ
end

function Base.:∘(f::TmpCellField,ϕ::SkeletonFaceMap)
  f.cf∘ϕ
end

# TODO use this definition in Gridap
function jump(sf)
  sf.⁺ - sf.⁻
end

# TODO use this definition in Gridap
function mean(sf)
  operate(CellData._mean,sf.⁺,sf.⁻)
end

struct LebesgueMeasure <: GridapType
  trian::Triangulation
  quad::CellQuadrature
end

function LebesgueMeasure(trian::Triangulation,degree::Integer)
  quad = CellQuadrature(trian,degree)
  LebesgueMeasure(trian,quad)
end

struct Integral
  integrand
end

const ∫ = Integral

function Base.:*(a::Integral,b::LebesgueMeasure)
  f = a.integrand
  ϕ = get_new_cell_map(b.trian)
  cell_value = integrate( f∘ϕ, b.trian, b.quad)
  cell_id = get_cell_id(b.trian)
  CellContribution(cell_value,cell_id,b)
end

struct CellContribution{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  cell_values::A
  cell_id
  object
end

Base.size(a::CellContribution) = size(get_array(a))
Arrays.array_cache(a::CellContribution) = array_cache(get_array(a))
@inline Arrays.getindex!(cache,a::CellContribution,i::Integer...) = getindex!(cache,get_array(a),i...)
@inline Base.getindex(a::CellContribution,i::Integer...) = get_array(a)[i...]
Base.IndexStyle(::Type{CellContribution{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.sum(a::CellContribution) = sum(get_array(a))

Arrays.get_array(a::CellContribution) = a.cell_values
Geometry.get_cell_id(a::CellContribution) = a.cell_id
get_object(a::CellContribution) = a.object

struct CollectionOfCellContribution <: GridapType
  dict::Dict{UInt,CellContribution}
end

function get_cell_contribution(a::CollectionOfCellContribution,object)
  id = objectid(object)
  a.dict[id]
end

function Base.:+(a::CellContribution,b::CellContribution)
  if get_object(a) === get_object(b)
    array = apply(bcast(+),get_array(a),get_array(b))
    CellContribution(array,get_cell_id(a),get_object(a))
  else
    dict = Dict{UInt,CellContribution}()
    dict[objectid(get_object(a))] = a
    dict[objectid(get_object(b))] = b
    CollectionOfCellContribution(dict)
  end
end

function Base.:+(a::CollectionOfCellContribution,b::CellContribution)
  id = objectid(get_object(b))
  if haskey(a.dict,id)
    c = a.dict[id]
    a.dict[id] = c + b
  else
    a.dict[id] = b
  end
  a
end

function Base.:+(a::CellContribution,b::CollectionOfCellContribution)
  b + a
end

# TODO this function will be in CellData
function FESpaces.collect_cell_matrix(a::CollectionOfCellContribution)
  w = []
  r = []
  for cell_contribution in values(a.dict)
    push!(w,get_array(cell_contribution))
    push!(r,get_cell_id(cell_contribution))
  end
  (w,r,r)
end

function FESpaces.collect_cell_matrix(a::CellContribution)
  w = [get_array(a)]
  r = [get_cell_id(a)]
  (w,r,r)
end

# TODO this function will be in CellData
function FESpaces.collect_cell_vector(a::CollectionOfCellContribution)
  w = []
  r = []
  for cell_contribution in values(a.dict)
    push!(w,get_array(cell_contribution))
    push!(r,get_cell_id(cell_contribution))
  end
  (w,r)
end

function FESpaces.collect_cell_vector(a::CellContribution)
  w = [get_array(a)]
  r = [get_cell_id(a)]
  (w,r)
end

# TODO this function will be in CellData
function FESpaces.collect_cell_matrix_and_vector(
  biform::CollectionOfCellContribution,liform::CollectionOfCellContribution)

  matvec, mat, vec = _pair_contribution_when_possible(biform,liform)
  matvecdata = collect_cell_matrix(matvec)
  matdata = collect_cell_matrix(mat)
  vecdata = collect_cell_matrix(mat)
  (matvecdata, matdata, vecdata)
end

function FESpaces.collect_cell_matrix_and_vector(
  biform::CollectionOfCellContribution,liform::CollectionOfCellContribution,uhd)

  @assert is_a_fe_function(uhd)
  matvec, mat, vec = _pair_contribution_when_possible(biform,liform,uhd)
  matvecdata = collect_cell_matrix(matvec)
  matdata = collect_cell_matrix(mat)
  vecdata = collect_cell_matrix(mat)
  (matvecdata, matdata, vecdata)
end

function _pair_contribution_when_possible(biform,liform)
  matvec = Dict{UInt,CellContribution}()
  mat = Dict{UInt,CellContribution}()
  vec = Dict{UInt,CellContribution}()
  for (id,t) in biform.dict
    if haskey(liform.dict,id)
      a = pair_arrays(get_array(biform.dict[id]),get_array(liform.dict[id]))
      matvec[id] = CellContribution(a,t.cell_id,t.object)
    else
      mat[id] = t
    end
  end
  for (id,t) in liform.dict
    if ! haskey(biform.dict,id)
      vec[id] = t
    end
  end
  map(CollectionOfCellContribution, (matvec, mat, vec))
end

function _pair_contribution_when_possible(biform,liform,uhd)
  matvec = Dict{UInt,CellContribution}()
  mat = Dict{UInt,CellContribution}()
  vec = Dict{UInt,CellContribution}()
  _matvec, _mat, _vec = _pair_contribution_when_possible(biform,liform)
  for (id,t) in _matvec.dict
    cellvals = get_cell_values(uhd,t.cell_id)
    a = attach_dirichlet(get_array(t),cellvals)
    matvec[id] = CellContribution(a,t.cell_id,t.object)
  end
  for (id,t) in _mat.dict
    cellvals = get_cell_values(uhd,t.cell_id)
    a = attach_dirichlet(get_array(t),cellvals)
    matvec[id] = CellContribution(a,t.cell_id,t.object)
  end
  map(CollectionOfCellContribution, (matvec, mat, vec))
end

#u(x) = x[1]^2 + x[2]
u(x) = x[2] + x[1]
ud = u
f(x) = Δ(u)(x)

domain = (0,1,0,1)
cells = (10,10)
model = CartesianDiscreteModel(domain,cells)

trian_Ω = Triangulation(model)
trian_Γ = BoundaryTriangulation(model)
trian_Λ = SkeletonTriangulation(model)

ϕ_Ω = get_new_cell_map(trian_Ω)
ϕ_Γ = get_new_cell_map(trian_Γ)
ϕ_Λ = get_new_cell_map(trian_Λ)

Vh = FESpace(model=model,reffe=:Lagrangian,order=1,valuetype=Float64,conformity=:H1)
Uh = TrialFESpace(Vh)
vh = TmpCellField(FEFunction(Vh,rand(num_free_dofs(Vh))))
uh = TmpCellField(FEFunction(Uh,rand(num_free_dofs(Uh))))
dvh = TmpCellField(get_cell_basis(Vh))
duh = TmpCellField(get_cell_basis(Uh))
n_Γ = TmpCellField(get_new_normal_vector(trian_Γ))
n_Λ = TmpCellField(get_new_normal_vector(trian_Λ))

vh ∘ ϕ_Ω
vh ∘ ϕ_Γ
n_Γ ∘ ϕ_Γ
n_Λ ∘ ϕ_Λ

uh.⁺∘ ϕ_Λ
uh.⁻∘ ϕ_Λ

cfop = (uh.⁺ - uh.⁻) ∘ ϕ_Λ
#@show typeof(cfop)
#@show typeof(cfop.args[1].args[1])
#@show typeof(cfop.args[2].args[1])

cfop = (uh*n_Λ).⁺ ∘ ϕ_Λ

( vh*(n_Γ⋅∇(uh)) )∘ϕ_Γ
writevtk(trian_Ω,"trian_Ω",cellfields=["v"=>vh∘ϕ_Ω,"vu"=>( ∇(vh)⋅∇(uh) )∘ϕ_Ω])

writevtk(trian_Γ,"trian_Γ",cellfields=
  ["v"=>vh∘ϕ_Γ,
   "n_Γ"=>n_Γ∘ϕ_Γ,
   "n_Γ⋅∇v"=>(n_Γ∘ϕ_Γ)⋅(∇(vh)∘ϕ_Γ),
   "v*n_Γ⋅∇u"=>( vh*(n_Γ⋅∇(uh)) )∘ϕ_Γ])

writevtk(trian_Λ,"trian_Λ",nsubcells=3,cellfields= [
  "n_Λ"=>n_Λ∘ϕ_Λ,
  "u1"=>uh.⁺∘ ϕ_Λ,
  "u2"=>uh.⁻∘ ϕ_Λ,
  "j"=>jump(uh)∘ϕ_Λ,
  "gu+" => ∇(uh.⁺) ∘ ϕ_Λ,
  "gu-" => ∇(uh.⁻) ∘ ϕ_Λ,
 ])

degree = 3
dΩ = LebesgueMeasure(trian_Ω,degree)
dΓ = LebesgueMeasure(trian_Γ,degree)
dΛ = LebesgueMeasure(trian_Λ,degree)

a_Ω(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
a_Γ(u,v) = ∫( v*(n_Γ⋅∇(u)) )*dΓ

const h = .25
const γ = 10

a_Λ(u,v) = ∫( (γ/h)*jump(v*n_Λ)⋅jump(u*n_Λ) )*dΛ

a(u,v) = a_Ω(u,v) + a_Γ(u,v) + a_Λ(u,v)

cellval = a_Λ(uh,vh)
cellvec = a_Λ(uh,dvh)
cellmat = a_Λ(duh,dvh)

#display(cellmat[1])

cellval = a_Ω(uh,vh)
cellvec = a_Ω(uh,dvh)
cellmat = a_Ω(duh,dvh)

cellval = a(uh,vh)
cellvec = a(uh,dvh)
cellmat = a(duh,dvh)

#display( a_Ω(uh,vh) + a_Ω(uh,vh) )
#display( a(uh,vh) )
#display(cellvec)

#(n_Γ∘ϕ_Γ)⋅(∇(vh)∘ϕ_Γ)
#
#cf = o(*,o(⋅,n_Γ,∇(uh)),vh)∘ϕ_Γ
#@show typeof(cf)
#
##(n_Γ⋅∇(vh))∘ϕ_Γ
#
#writevtk(trian_Ω,"trian_Ω",cellfields=["v"=>vh∘ϕ_Ω])
#writevtk(trian_Γ,"trian_Γ",cellfields=["v"=>vh∘ϕ_Γ,"n_Γ"=>n_Γ∘ϕ_Γ,"n_Γ⋅∇v"=>(n_Γ∘ϕ_Γ)⋅(∇(vh)∘ϕ_Γ)])

a(u,v) =
  ∫( ∇(v)⋅∇(u) )*dΩ +
  ∫( (γ/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )*dΓ +
  ∫(
    (γ/h)*jump(v*n_Λ)⋅jump(u*n_Λ) -
    jump(v*n_Λ)⋅mean(∇(u)) -
    mean(∇(v))⋅jump(u*n_Λ) )*dΛ

biform = a(duh,dvh)

display( get_cell_contribution(biform,dΩ)[1] )
display( get_cell_contribution(biform,dΓ)[1] )
display( get_cell_contribution(biform,dΛ)[1] )

l(v) =
  ∫( v*f )*dΩ + 
  ∫( (γ/h)*v*ud - (n_Γ⋅∇(v))*ud )*dΓ

liform = l(dvh)
display( get_cell_contribution(liform,dΩ)[1] )
display( get_cell_contribution(liform,dΓ)[1] )

assem = SparseMatrixAssembler(Vh,Uh)
uhd = zero(Uh)
data = collect_cell_matrix_and_vector(biform,liform,uhd)
A,b = assemble_matrix_and_vector(assem,data)

x = A\b 
uh = TmpCellField(FEFunction(Uh,x))

e = u - uh

e_l2 = sqrt(sum( ∫( e*e )*dΩ ))
e_h1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))
@show e_l2
@show e_h1

e_l2 = sqrt(sum( ∫( e*e )*dΓ ))
e_h1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΓ ))
@show e_l2
@show e_h1

e_l2 = sqrt(sum( ∫( mean(e)*mean(e) )*dΛ ))
e_h1 = sqrt(sum( ∫( mean(e)*mean(e) + mean(∇(e))⋅mean(∇(e)) )*dΛ ))
@show e_l2
@show e_h1





end # mdoule
