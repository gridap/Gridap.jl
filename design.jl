module KK

using Gridap
using Gridap.CellData
using Gridap.Arrays
using Gridap.Helpers
using Gridap.Geometry
using Gridap.Fields

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

struct FaceMap{T<:CellField} <: CellField
  face_map::T
  cell_map::AbstractArray
  face_to_cell::AbstractArray
  refface_to_refcell_map::AbstractArray
end

Arrays.get_array(a::FaceMap) = get_array(a.face_map)
CellData.get_cell_map(a::FaceMap) = get_cell_map(a.face_map)
CellData.get_cell_axes(a::FaceMap) = get_cell_axes(a.face_map)
CellData.get_memo(a::FaceMap) = get_memo(a.face_map)
CellData.RefStyle(::Type{FaceMap{T}}) where T = RefStyle(T)
CellData.MetaSizeStyle(::Type{FaceMap{T}}) where T = MetaSizeStyle(T)

#TODO equalty of cell maps is the crucial ingredient, refine this part (overload == with efficient code)

function Base.:∘(f::CellField,ϕ::CellField)
  if get_cell_map(f) === get_array(ϕ) || get_cell_map(f) == get_array(ϕ)
    to_ref_space(f)
  else
    @notimplemented
  end
end

function Base.:∘(f::CellField,ϕ::FaceMap)
  if get_cell_map(f) === get_array(ϕ) || get_cell_map(f) == get_array(ϕ) || typeof(get_cell_map(f)) == typeof(get_array(ϕ)) # TODO the last part is an inadmissible hack
    to_ref_space(f)
  elseif get_cell_map(f) === ϕ.cell_map || get_cell_map(f) == ϕ.cell_map
    _to_ref_face_space(f,ϕ.face_to_cell,ϕ.refface_to_refcell_map,get_array(ϕ.face_map))
  else
    @show typeof(get_cell_map(f))
    @show typeof(get_array(ϕ))
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

struct CellFieldFromOperation
  op
  args
end

function o(op,args...)
  CellFieldFromOperation(op,args)
end

function Base.:∘(f::CellFieldFromOperation,ϕ::CellField)
  operate(f.op,map(i->i∘ϕ,f.args)...)
end

## Just a temporary wrapper to change the behavior of operate()
struct TmpCellField{T<:CellField} <: CellField
  cf::T
end

Arrays.get_array(a::TmpCellField) = get_array(a.cf)
CellData.get_cell_map(a::TmpCellField) = get_cell_map(a.cf)
CellData.get_cell_axes(a::TmpCellField) = get_cell_axes(a.cf)
CellData.get_memo(a::TmpCellField) = get_memo(a.cf)
CellData.RefStyle(::Type{TmpCellField{T}}) where T = RefStyle(T)
CellData.MetaSizeStyle(::Type{TmpCellField{T}}) where T = MetaSizeStyle(T)

function Helpers.operate(op,a::Union{TmpCellField,CellFieldFromOperation}...)
  CellFieldFromOperation(op,a)
end

function Fields.gradient(a::TmpCellField)
  TmpCellField(gradient(a.cf))
end

domain = (0,1,0,1)
cells = (10,10)
model = CartesianDiscreteModel(domain,cells)

trian_Ω = Triangulation(model)
trian_Γ = BoundaryTriangulation(model)

ϕ_Ω = get_new_cell_map(trian_Ω)
ϕ_Γ = get_new_cell_map(trian_Γ)

Vh = FESpace(model=model,reffe=:Lagrangian,order=1,valuetype=Float64,conformity=:H1)
Uh = TrialFESpace(Vh)
vh = TmpCellField(FEFunction(Vh,rand(num_free_dofs(Vh))))
uh = TmpCellField(FEFunction(Uh,rand(num_free_dofs(Uh))))
n = TmpCellField(get_new_normal_vector(trian_Γ))

vh ∘ ϕ_Ω
vh ∘ ϕ_Γ
n ∘ ϕ_Γ

( vh*(n⋅∇(uh)) )∘ϕ_Γ

writevtk(trian_Ω,"trian_Ω",cellfields=["v"=>vh∘ϕ_Ω])
writevtk(trian_Γ,"trian_Γ",cellfields=
  ["v"=>vh∘ϕ_Γ,
   "n"=>n∘ϕ_Γ,
   "n⋅∇v"=>(n∘ϕ_Γ)⋅(∇(vh)∘ϕ_Γ),
   "v*n⋅∇u"=>( vh*(n⋅∇(uh)) )∘ϕ_Γ])

#(n∘ϕ_Γ)⋅(∇(vh)∘ϕ_Γ)
#
#cf = o(*,o(⋅,n,∇(uh)),vh)∘ϕ_Γ
#@show typeof(cf)
#
##(n⋅∇(vh))∘ϕ_Γ
#
#writevtk(trian_Ω,"trian_Ω",cellfields=["v"=>vh∘ϕ_Ω])
#writevtk(trian_Γ,"trian_Γ",cellfields=["v"=>vh∘ϕ_Γ,"n"=>n∘ϕ_Γ,"n⋅∇v"=>(n∘ϕ_Γ)⋅(∇(vh)∘ϕ_Γ)])



end # mdoule
