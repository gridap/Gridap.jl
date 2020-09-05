module AttachConstraintsTests

using Test
using FillArrays
using LinearAlgebra
using Gridap.Helpers
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Integration
using Gridap.CellData
using Gridap.TensorValues

using Gridap.Fields: MockField, MockBasis, OtherMockBasis

ncells = 10
ndofs = 3
ndofs_c = 4
cellmat = [ rand(ndofs,ndofs) for cell in 1:ncells]
cellvec = [ rand(ndofs) for cell in 1:ncells]
cellmatvec = pair_arrays(cellmat,cellvec)
cellconstr = [ rand(ndofs_c,ndofs) for cell in 1:ncells]

a = attach_constraints_rows(cellvec,cellconstr)
r = map( (vec,constr) -> constr*vec ,cellvec,cellconstr)
test_array(a,r)

a = attach_constraints_rows(cellmat,cellconstr)
r = map( (mat,constr) -> constr*mat ,cellmat,cellconstr)
test_array(a,r)

a = attach_constraints_rows(cellmatvec,cellconstr)
r = map( (mat,vec,constr) -> (constr*mat,constr*vec),cellmat,cellvec,cellconstr)
test_array(a,r)

a = attach_constraints_cols(cellmat,cellconstr)
r = map( (mat,constr) -> mat*transpose(constr) ,cellmat,cellconstr)
test_array(a,r)

a = attach_constraints_cols(cellmatvec,cellconstr)
r = map( (mat,vec,constr) -> (mat*transpose(constr),vec),cellmat,cellvec,cellconstr)
test_array(a,r)

#cache = array_cache(a)
#using BenchmarkTools
#@btime getindex!($cache,$a,2)

ndofs = 4

z = 2.0
v = VectorValue(3.0,1.5)
w = VectorValue(3.4,3.5)
a = MockBasis{2}(v,ndofs)
b = MockBasis{2}(w,ndofs)
c = fill(1.0,ndofs)
f = OtherMockBasis{2}(ndofs)
g = MockField{2}(v)

l = 10
zl = [ z for  i in 1:l]
cl = fill(c,l)
fl = Fill(f,l)
ϕ = GenericCellField(lincomb(fl,cl))
gl = fill(g,l)
al = Fill(a,l)
bl = fill(b,l)

cell_axes = Fill((Base.OneTo(ndofs),),l)
gf = GenericCellField(gl)∘inverse_map(ϕ)
af = GenericCellField(al,cell_axes,Val((:,)))∘inverse_map(ϕ)
bf = GenericCellField(bl,cell_axes,Val((:,)))∘inverse_map(ϕ)
#zf = convert_to_cell_field(zl,ϕ)∘inverse_map(ϕ)
zf = convert_to_cell_field(zl,length(ϕ))
df = af*zf
dft = trialize_cell_basis(df)

degree = 3
quad_ref = CellQuadrature(degree,[QUAD,],ones(Int,l))
quad = ϕ(quad_ref)

q = get_coordinates(quad_ref)

r1 = (∇(2*af)∘ϕ)(q)
r2 = (2*∇(af)∘ϕ)(q)
test_array(r1,r2)

r1 = (∇(af*2)∘ϕ)(q)
r2 = (2*∇(af)∘ϕ)(q)
test_array(r1,r2)

ndofs_c = 5
cellconstr = [ rand(ndofs_c,ndofs) for i in 1:l]
cellvals = [ rand(ndofs) for i in 1:l]
cellvec = integrate(bf⋅v,quad)
cellmat = integrate(∇(bf)⊙∇(dft),quad)
cellmatvec = pair_arrays(cellmat,cellvec)
cellmatvec = attach_dirichlet(cellmatvec,cellvals)
cellmatvec = attach_constraints_rows(cellmatvec,cellconstr)
cellmatvec = attach_constraints_cols(cellmatvec,cellconstr)

@test size(cellmatvec[1][1]) == (ndofs_c,ndofs_c)
@test size(cellmatvec[1][2]) == (ndofs_c,)

cellconstr = identity_constraints(cell_axes)
test_array(cellconstr,Fill(Matrix(I,ndofs,ndofs),l))
@test isa(cellconstr,Fill)

# Test at skeleton


axesL = Fill((Base.OneTo(ndofs),),l)
axesR = axesL
cellvalsL = [ rand(ndofs) for i in 1:l]
cellvalsR = cellvalsL
cellvals = merge_cell_dofs_at_skeleton(cellvalsL,cellvalsR,axesL,axesR)
cellconstrL = [ rand(ndofs_c,ndofs) for i in 1:l]
cellconstrR = cellconstrL
axesL_rows = Fill((Base.OneTo(ndofs_c),),l)
axesR_rows = axesL_rows
axesL_cols = axesL
axesR_cols = axesR
cellconstr = merge_cell_constraints_at_skeleton(cellconstrL,cellconstrR,axesL_rows,axesR_rows,axesL_cols,axesR_cols)

ϕ_Γ = SkeletonFaceMap(ϕ,ϕ)
quad_Γ = ϕ_Γ(quad_ref)

cellvec = integrate(jump(af⋅v), quad_Γ)
cellmat = integrate( jump(af⋅v)*(w⋅dft.⁻),quad_Γ)

cellmatvec = pair_arrays(cellmat,cellvec)
cellmatvec = attach_dirichlet(cellmatvec,cellvals)
cellmatvec = attach_constraints_rows(cellmatvec,cellconstr)
cellmatvec = attach_constraints_cols(cellmatvec,cellconstr)
test_array(cellmatvec,collect(cellmatvec))

cellconstr = identity_constraints(cell_axes)
cellconstr = merge_cell_constraints_at_skeleton(cellconstr,cellconstr,cell_axes,cell_axes,cell_axes,cell_axes)
@test isa(cellconstr,VectorOfBlockArrayCoo)

##a = cellmatvec
#a = cellconstr
#cache = array_cache(a)
#using BenchmarkTools
#@btime getindex!($cache,$a,2)

end # module
