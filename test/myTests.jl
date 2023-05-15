
using Gridap
using Gridap.Geometry
using Gridap.Adaptivity
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.CellData

using FillArrays

############################################################################################
# Fine model
domain = (0,1,0,1,0,1)
partition = (2,2,3)
model = CartesianDiscreteModel(domain,partition)

Ω  = Triangulation(model)
Γ  = BoundaryTriangulation(model,tags=["tag_26"]) # y-z, x = 1

############################################################################################
# Coarse model

Dc      = num_cell_dims(model)
topo    = get_grid_topology(model)
f2n_map = Geometry.get_faces(topo,Dc-1,0)
coords  = Geometry.get_vertex_coordinates(topo)
face_coords = map(Reindex(coords),f2n_map)

tol = 1.e-3
com_coords = lazy_map(N->sum(N)/length(N),face_coords) # center of mass
int_coords = map(N->VectorValue(Int(floor(N[1]/tol)),Int(floor(N[2]/tol)),Int(floor(N[3]/tol))),com_coords)

face_x = maximum(lazy_map(c->c[1],int_coords))
interface_faces = findall(c->c[1]==face_x,int_coords)
interface_coords = view(int_coords,interface_faces)

y_coords = map(c->c[2],interface_coords)
y_unique = sort(unique(map(c->c[2],interface_coords)))
y_counts = [count(==(y),y_coords) for y in sort(unique(y_coords))]
y_ptrs = Gridap.Adaptivity.counts_to_ptrs(y_counts)

perm = sortperm(interface_coords,by=x->x[2])
data = lazy_map(Reindex(interface_faces),perm)
c2f_faces = Table(data,y_ptrs)

n2o_cells = zeros(Int, length(Γ.glue.face_to_bgface))
child_ids = zeros(Int, length(Γ.glue.face_to_bgface))
for (islide, slide_list) in enumerate(c2f_faces)
    for (izpos, num) in enumerate(slide_list)
        ind = findall(x->x==num, Γ.glue.face_to_bgface)[1]
        n2o_cells[ind] = islide
        child_ids[ind] = izpos
    end  
end

n2o_faces = [Int[],Int[],n2o_cells]
rrules = Fill(RefinementRule(QUAD,(1,length(c2f_faces[1]))),length(c2f_faces))
glue = AdaptivityGlue(n2o_faces,child_ids,rrules) # From coarse to fine

cface_model = CartesianDiscreteModel((0,1,0,1),(length(c2f_faces),1))

Γc = Triangulation(cface_model)
Γf = Adaptivity.GluedTriangulation(Γ,Γc,glue)

############################################################################################
# FESpaces 

reffe_u = ReferenceFE(lagrangian,Float64,1)
reffe_λ = ReferenceFE(lagrangian,Float64,0)

Vu = FESpace(Ω,reffe_u)
Uu = TrialFESpace(Vu)

Vλ = FESpace(Γc,reffe_λ,conformity=:L2)
Uλ = TrialFESpace(Vλ)

Y = MultiFieldFESpace([Vu,Vλ])
X = MultiFieldFESpace([Uu,Uλ])

dΩ = Measure(Ω, 3)
dΓ = Measure(Γf, 3)

tr_Γf(λ) = change_domain(λ,Γf,DomainStyle(λ))

α = 1.e-2
aΩ((u,λ),(v,μ)) = ∫(u*v)*dΩ
aΓ((u,λ),(v,μ)) = ∫(tr_Γf(λ*μ) - α*(u*tr_Γf(μ) + v*tr_Γf(λ)))*dΓ
a((u,λ),(v,μ))  = aΩ((u,λ),(v,μ)) + aΓ((u,λ),(v,μ))

A = assemble_matrix(a,X,Y)

