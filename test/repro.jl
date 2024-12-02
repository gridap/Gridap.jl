module reproduction_file
using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.Fields
using GridapEmbedded
using Gridap.Arrays
using Gridap.FESpaces

model = CartesianDiscreteModel(Point(0.0, 0.0), Point(1.0, 1.0), (5,5))
Ω = Interior(model)
dΩ = Measure(Ω,2)
geo = disk(0.1, x0=Point(0.5,0.5))
cutgeo = cut(model, !geo)

W1 = FESpace(Ω, ReferenceFE(lagrangian, Float64, 1), vector_type=Vector{ComplexF64})
@show get_dof_value_type(W1)

W2 = AgFEMSpace(W1,aggregate(AggregateCutCellsByThreshold(1.0), cutgeo, geo, OUT))  # AgFEMSpace returns a FESpaceWithLinearConstraints
@show get_dof_value_type(W2)
Φ = TrialFESpace(W2)
@show get_dof_value_type(Φ)

a_wϕ(ϕ, w) = ∫( im*∇(ϕ)⋅∇(w) )dΩ        # if we do NOT include the imaginary operator im here, then the code will pass, but the FE Spaces will be of type Float64
A_wϕ1 = assemble_matrix(a_wϕ, Φ, W1)    # function will NOT fail here, due to get_dof_value_type function being implemented
A_wϕ2 = assemble_matrix(a_wϕ, Φ, W2)    # function will fail here, because W2 should be ComplexF64, but is Float64 due to missing get_dof_value_type function

end # module