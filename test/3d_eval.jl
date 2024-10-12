
using Gridap
using Gridap.FESpaces, Gridap.CellData, Gridap.Fields, Gridap.Arrays

using FillArrays

# Function that goes from 2D points to 2-component vectors
u2D(x::VectorValue{2}) = VectorValue(x[1],x[2])

function point_3D_to_2D(x::VectorValue{3})
  return VectorValue(x[1],x[2])
end

function result_2D_to_3D(v::VectorValue{2})
  return VectorValue(v[1],v[2],0.0)
end

function embed_3D(u2D::CellField,trian_3D::Triangulation)
  u22_data = CellData.get_data(u2D)
  ncells = length(u22_data)

  u32_data = lazy_map(Broadcasting(∘),u22_data,Fill(GenericField(point_3D_to_2D),ncells))
  u33_data = lazy_map(Broadcasting(∘),Fill(GenericField(result_2D_to_3D),ncells),u32_data)
  u3D = CellData.similar_cell_field(u2D,u33_data,trian_3D,PhysicalDomain())
  return u3D
end

model_2D = CartesianDiscreteModel((0,1,0,1),(2,2))
model_3D = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))

reffe_2D = ReferenceFE(lagrangian,VectorValue{2,Float64},1)
reffe_3D = ReferenceFE(lagrangian,VectorValue{3,Float64},1)

V_2D = TestFESpace(model_2D,reffe_2D)
V_3D = TestFESpace(model_3D,reffe_3D)

uh_2D = interpolate(u2D,V_2D)
uh_3D = embed_3D(uh_2D,get_triangulation(V_3D))

pts = get_cell_points(get_triangulation(V_3D))
uh_3D(pts)
