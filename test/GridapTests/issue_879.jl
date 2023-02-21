module Issue873

using Gridap

# Create a Cartesian grid
domain =  (0, 1, 0, 0.5)
partition = (10, 5)
model = CartesianDiscreteModel(domain, partition)

# Assign the label "surface" to the entity 3,4 and 6 
# (top corners and top side)
labels_Ω = get_face_labeling(model)
add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   

# Triangulations
Ω = Interior(model) 
Γ = Boundary(model,tags="surface") 

# Define a mock function
fa(x) = x[1] + x[2]*2

# Interpolate the mock function onto the "interior"
reffe = ReferenceFE(lagrangian, Float64, 2)
fΩ = interpolate_everywhere(fa, FESpace(Ω,reffe))

# "evaluate" at an arbitrary point within the "interior" FESpace
@show fΩ( Point(0.52, 0.5) ) 

# Interpolate the mock function onto the "surface"
fΓ = interpolate_everywhere(fa, FESpace(Γ,reffe))

# "evaluate" at an arbitrary point on the "surface" FESpace
@show fΓ( Point(0.52, 0.5) ) 


# The above evalute call fails in Gridap v0.17.16, 
# with the following error
# ERROR: LoadError: DimensionMismatch: matrix is not square: dimensions are (1, 2)

# Corrections
# Modified src/Fields/AffineMaps.jl

end