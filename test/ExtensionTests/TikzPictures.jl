module TikzPicturesTests

using TikzPictures
using Gridap
using Gridap.Geometry: DiscreteModelMock

d = mktempdir()

model = DiscreteModelMock()
tp = TikzPicture(model)
# TikzPictures.save(PDF(joinpath(d,"tikz_model_2d")), tp)

model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
tp = TikzPicture(model)
# TikzPictures.save(PDF(joinpath(d,"tikz_model_3d")), tp)

poly = QUAD
tp = TikzPicture(poly)
# TikzPictures.save(PDF(joinpath(d,"tikz_poly_2d")), tp)

end