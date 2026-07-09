
using Gridap
using TikzPictures

using Gridap.ReferenceFEs, Gridap.Adaptivity

function export_tikz_poly(name, p, d, label_offset)
  options = "background rectangle/.style={draw=white,fill=white,rounded corners=1ex}, show background rectangle, color=black, scale=2"
  preamble = "\\usetikzlibrary {backgrounds}"
  tp = TikzPicture(p;options,preamble,label_dim=d,label_offset)
  TikzPictures.save(SVG(joinpath("src/assets/polytopes/",name)), tp)
end

export_tikz_poly("SEGMENT_0", Polygon([Point(0.0,0.0),Point(1.0,0.0)]), 0, Point(0.0,0.1))

for d in 0:1
  offset = iszero(d) ? 0.1 : 0.0
  export_tikz_poly("QUAD_$d", QUAD, d, offset)
  export_tikz_poly("TRI_$d" , TRI, d, offset)
end

for d in 0:2
  offset = iszero(d) ? 0.1 : 0.0
  export_tikz_poly("HEX_$d", HEX, d, offset)
  export_tikz_poly("TET_$d" , TET, d, offset)
  export_tikz_poly("PYRAMID_$d", PYRAMID, d, offset)
  export_tikz_poly("WEDGE_$d", WEDGE, d, offset)
end
