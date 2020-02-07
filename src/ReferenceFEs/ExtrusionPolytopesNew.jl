
struct NFace
  anchor::Vector{Int}
  extrusion::Vector{Int}
end

struct DFace
  extrusion::Vector{Int}
  nfaces::Vector{NFace}
  nf_nfs::Vector{Vector{Int64}}
  nf_dim::Vector{Vector{UnitRange{Int64}}}
end

function _polytopenfaces(anchor, extrusion)
  D = length(extrusion)
  dnf = _nfdim(extrusion)
  zerop = Point{D,Int}(zeros(Int64, D))
end

function _nfdim(extrusion)
  z = zero(eltype(extrusion))
  for e in extrusion
    if e > 0
      z +=1
    end
  end
  z
end


