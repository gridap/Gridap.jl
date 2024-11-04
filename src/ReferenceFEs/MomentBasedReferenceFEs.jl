
"""
For each moment, I need: 
 - its dimension
 - the basis we are integrating against
 - a function (φ,ϕ) -> the functional to integrate
 - a face filter/mask
"""
function MomentBasedReferenceFE(
  p::Polytope{D},
  face_moments::AbstractVector{<:AbstractVector{<:Function}},
  face_basis::AbstractVector{<:AbstractVector{<:AbstractVector{<:Field}}};
  face_mask::AbstractVector{<:AbstractVector{Bool}} = [fill(true,num_faces(p, d)) for d in 0:D],
) where D
  @assert length(face_moments) == length(face_mask) == length(face_basis) == D+1
  @assert all(length(face_moments[d+1]) == length(face_basis[d+1]) for d in 0:D)
  @assert all(length(face_mask[d+1]) == num_faces(p, d) for d in 0:D)

end



