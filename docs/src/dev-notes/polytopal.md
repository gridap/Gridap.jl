
# Development notes for polytopal methods

## Overlapping boundary triangulations

- At the moment, boundary triangulations cannot have repeated faces. This is enforced by the constructor, which takes as input the `bgface_to_lcell` array which ties a single cell to each model face.

- The reality, however, is that `bgface_to_lcell` is only used within the constructor. So we can replace it (with some minor changes to the constructor) by a new array `face_to_lcell`, which allows repeated faces with different neighboring cells attached to them.

- Allowing for overlapping boundary triangulations would be useful in many situations other that polytopal methods. This is something we have talked bout with Alberto, and I actually already use it for patch-based solvers.

- So I think it would be nice to just re-work the `FaceToCellGlue` to allow for repeated faces. This would be a breaking change, but I think it is worth it.

- For now, I will simply work around it until we reach an understanding.
