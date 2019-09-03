module FaceLabelsModule

using Gridap
using Gridap.CellValuesGallery

export FaceLabels
export labels_on_dim
export labels_on_tag
export ntags
export tag_from_name
export name_from_tag
export first_tag_on_face

"""
Classification of nfaces into geometrical and physical labels
"""
struct FaceLabels
  dim_to_nface_to_label::Vector{IndexCellValue}
  tag_to_labels::Vector{Vector{Int}}
  tag_to_name::Vector{String}
end

function FaceLabels(
  dim_to_nface_to_label::Vector{Vector{Int}},
  physlabel_to_labels::Vector{Vector{Int}},
  tag_to_name::Vector{String})

  cv = [ CellValueFromArray(v) for v in dim_to_nface_to_label ]
  FaceLabels(cv, physlabel_to_labels, tag_to_name)
end

"""
Returns an AbstractVector{Int} that represent the label for
each nface of dimension dim
"""
labels_on_dim(l::FaceLabels,dim::Integer) = l.dim_to_nface_to_label[dim+1]

"""
Returns a Vector{Int} with the labels associated with a given phystag
"""
labels_on_tag(l::FaceLabels,tag::Integer) = l.tag_to_labels[tag]

ntags(l::FaceLabels) = length(l.tag_to_labels)

function tag_from_name(l::FaceLabels,name::String)
  for tag in 1:ntags(l)
    if l.tag_to_name[tag] == name
      return tag
    end
  end
  0
end

name_from_tag(l::FaceLabels,tag::Integer) = l.tag_to_name[tag]

function first_tag_on_face(labels::FaceLabels,dim::Integer)
  face_to_label = labels_on_dim(labels,dim)
  nfaces = length(face_to_label)
  face_to_tag = zeros(Int,nfaces)
  _first_tag_on_face!(face_to_tag,face_to_label,labels.tag_to_labels)
  CellValueFromArray(face_to_tag)
end

function _first_tag_on_face!(face_to_tag,face_to_label,tag_to_labels)
  nfaces = length(face_to_tag)
  for face in 1:nfaces
    label = face_to_label[face]
    for (tag,labels) in enumerate(tag_to_labels)
      if label in labels
        face_to_tag[face] = tag
      end
    end
  end
end

end #module
