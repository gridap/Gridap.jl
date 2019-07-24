module FaceLabelsModule

using Gridap
using Gridap.CellValuesGallery

export FaceLabels
export labels_on_dim
export labels_on_tag
export ntags
export tag_from_name
export name_from_tag

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

end #module
