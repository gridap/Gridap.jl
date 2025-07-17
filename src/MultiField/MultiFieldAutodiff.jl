
function CellData.SkeletonCellFieldPair(
  cf_plus::Union{MultiFieldFEFunction,MultiFieldCellField},
  cf_minus::Union{MultiFieldFEFunction,MultiFieldCellField}
)
  cfs = map(SkeletonCellFieldPair,cf_plus,cf_minus)
  MultiFieldCellField(cfs)
end

#################################
## Split and monolithic multifield autodiff
# - split: compute the gradient for each field separately
# - monolithic: compute the gradient for all fields together (original)
#
# For most problems, the split version is faster because the ForwardDiff
# chunk size is smaller. In addition, the split version allows for different
# triangulations.
#
# TODO: Currently, this is only implemented for the gradient and jacobian.
#  The Hessian is slightly proplematic because the off-diagonal blocks are
#  missed. This is because the basis isn't baked into f as it is in jacobian.

grad_ops = [
  (;op=:(FESpaces.gradient),split=:_mf_grad_split,mono=:(FESpaces._gradient)),
  (;op=:(FESpaces.jacobian),split=:_mf_jac_split,mono=:(FESpaces._jacobian )),
  # (;op=:(FESpaces.hessian ),split=:_mf_hes_split,mono=:(FESpaces._hessian  )),
]
for op in grad_ops
  @eval begin
    function $(op.op)(f::Function,uh::MultiFieldFEFunction;ad_type=:split)
      fuh = f(uh)
      if ad_type == :split
        $(op.split)(f,uh,fuh)
      elseif ad_type == :monolithic
        $(op.mono)(f,uh,fuh)
      else
        @notimplemented """Unknown ad_type = $ad_type
          Options:
          - :split      -- compute the gradient for each field separately
          - :monolithic -- compute the gradient for all fields together
          """
      end
    end

    function $(op.split)(f,uh,fuh)
      nfields = num_fields(uh)
      terms = map(Base.OneTo(nfields)) do k
        f_k = restrict_function(f,uh,k)
        $(op.mono)(f_k,uh[k],f_k(uh[k]))
      end
      return _combine_contributions($(op.op),terms,fuh)
    end
  end
end

# Helpers
function restrict_function(f,uh,k)
  uk->f((uh[1:k-1]...,uk,uh[k+1:end]...))
end

function concat_contribs_vec(size,contribs...)
  ArrayBlock([contribs...],size)
end

function concat_contribs_mat(touched,contribs...)
  mat = map(ij->contribs[ij[2]][ij[1]],CartesianIndices(touched))
  ArrayBlock(mat,touched)
end

GetIndex(k) = i->getindex(i,k)

function _combine_contributions(::typeof(gradient),terms::Vector{DomainContribution},fuh::DomainContribution)
  contribs = DomainContribution()
  nfields = length(terms)
  for trian in get_domains(fuh)
    trian_to_contrib = lazy_map(GetIndex(trian),terms)
    contrib_to_touched = fill([true for _ in 1:nfields],length(first(trian_to_contrib)));
    mf_cell_grad = lazy_map(concat_contribs_vec,contrib_to_touched,trian_to_contrib...);
    add_contribution!(contribs,trian,mf_cell_grad)
  end
  contribs
end

function _combine_contributions(::Union{typeof(jacobian),typeof(hessian)},terms::Vector{DomainContribution},fuh::DomainContribution)
  contribs = DomainContribution()
  nfields = length(terms)
  for trian in get_domains(fuh)
    trian_to_contrib = lazy_map(GetIndex(trian),terms)
    contrib_to_touched = fill(ones(Bool,nfields,nfields),length(first(trian_to_contrib)));
    mf_cell_grad = lazy_map(concat_contribs_mat,contrib_to_touched,trian_to_contrib...);
    add_contribution!(contribs,trian,mf_cell_grad)
  end
  contribs
end