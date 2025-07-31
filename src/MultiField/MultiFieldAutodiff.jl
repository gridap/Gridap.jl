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
# chunk size is smaller. In addition, the split version allows for fields to
# ve defined on different triangulations.
#
# TODO: Currently, this is only implemented for the gradient and jacobian.
#  The Hessian is proplematic because the off-diagonal blocks are missed.
for (op,_op) in ((:gradient,:_gradient),(:jacobian,:_jacobian))
  @eval begin
    function FESpaces.$(op)(f::Function,uh::MultiFieldFEFunction;ad_type=:split)
      fuh = f(uh)
      if ad_type == :split
        multifield_autodiff_split($op,f,uh,fuh)
      elseif ad_type == :monolithic
        FESpaces.$(_op)(f,uh,fuh)
      else
        @notimplemented UNKNOWN_AD_TYPE(ad_type)
      end
    end
  end
end

# Defaults to monolithic for the Hessian
function FESpaces.hessian(f::Function,uh::MultiFieldFEFunction;ad_type=:monolithic)
  fuh = f(uh)
  if ad_type == :split
    @notimplemented "Hessian AD with ad_type = :split is not currently implemented. Use ad_type = :monolithic instead."
  elseif ad_type == :monolithic
    FESpaces._hessian(f,uh,fuh)
  else
    @notimplemented UNKNOWN_AD_TYPE(ad_type)
  end
end

for (op,_op) in ((:gradient,:_gradient),(:jacobian,:_jacobian))
  @eval begin

    function multifield_autodiff_split(::typeof($op),f,uh,fuh)
      nfields = num_fields(uh)
      terms = map(Base.OneTo(nfields)) do k
        # Although technically wrong, we can reuse fuh for each field since
        # we only use it to extract the triangulations
        f_k = uk -> f((uh[1:k-1]...,uk,uh[k+1:end]...))
        FESpaces.$(_op)(f_k,uh[k],fuh)
      end
      return _combine_contributions($op,terms,fuh)
    end

  end
end

function _combine_contributions(::typeof(gradient),terms,fuh::DomainContribution)
  contribs = DomainContribution()
  nfields = length(terms)
  block_map = BlockMap(nfields,collect(Base.OneTo(nfields)))
  for trian in get_domains(fuh)
    sf_contributions = map(Base.Fix2(get_contribution,trian),terms)
    mf_contribution = lazy_map(block_map,sf_contributions...)
    add_contribution!(contribs,trian,mf_contribution)
  end
  contribs
end

function _combine_contributions(::typeof(jacobian),terms,fuh::DomainContribution)
  contribs = DomainContribution()
  nfields = length(terms)
  I = [[(CartesianIndex(j),CartesianIndex(j,i)) for j in 1:nfields] for i in 1:nfields]
  block_map = Arrays.MergeBlockMap((nfields,nfields),I)
  for trian in get_domains(fuh)
    sf_contributions = map(Base.Fix2(get_contribution,trian),terms)
    mf_contribution = lazy_map(block_map,sf_contributions...)
    add_contribution!(contribs,trian,mf_contribution)
  end
  contribs
end

UNKNOWN_AD_TYPE(ad_type) = """Unknown ad_type = $ad_type
  Options:
  - :split      -- compute the gradient for each field separately, then merge
  - :monolithic -- compute the gradient for all fields together
"""