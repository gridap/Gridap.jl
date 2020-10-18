# This Map has non-trivial domain, thus we need the define testargs
"""
    Reindex(values) -> Map
"""
struct Reindex{A} <: Map
  values::A
end

function testargs(k::Reindex,i)
  @check length(k.values) !=0 "This map has empty domain"
  (one(i),)
end
function testargs(k::Reindex,i::Integer...)
  @check length(k.values) !=0 "This map has empty domain"
  map(one,i)
end
function return_value(k::Reindex,i...)
  length(k.values)!=0 ? evaluate(k,testargs(k,i...)...) : testitem(k.values)
end
return_cache(k::Reindex,i...) = array_cache(k.values)
@inline evaluate!(cache,k::Reindex,i...) = getindex!(cache,k.values,i...)

#"""
#    reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
#"""
#function reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
#  lazy_map(Reindex(i_to_v),j_to_i)
#end

function lazy_map(k::Reindex{<:Fill}, j_to_i::AbstractArray)
  v = k.values.value
  Fill(v,size(j_to_i)...)
end

function lazy_map(k::Reindex{<:CompressedArray}, j_to_i::AbstractArray)
  i_to_v = k.values
  values = i_to_v.values
  ptrs = lazy_map(Reindex(i_to_v.ptrs),j_to_i)
  CompressedArray(values,ptrs)
end

# This optimization is important for surface-coupled problems
function lazy_map(k::Reindex{<:LazyArray{<:Fill{<:PosNegReindex}}}, j_to_i::AbstractArray)
  i_to_iposneg = k.values.f[1]
  ipos_to_value = k.values.g.value.values_pos
  ineg_to_value = k.values.g.value.values_neg
  if aligned_with_pos(i_to_iposneg,j_to_i,length(ipos_to_value))
    ipos_to_value
  elseif aligned_with_neg(i_to_iposneg,j_to_i,length(ineg_to_value))
    ineg_to_value
  elseif all_pos(i_to_iposneg,j_to_i)
    j_to_ipos = lazy_map(Reindex(get_array(i_to_iposneg)),j_to_i)
    j_to_value = lazy_map(Reindex(ipos_to_value),j_to_ipos)
  elseif all_neg(i_to_iposneg,j_to_i)
    j_to_ineg = lazy_map(Reindex(get_array(i_to_iposneg)),j_to_i)
    j_to_value = lazy_map(Reindex(ineg_to_value),lazy_map(ineg->-ineg,j_to_ineg))
  else
    j_to_iposneg = lazy_map(Reindex(get_array(i_to_iposneg)),j_to_i)
    j_to_value = lazy_map(PosNegReindex(ipos_to_value,ineg_to_value),j_to_iposneg)
  end
end

function lazy_map(k::Reindex{<:LazyArray{<:PosNegReindex}}, j_to_i::IdentityVector)
  @check length(k.values) == length(indices)
  k.values
end

function lazy_map(k::Reindex{<:AbstractArray}, indices::IdentityVector)
  @check length(k.values) == length(indices)
  k.values
end

function lazy_map(k::Reindex{<:Fill},b::IdentityVector)
  @check length(k.values) == length(indices)
  k.values
end

function lazy_map(k::Reindex{<:CompressedArray},b::IdentityVector)
  @check length(k.values) == length(indices)
  k.values
end

@propagate_inbounds function Base.setindex!(a::LazyArray{<:Fill{<:Reindex}},v,j::Integer)
  k = a.g.value
  i_to_v = k.values
  j_to_i, = a.f
  i = j_to_i[j]
  i_to_v[i]=v
end

