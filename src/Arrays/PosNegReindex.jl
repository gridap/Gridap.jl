# This Map has non-trivial domain, thus we need the define testargs
"""
    PosNegReindex(values_pos,values_neg)
"""
struct PosNegReindex{A,B} <: Map
  values_pos::A
  values_neg::B
end

function testargs(k::PosNegReindex,i::Integer)
  @check length(k.values_pos) !=0 || length(k.values_neg) != 0 "This map has empty domain"
  @check eltype(k.values_pos) == eltype(k.values_neg) "This map is type-instable"
  length(k.values_pos) !=0 ? (one(i),) : (-one(i))
end

function return_value(k::PosNegReindex,i::Integer)
  if length(k.values_pos)==0 && length(k.values_neg)==0
    @check eltype(k.values_pos) == eltype(k.values_neg) "This map is type-instable"
    testitem(k.values_pos)
  else
    evaluate(k,testargs(k,i)...)
  end
end

function return_cache(k::PosNegReindex,i::Integer)
  c_p = array_cache(k.values_pos)
  c_n = array_cache(k.values_neg)
  c_p, c_n
end

@inline function evaluate!(cache,k::PosNegReindex,i::Integer)
  c_p, c_n = cache
  i>0 ? getindex!(c_p,k.values_pos,i) : getindex!(c_n,k.values_neg,-i)
end

@inline function evaluate(k::PosNegReindex,i::Integer)
  i>0 ? k.values_pos[i] : k.values_neg[-i]
end

function lazy_map(k::PosNegReindex,::Type{T},i_to_iposneg::AbstractArray) where T
  if all_pos(i_to_iposneg)
    lazy_map(Reindex(k.values_pos),get_array(i_to_iposneg))
  elseif all_neg(i_to_iposneg)
    lazy_map(Reindex(k.values_neg),lazy_map(iposneg->-iposneg,get_array(i_to_iposneg)))
  else
    s = size(i_to_iposneg)
    lazy_map(evaluate,T,Fill(k, s), i_to_iposneg)
  end
end

function lazy_map(::typeof(evaluate),::Type{T},a::LazyArray{<:Fill{<:PosNegReindex}}...) where T
  i_to_iposneg = a[1].args[1]
  if all(map( ai-> is_exhaustive(a[1].args[1]),a)) && all( map( ai-> i_to_iposneg==a[1].args[1],a) )
    bpos = map(ai->ai.maps.value.values_pos,a)
    bneg = map(ai->ai.maps.value.values_neg,a)
    cpos = lazy_map(evaluate,T,bpos...)
    cneg = lazy_map(evaluate,T,bneg...)
    lazy_map(PosNegReindex(cpos,cneg),T,i_to_iposneg)
  else
    LazyArray(T,a...)
  end
end

function lazy_map(::typeof(evaluate),::Type{T},a::LazyArray{<:Fill{<:PosNegReindex}},x::AbstractArray) where T
  apos = a.maps.value.values_pos
  aneg = a.maps.value.values_neg
  i_to_iposneg = a.args[1]
  ipos_to_i, ineg_to_i = pos_and_neg_indices(i_to_iposneg)
  xpos = lazy_map(Reindex(x),ipos_to_i)
  xneg = lazy_map(Reindex(x),ineg_to_i)
  cpos = lazy_map(evaluate,apos,xpos)
  cneg = lazy_map(evaluate,aneg,xneg)
  lazy_map(PosNegReindex(cpos,cneg),T,i_to_iposneg)
end

function lazy_map(::typeof(evaluate),::Type{T},b::Fill,a::LazyArray{<:Fill{<:PosNegReindex}}...) where T
  i_to_iposneg = a[1].args[1]
  if all(map( ai-> is_exhaustive(ai.args[1]),a)) && all( map( ai-> i_to_iposneg==ai.args[1],a) )
    k = b.value
    bpos = map(ai->ai.maps.value.values_pos,a)
    bneg = map(ai->ai.maps.value.values_neg,a)
    cpos = lazy_map(k,T,bpos...)
    cneg = lazy_map(k,T,bneg...)
    lazy_map(PosNegReindex(cpos,cneg),T,i_to_iposneg)
  else
    LazyArray(T,b,a...)
  end
end

# Helper functions to work with arrays representing a binary partition

function pos_and_neg_length(i_to_iposneg)
  Npos = maximum(i_to_iposneg)
  Nneg = -minimum(i_to_iposneg)
  Npos, Nneg
end

function is_exhaustive(i_to_iposneg)
  Npos, Nneg = pos_and_neg_length(i_to_iposneg)
  if length(i_to_iposneg) != Npos+Nneg
    return false
  end
  ipos_to_touched = fill(false,Npos)
  ineg_to_touched = fill(false,Nneg)
  for iposneg in i_to_iposneg
    iposneg>0 ? ipos_to_touched[iposneg] = true : ineg_to_touched[-iposneg] = true
  end
  all(ipos_to_touched) && all(ineg_to_touched)
end

function pos_and_neg_indices(i_to_iposneg)
  @check is_exhaustive(i_to_iposneg)
  Npos, Nneg = pos_and_neg_length(i_to_iposneg)
  ipos_to_i = zeros(Int,Npos)
  ineg_to_i = zeros(Int,Nneg)
  @inbounds for (i,iposneg) in enumerate(i_to_iposneg)
    iposneg>0 ? ipos_to_i[iposneg] = i : ineg_to_i[-iposneg] = i
  end
  ipos_to_i, ineg_to_i
end

function aligned_with_pos(i_to_iposneg,j_to_i,npos)
  j_to_iposneg = lazy_map(Reindex(i_to_iposneg),j_to_i)
  j_to_iposneg == 1:npos
end

function aligned_with_neg(i_to_iposneg,j_to_i,nneg)
  j_to_iposneg = lazy_map(Reindex(i_to_iposneg),j_to_i)
  j_to_iposneg == -(1:nneg)
end

function all_in_pos(i_to_iposneg,j_to_i)
  j_to_iposneg = lazy_map(Reindex(i_to_iposneg),j_to_i)
  all_pos(j_to_iposneg)
end

function all_in_neg(i_to_iposneg,j_to_i)
  j_to_iposneg = lazy_map(Reindex(i_to_iposneg),j_to_i)
  all_neg(j_to_iposneg)
end

all_pos(i_to_iposneg) = all( lazy_map(iposneg->iposneg>0, i_to_iposneg))
all_neg(i_to_iposneg) = all( lazy_map(iposneg->iposneg<0, i_to_iposneg))

# This is important to do optimizations associated with ExtendedFESpace
"""
struct representing a binary partition of a range of indices

Using this allows one to do a number of important optimizations when working with `PosNegReindex`
"""
struct PosNegPartition{T,V,Vp,Vn} <: AbstractVector{T}
  i_to_iposneg::V
  ipos_to_i::Vp
  ineg_to_i::Vn
  function PosNegPartition(
    i_to_iposneg::AbstractVector,
    ipos_to_i::AbstractVector,
    ineg_to_i::AbstractVector)

    @check eltype(ipos_to_i) == eltype(ineg_to_i)
    @check all((ipos_to_i, ineg_to_i) .== pos_and_neg_indices(i_to_iposneg))
    T = eltype(i_to_iposneg)
    @check T <: Integer
    V = typeof(i_to_iposneg)
    Vp = typeof(ipos_to_i)
    Vn = typeof(ineg_to_i)
    new{T,V,Vp,Vn}(i_to_iposneg,ipos_to_i,ineg_to_i)
  end
end

function PosNegPartition(ipos_to_i::AbstractArray,Ni::Integer)
  i_to_iposneg = zeros(typeof(Ni),Ni)
  i_to_iposneg[ipos_to_i] .= 1:length(ipos_to_i)
  Nneg = count(k->k==0,i_to_iposneg)
  ineg_to_i = similar(ipos_to_i,eltype(ipos_to_i),(Nneg,))
  ineg = 1
  for (i,iposneg) in enumerate(i_to_iposneg)
    if iposneg==0
      i_to_iposneg[i] = -ineg
      ineg_to_i[ineg] = i
      ineg += 1
    end
  end
  PosNegPartition(i_to_iposneg,ipos_to_i,ineg_to_i)
end

Base.IndexStyle(::Type{<:PosNegPartition}) = IndexLinear()
@propagate_inbounds Base.getindex(a::PosNegPartition,i::Integer) = a.i_to_iposneg[i]
Base.size(a::PosNegPartition) = size(a.i_to_iposneg)
get_array(a::PosNegPartition) = a.i_to_iposneg

pos_and_neg_length(a::PosNegPartition) = (length(a.ipos_to_i),length(a.ineg_to_i))
is_exhaustive(a::PosNegPartition) = true
pos_and_neg_indices(a::PosNegPartition) = (a.ipos_to_i,a.ineg_to_i)
function aligned_with_pos(a::PosNegPartition,j_to_i,npos)
  @check length(a.ipos_to_i) == npos
  a.ipos_to_i === j_to_i || a.ipos_to_i == j_to_i
end
function aligned_with_neg(a::PosNegPartition,j_to_i,nneg)
  @check length(a.ineg_to_i) == nneg
  a.ineg_to_i === j_to_i || a.ineg_to_i == j_to_i
end
all_pos(a::PosNegPartition) = length(a.ineg_to_i)==0
all_neg(a::PosNegPartition) = length(a.ipos_to_i)==0

# @propagate_inbounds function Base.setindex!(a::LazyArray{<:Fill{<:PosNegReindex}},v,j::Integer)
#   k = a.map.value
#   i_to_v = k.values
#   j_to_i, = a.f
#   i = j_to_i[j]
#   i_to_v[i]=v
# end
