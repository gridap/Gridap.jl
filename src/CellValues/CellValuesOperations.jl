module CellValuesOperations

using Gridap

using Gridap.Kernels: CellSumKernel
using Gridap.Kernels: CellNewAxisKernel

export cellsum
export cellmean
export cellnewaxis
import Base: +,-,*,/,\, ==, ≈
import LinearAlgebra: inv, det
import TensorValues: inner, outer, meas

for op in (:+,:-,:inv,:det,:meas)
  @eval begin

    function ($op)(m::CellNumber)
      apply($op,m)
    end

    function ($op)(m::CellArray)
      apply($op,m,broadcast=true)
    end

    function ($op)(m::CellMap)
      apply($op,m,broadcast=true)
    end

    function ($op)(m::Map)
      apply($op,m,broadcast=true)
    end

  end
end

for op in (:+,:-,:*,:/,:\,:inner,:outer)
  @eval begin

    function ($op)(a::CellNumber,b::CellNumber)
      apply($op,a,b)
    end

    function ($op)(a::CellArray,b::CellArray)
      apply($op,a,b,broadcast=true)
    end

    function ($op)(a::CellArray,b::CellNumber)
      apply($op,a,b,broadcast=true)
    end

    function ($op)(a::CellNumber,b::CellArray)
      apply($op,a,b,broadcast=true)
    end

    function ($op)(a::CellMap,b::CellMap)
      apply($op,a,b,broadcast=true)
    end

    function ($op)(a::Map,b::Map)
      apply($op,a,b,broadcast=true)
    end

  end
end

for op in (:(==),:≈)
  @eval begin

    function ($op)(a::CellNumber,b::CellNumber)
      _eq_kernel($op,a,b)
    end

    function ($op)(a::CellArray,b::CellArray)
      _eq_kernel($op,a,b)
    end

  end
end

@inline function _eq_kernel(op,a,b)
  length(a) != length(b) && return false
  for (ai,bi) in zip(a,b)
    !(op(ai,bi)) && return false
  end
  return true
end

function cellsum(self::CellArray{T,N};dim::Int) where {T,N}
  k = CellSumKernel{dim}()
  apply(k,self)
end

function cellsum(self::CellArray{T,1};dim::Int) where T
  apply(sum,self)
end

_mean(a) = sum(a)/length(a)

function cellmean(self::CellArray)
  apply(_mean,self)
end

function cellnewaxis(ca::CellArray;dim::Int)
  k = CellNewAxisKernel{dim}()
  apply(k,ca)
end

function cellnewaxis(ca::CellMap;dim::Int)
  k = CellNewAxisKernel{dim}()
  apply(k,ca)
end

end # module
