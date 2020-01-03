
"""
"""
struct SkeletonPair{L,R} <: GridapType
  left::L
  right::R
end

struct SkeletonField <: Field end

struct SkeletonReduction{A,B,C} <: Kernel
  funleft::A
  funright::B
  fun::C
end

function apply_kernel!(cache,k::SkeletonReduction,left,right)
  SkeletonField()
end

"""
"""
function jump(a::SkeletonPair)
  k = SkeletonReduction(identity,-,-)
  apply(k,a.left,a.right)
end

"""
"""
function mean(a::SkeletonPair)
  k = SkeletonReduction(_mean,_mean,_mean)
  apply(k,a.left,a.right)
end

_mean(x) = 0.5*x
_mean(x,y) = 0.5*x + 0.5*y

function kernel_evaluate(k::SkeletonReduction,x,left,right)
  lx = evaluate_field_array(left,x)
  rx = evaluate_field_array(right,x)
  _apply_skeleton_reduction(k,lx,rx,testitem(lx))
end

function _apply_skeleton_reduction(k,lx,rx,::AbstractVector)
  apply(bcast(k.fun),lx,rx)
end

function _apply_skeleton_reduction(k,lx,rx,::AbstractMatrix)
  left = apply(bcast(k.funleft),lx)
  right = apply(bcast(k.funright),rx)
  SkeletonPair(left,right)
end

