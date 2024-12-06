
struct DorflerMarking
  θ :: Float64
  ν :: Float64
  strategy :: Symbol
  function DorflerMarking(
    θ::Float64;
    ν::Float64 = 0.5,
    strategy::Symbol = :sort
  )
    @assert 0 < θ < 1
    @assert strategy ∈ (:sort,:binsort,:quickmark) "Strategy not recognized. Available values are (:sort)"
    new(θ,ν,strategy)
  end
end

mark(m::DorflerMarking, η::Vector{<:Real}) = mark(Val(m.strategy), m, η)

function mark(::Val{:sort}, m::DorflerMarking, η::Vector{<:Real})
  target = m.θ * sum(η)
  perm = sortperm(η, rev=true, alg=QuickSort)
  s = zero(eltype(η))
  k = 0
  while s < target
    k += 1
    s += η[perm[k]]
  end
  return perm[1:k]
end

function mark(::Val{:binsort}, m::DorflerMarking, η::Vector{<:Real})
  target = m.θ * sum(η)
  M = maximum(η)
  N = length(η)
  
  # Find minimal K such that
  #   νᴷ⁺¹ M ≤ (1 - θ) / (θ * N) * target
  K = 0
  a = M*m.ν
  b = (1 - m.θ) / (m.θ * N * M) * target
  while a > b
    K += 1
    a *= m.ν
  end

  # Sort into bins Bk such that ηi ∈ Bk if
  #    νᴷ⁺¹ M ≤ ηi < νᴷ M
  bins = zeros(Int,N)
  sums = zeros(Float64,K+1)
  lbs = zeros(Float64,K+1)
  lbs[1] = M * m.ν
  for i in 2:K
    lbs[i] = lbs[i-1] * m.ν
  end
  lbs[end] = 0.0

  for (i,ηi) in enumerate(η)
    k = 1
    while ηi < lbs[k]
      k += 1
    end
    bins[i] = k
    sums[k] += ηi
  end

  # Find minimal set of bins that gets over target
  k = 0
  s = 0.0
  while s < target
    k += 1
    s += sums[k]
  end
  
  return findall(i -> i <= k, bins)
end

function mark(::Val{:quickmark}, m::DorflerMarking, η::Vector{<:Real})
  function quickmark!(η, perm, l, u, target) :: Int
    m = (u - l) ÷ 2
    sort!(view(perm,l:u), by=i->η[i], rev=true, alg=PartialQuickSort(m))

    p = l + m
    t = η[perm[p]]
    σ = sum(η[perm[l:p-1]])
    
    (σ >= target) && return quickmark!(η, perm, l, p, target)
    (σ + t >= target) && return p
    return quickmark!(η, perm, p + 1, u, target - σ - t)
  end

  N = length(η)
  l, u = 1, N
  perm = collect(1:N)
  target = m.θ * sum(η)
  m = quickmark!(η, perm, l, u, target)
  return perm[1:m]
end
