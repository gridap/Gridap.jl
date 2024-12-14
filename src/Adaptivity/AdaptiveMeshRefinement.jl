
# Marking strategies

"""
    struct DorflerMarking
      θ :: Float64
      ν :: Float64
      strategy :: Symbol
    end

    DorflerMarking(θ::Float64; ν::Float64 = 0.5, strategy::Symbol = :quickmark)

Implements the Dorfler marking strategy. Given a vector `η` of real positive numbers, 
the marking strategy find a subset of indices `I` such that 

  sum(η[I]) > θ * sum(η)

where `0 < θ < 1` is a threshold parameter. 

For more details, see the following reference: 

"Dörfler marking with minimal cardinality is a linear complexity problem", Pfeiler et al. (2020)

The marking algorithm is controlled by the `strategy` parameter, which can take 
the following values: 

- `:sort`: Optimal cardinality, O(N log N) complexity. See Algorithm 2 in the reference.
- `:binsort`: Quasi-optimal cardinality, O(N) complexity. See Algorithm 7 in the reference.
- `:quickmark`: Optimal cardinality, O(N) complexity.  See Algorithm 10 in the reference.

# Arguments

- `θ::Float64`: The threshold parameter. Between 0 and 1.
- `ν::Float64`: Extra parameter for `:binsort`. Default is 0.5.
- `strategy::Symbol`: The marking strategy. Default is `:quickmark`.

# Usage

```julia
η = abs.(randn(1000))
m = DorflerMarking(0.5)
I = mark(m,η)
```

"""
struct DorflerMarking
  θ :: Float64
  ν :: Float64
  strategy :: Symbol
  function DorflerMarking(
    θ::Float64;
    ν::Float64 = 0.5,
    strategy::Symbol = :quickmark
  )
    @assert 0 < θ < 1
    @assert strategy ∈ (:sort,:binsort,:quickmark) "Strategy not recognized. Available values are (:sort)"
    new(θ,ν,strategy)
  end
end

"""
    mark(m::DorflerMarking, η::Vector{<:Real}) -> Vector{Int}

Given a vector `η` of real positive numbers, returns a subset of indices `I` such that 
satisfying the Dorfler marking condition.
"""
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

# Estimators

"""
    estimate(f::Function, uh::Function) -> Vector{Float64}

Given a functional `f` and a function `uh`, such that `f(uh)` produces a 
scalar-valued `DomainContribution`, collects the estimator values for 
each cell in the background model.
"""
function estimate(f::Function, uh)
  collect_estimator(f(uh))
end

function collect_estimator(c::DomainContribution)
  trians = get_domains(c)
  bgmodel = get_background_model(first(trians))
  msg = "Estimator not implemented for mixed background models"
  @notimplementedif !all([bgmodel == get_background_model(trian) for trian in trians]) msg

  Dc = num_cell_dims(bgmodel)
  η = zeros(Float64,num_cells(bgmodel))
  for trian in trians
    glue = get_glue(trian,Val(Dc))
    collect_estimator!(η,glue,get_contribution(c,trian))
  end

  return η
end

function collect_estimator!(η, glue::FaceToFaceGlue, c)
  cache = array_cache(c)
  for (face,bgcell) in enumerate(glue.tface_to_mface)
    η[bgcell] += getindex!(cache,c,face)
  end
end

function collect_estimator!(η, glue::SkeletonPair, c)
  collect_estimator!(η,glue.plus,c)
  collect_estimator!(η,glue.minus,c)
end
