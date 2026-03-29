
# ---------------------------------------------------------------------------
# ConstraintLine
# ---------------------------------------------------------------------------

# One constrained (slave) DoF and its linear relation to master DoFs.
struct ConstraintLine{T<:Number}
  dof::Int
  masters::Vector{Int}
  coeffs::Vector{T}
  offset::T
end

function Base.show(io::IO, l::ConstraintLine)
  terms = join(["$(c)*x[$(d)]" for (d, c) in zip(l.masters, l.coeffs)], " + ")
  rhs = isempty(terms) ? "$(l.offset)" :
      iszero(l.offset) ? terms : "$terms + $(l.offset)"
  print(io, "x[$(l.dof)] = $rhs")
end

# ---------------------------------------------------------------------------
# ConstraintHandler
# ---------------------------------------------------------------------------

mutable struct ConstraintHandler{T<:Number}
  n_dofs::Int
  constraints::Vector{ConstraintLine{T}}
  dof_to_constraint::Vector{Int} # 0 if unconstrained
  closed::Bool
end

is_constrained(ch::ConstraintHandler, dof::Int) = !iszero(ch.dof_to_constraint[dof])

num_constrained_dofs(ch::ConstraintHandler) = length(ch.constraints)
num_free_dofs(ch::ConstraintHandler)        = ch.n_dofs - num_constrained_dofs(ch)

function Base.show(io::IO, ch::ConstraintHandler{T}) where T
  status = ch.closed ? "closed" : "open"
  print(io, "ConstraintHandler{$T}($(ch.n_dofs) dofs, " *
        "$(num_constrained_dofs(ch)) constrained, $status)")
end

"""
  ConstraintHandler(n_dofs, T=Float64)

Create an empty constraint handler for a space with `n_dofs` degrees of freedom.
`T` is the coefficient type.
"""
ConstraintHandler(n_dofs::Int, ::Type{T}=Float64) where T =
  ConstraintHandler{T}(n_dofs, ConstraintLine{T}[], zeros(Int, n_dofs), false)

"""
    ConstraintHandler(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V::FESpace, T=Float64)

Reconstruct an open `ConstraintHandler` from a Gridap Constructor-2 constraint
table triple defined on the `SingleFieldFESpace` `V`. This is the inverse of
`constraint_tables`.

The returned handler is open; call `close!` to chain-resolve and finalise it.
"""
function ConstraintHandler(
  sDOF_to_dof, sDOF_to_dofs::Table, sDOF_to_coeffs::Table,
  V::SingleFieldFESpace, ::Type{T}=Float64
) where T
  n_fdofs = num_free_dofs(V)
  n_dofs  = n_fdofs + num_dirichlet_dofs(V)
  ch      = ConstraintHandler(n_dofs, T)
  for (i, dof) in enumerate(sDOF_to_dof)
    slave_DOF = _dof_to_DOF(dof, n_fdofs)
    rd = dataview(sDOF_to_dofs, i)
    rc = dataview(sDOF_to_coeffs, i)
    masters = [_dof_to_DOF(d, n_fdofs) for d in rd]
    coeffs  = collect(T, rc)
    add_constraint!(ch, slave_DOF, masters, coeffs)
  end
  return ch
end

# ---------------------------------------------------------------------------
# Constraint addition
# ---------------------------------------------------------------------------

"""
    add_constraint!(ch, line::ConstraintLine)

Register a pre-built `ConstraintLine`. The other `add_constraint!` methods
delegate to this one.
"""
function add_constraint!(ch::ConstraintHandler{T}, line::ConstraintLine{T}) where T
  @check !ch.closed "Cannot add constraints after close!()"
  @check 1 <= line.dof <= ch.n_dofs "DoF $(line.dof) out of range [1, $(ch.n_dofs)]"
  @check !is_constrained(ch, line.dof) "DoF $(line.dof) is already constrained"
  @check length(line.masters) == length(line.coeffs) "masters and coeffs must have the same length"
  for master in line.masters
    @check 1 <= master <= ch.n_dofs "Master DoF $master out of range [1, $(ch.n_dofs)]"
    @check master != line.dof "Self-referential constraint on DoF $(line.dof)"
  end
  push!(ch.constraints, line)
  ch.dof_to_constraint[line.dof] = length(ch.constraints)
  return ch
end

"""
    add_constraint!(ch, dof, masters, coeffs, offset=zero(T))

Register the linear constraint  x[dof] = Σ coeffs[i] * x[masters[i]] + offset.

Both `masters` and `coeffs` are copied internally, so the caller can safely reuse 
or modify them after the call.
"""
function add_constraint!(
  ch::ConstraintHandler{T}, dof::Int,
  masters::AbstractVector{<:Integer}, 
  coeffs::AbstractVector{<:Number}, 
  offset::Number=zero(T)
) where T
  # collect ensures ConstraintLine owns its data regardless of what the caller passes in
  add_constraint!(ch, ConstraintLine{T}(dof, collect(Int, masters), collect(T, coeffs), T(offset)))
end

"""
    add_constraint!(ch, dof, offset)

Register a Dirichlet constraint  x[dof] = offset  (no master DoFs).
"""
function add_constraint!(ch::ConstraintHandler{T}, dof::Int, offset::Number) where T
  add_constraint!(ch, ConstraintLine{T}(dof, Int[], T[], T(offset)))
end

"""
    set_offset!(ch, dof, value)

Update the offset of an already-registered constraint on `dof`.
"""
function set_offset!(ch::ConstraintHandler{T}, dof::Int, value::Number) where T
  @check !ch.closed "Cannot modify constraints after close!()"
  @check is_constrained(ch, dof) "DoF $dof has no constraint to update"
  idx = ch.dof_to_constraint[dof]
  line = ch.constraints[idx]
  ch.constraints[idx] = ConstraintLine{T}(line.dof, line.masters, line.coeffs, T(value))
  return ch
end

"""
    merge!(a, b; on_conflict=nothing)

Absorb all constraints from `b` into `a`. Both handlers must be open.
Conflict resolution is optional and controlled by the user-defined function `on_conflict`,
which is called when the same DoF is constrained in both handlers.
Its signature should be

  `on_conflict(existing::ConstraintLine, incoming::ConstraintLine) → ConstraintLine`

Example resolvers:
  - keep existing:  `(e, _) -> e`
  - keep incoming:  `(_, i) -> i`
  - average:        `(e, i) -> ConstraintLine(e.dof, vcat(e.masters, i.masters), vcat(e.coeffs, i.coeffs) .* 0.5, 0.0)`

Chain resolution and cycle detection are deferred to `close!(a)`.
"""
function Base.merge!(a::ConstraintHandler{T}, b::ConstraintHandler{T}; on_conflict=nothing) where T
  @check !a.closed && !b.closed "Both handlers must be open to merge"
  @check a.n_dofs == b.n_dofs "Handlers have incompatible sizes: $(a.n_dofs) vs $(b.n_dofs)"
  for line in b.constraints
    dof = line.dof
    cid = a.dof_to_constraint[dof]
    if !is_constrained(a, dof)
      push!(a.constraints, line)
      a.dof_to_constraint[dof] = length(a.constraints)
    elseif isnothing(on_conflict)
      existing = a.constraints[cid]
      error("Constraint conflict on DoF $(dof):\n  existing: $existing\n  incoming: $line")
    else
      existing = a.constraints[cid]
      a.constraints[cid] = on_conflict(existing, line)
    end
  end
  return a
end

"""
    close!(ch; tol=1e-13)

Finalize the constraint handler:
1. Sort constraints by DoF index; rebuild dof_to_constraint.
2. Detect dependency cycles via topological sort (Kahn's algorithm).
3. Resolve constraint chains by iterative substitution.
4. Merge duplicate master entries, strip near-zero coefficients.
"""
function close!(ch::ConstraintHandler{T}; tol::Float64=1e-13) where T
  ch.closed && return ch

  # Step 1: sort, rebuild dof_to_constraint.
  sort!(ch.constraints, by=l -> l.dof)
  fill!(ch.dof_to_constraint, 0)
  for (i, line) in enumerate(ch.constraints)
    ch.dof_to_constraint[line.dof] = i
  end

  # Step 2: cycle detection.
  _detect_cycles(ch)

  # Step 3: chain resolution (multi-pass substitution).
  # After cycle detection we know the dependency graph is a DAG, so this
  # terminates in at most n_slaves passes.
  masters, coeffs = Int[], T[]
  isconst = Base.Fix1(is_constrained, ch)
  for _ in 1:length(ch.constraints)
    changed = false
    for i in eachindex(ch.constraints)
      line = ch.constraints[i]
      any(isconst, line.masters) || continue
      empty!(masters); empty!(coeffs)
      new_offset = line.offset
      for (master, coeff) in zip(line.masters, line.coeffs)
        if isconst(master)
          sub = ch.constraints[ch.dof_to_constraint[master]]
          append!(masters, sub.masters)
          append!(coeffs, coeff .* sub.coeffs)
          new_offset += coeff * sub.offset
        else
          push!(masters, master); push!(coeffs, coeff)
        end
      end
      resize!(line.masters, length(masters)); copyto!(line.masters, masters)
      resize!(line.coeffs,  length(coeffs)); copyto!(line.coeffs,  coeffs)
      # ConstraintLine is immutable so we must rebuild only when offset changes.
      new_offset != line.offset &&
        (ch.constraints[i] = ConstraintLine{T}(line.dof, line.masters, line.coeffs, new_offset))
      changed = true
    end
    !changed && break
  end

  # Step 4: merge duplicate masters, strip near-zeros.
  # Masters and coeffs are modified in-place — no new ConstraintLine needed.
  acc       = Dict{Int,T}()
  buf_pairs = Tuple{Int,T}[]
  for line in ch.constraints
    empty!(acc)
    for (master, coeff) in zip(line.masters, line.coeffs)
      acc[master] = get(acc, master, zero(T)) + coeff
    end
    empty!(buf_pairs)
    for (k, v) in acc
      abs(v) > tol && push!(buf_pairs, (k, v))
    end
    sort!(buf_pairs, by=first)
    resize!(line.masters, length(buf_pairs))
    resize!(line.coeffs,  length(buf_pairs))
    for (j, (k, v)) in enumerate(buf_pairs)
      line.masters[j] = k; line.coeffs[j] = v
    end
  end

  ch.closed = true
  return ch
end

# Cycle detection via Kahn's topological sort on the slave→slave subgraph.
# Edge direction: t → s means "t must be resolved before s" (s depends on t).
function _detect_cycles(ch::ConstraintHandler)
  in_degree = zeros(Int, ch.n_dofs)
  adj       = Dict{Int,Vector{Int}}()  # only populated for slave→slave edges

  for line in ch.constraints
    for master in line.masters
      if is_constrained(ch, master)
        push!(get!(adj, master, Int[]), line.dof)  # master must precede line.dof
        in_degree[line.dof] += 1
      end
    end
  end

  queue     = [line.dof for line in ch.constraints if iszero(in_degree[line.dof])]
  processed = 0
  while !isempty(queue)
    s = pop!(queue)
    processed += 1
    for t in get(adj, s, Int[])
      in_degree[t] -= 1
      in_degree[t] == 0 && push!(queue, t)
    end
  end

  if processed < length(ch.constraints)
    cycle_dofs = sort!([line.dof for line in ch.constraints if in_degree[line.dof] > 0])
    error("Constraint cycle detected among DoFs: $cycle_dofs\n" *
        "Cycles are not resolvable by substitution. Check your constraint sources.")
  end
end

"""
    get_matrix(ch) → (C, g)

Return the sparse prolongation matrix `C` (n × n_free) and offset vector `g`
(length n) such that the full solution vector satisfies `x = C * u + g`,
where `u` is the reduced (free-dof) solution.

Columns of `C` correspond to free DoFs in natural index order.
Free-DoF rows are identity rows; slave-DoF rows carry the constraint coefficients.
"""
function Algebra.get_matrix(ch::ConstraintHandler{T}) where T
  @check ch.closed "The ConstraintHandler must be closed before building the constraint matrix."
  n      = ch.n_dofs
  n_free = num_free_dofs(ch)

  # Map each free DoF to its column index in C.
  col = 0
  dof_to_col = zeros(Int, n)
  for i in 1:n
    if !is_constrained(ch, i)
      col += 1
      dof_to_col[i] = col
    end
  end

  I, J, V = Int[], Int[], T[]
  g = zeros(T, n)

  for i in 1:n
    if !is_constrained(ch, i)
      push!(I, i)
      push!(J, dof_to_col[i])
      push!(V, one(T))
    else
      line = ch.constraints[ch.dof_to_constraint[i]]
      g[i] = line.offset
      for (master, coeff) in zip(line.masters, line.coeffs)
        push!(I, i)
        push!(J, dof_to_col[master])
        push!(V, coeff)
      end
    end
  end

  C = sparse(I, J, V, n, n_free)
  return C, g
end

"""
    apply_constraints(A, b, ch) → (Ã, b̃)

Return the reduced system `Ã = C'AC` and `b̃ = C'(b - Ag)` obtained by
eliminating all constrained DoFs via the prolongation matrix `C` and offset
vector `g` from `constraint_matrix(ch)`.

Works with any matrix type; sparsity is preserved when `A` is sparse.
After solving `Ã * u = b̃`, recover the full solution with `x = C * u + g`.
"""
function apply_constraints(A, b::AbstractVector, ch::ConstraintHandler)
  @check ch.closed "Must call close!(ch) before condense"
  C, g = get_matrix(ch)
  Ã = C' * A * C
  b̃ = C' * (b - A * g)
  return Ã, b̃
end

"""
    apply_constraints!(x, ch)

Overwrite `x[dof]` for every constrained DoF using the stored constraints:

  x[dof] = offset + Σ_j coefficient_j * x[master_j]

After `close!()`, all masters are free, so the iteration order is immaterial.
"""
function apply_constraints!(x::AbstractVector, ch::ConstraintHandler)
  @check ch.closed "Must call close!() before apply_constraints!()"
  @check length(x) == ch.n_dofs "Vector length $(length(x)) does not match n_dofs $(ch.n_dofs)"
  for line in ch.constraints
    val = line.offset
    for (master, coeff) in zip(line.masters, line.coeffs)
      val += coeff * x[master]
    end
    x[line.dof] = val
  end
  return x
end

# ---------------------------------------------------------------------------
# FESpaceWithLinearConstraints from a ConstraintHandler
# ---------------------------------------------------------------------------

function constraint_tables(V::SingleFieldFESpace, ch::ConstraintHandler{T}) where T
  @check ch.closed "Must call close!(ch) before building constraint tables"

  n_fdofs = num_free_dofs(V)
  n_s     = length(ch.constraints)

  ptrs = Vector{Int}(undef, n_s + 1)
  ptrs[1] = 1
  for (i, line) in enumerate(ch.constraints)
    ptrs[i+1] = ptrs[i] + length(line.masters)
  end
  nnz          = ptrs[end] - 1
  master_dofs  = Vector{Int}(undef, nnz)
  coefficients = Vector{T}(undef, nnz)
  for (i, line) in enumerate(ch.constraints)
    for j in eachindex(line.masters)
      master_dofs[ ptrs[i] + j - 1] = _DOF_to_dof(line.masters[j], n_fdofs)
      coefficients[ptrs[i] + j - 1] = line.coeffs[j]
    end
  end

  sDOF_to_dof    = [_DOF_to_dof(line.dof, n_fdofs) for line in ch.constraints]
  sDOF_to_dofs   = Table(master_dofs, ptrs)
  sDOF_to_coeffs = Table(coefficients, copy(ptrs))

  return sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs
end

"""
    FESpaceWithLinearConstraints(V::SingleFieldFESpace, ch::ConstraintHandler)

Build a `FESpaceWithLinearConstraints` from `V` and a closed `ConstraintHandler`.
"""
function FESpaceWithLinearConstraints(V::SingleFieldFESpace, ch::ConstraintHandler)
  sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = constraint_tables(V, ch)
  return FESpaceWithLinearConstraints(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V)
end
