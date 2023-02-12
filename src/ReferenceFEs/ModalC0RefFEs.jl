struct ModalC0 <: ReferenceFEName end

const modalC0 = ModalC0()

"""
  ModalC0RefFE(::Type{T},p::Polytope{D},orders) where {T,D}

Returns an instance of `GenericRefFE{ModalC0}` representing a ReferenceFE with
Modal C0-continuous shape functions (multivariate scalar-valued, vector-valued,
or tensor-valued, iso- or aniso-tropic).

For more details about the shape functions, see Section 1.1.5 in

Ern, A., & Guermond, J. L. (2013). Theory and practice of finite elements
(Vol. 159). Springer Science & Business Media.

and references therein.

The constructor is only implemented for for n-cubes and the minimum order in
any dimension must be greater than one. The DoFs are numbered by n-faces in the
same way as with CLagrangianRefFEs.
"""
function ModalC0RefFE(
  ::Type{T},
  p::Polytope{D},
  orders,
  a::Vector{Point{D,V}},
  b::Vector{Point{D,V}};
  space::Symbol=_default_space(p) ) where {T,D,V}

  @notimplementedif ! is_n_cube(p)
  @notimplementedif minimum(orders) < one(eltype(orders))

  shapefuns = ModalC0Basis{D}(T,orders,a,b)

  ndofs, predofs, lag_reffe, face_dofs = compute_reffe_data(T,p,orders,space=space)

  GenericRefFE{ModalC0}(
    ndofs,
    p,
    predofs,
    GradConformity(),
    lag_reffe,
    face_dofs,
    shapefuns)
end

function ModalC0RefFE(
  ::Type{T},
  p::Polytope{D},
  orders;
  space::Symbol=_default_space(p) ) where {T,D}

  @notimplementedif ! is_n_cube(p)
  @notimplementedif minimum(orders) < one(eltype(orders))

  shapefuns = ModalC0Basis{D}(T,orders)

  ndofs, predofs, lag_reffe, face_dofs = compute_reffe_data(T,p,orders,space=space)

  GenericRefFE{ModalC0}(
    ndofs,
    p,
    predofs,
    GradConformity(),
    lag_reffe,
    face_dofs,
    shapefuns)
end

function get_orders(reffe::GenericRefFE{ModalC0,D}) where{D}
  get_orders(reffe.prebasis)
end

function compute_reffe_data(::Type{T},
                            p::Polytope{D},
                            order::Int;
                            space::Symbol=_default_space(p)) where {T,D}
  orders = tfill(order,Val{D}())
  compute_reffe_data(T,p,orders,space=space)
end

function compute_reffe_data(::Type{T},
                            p::Polytope{D},
                            orders::NTuple{D,Int};
                            space::Symbol=_default_space(p)) where {T,D}
  lag_reffe = LagrangianRefFE(T,p,orders,space=space)
  reffe = lag_reffe.reffe
  reffe.ndofs, reffe.dofs, lag_reffe, reffe.face_dofs
end

function ReferenceFE(
  polytope::Polytope{D},
  ::ModalC0,
  ::Type{T},
  orders::Union{Integer,NTuple{D,Int}};
  kwargs...) where {T,D}
  ModalC0RefFE(T,polytope,orders;kwargs...)
end

function Conformity(::GenericRefFE{ModalC0},sym::Symbol)
  h1 = (:H1,:C0,:Hgrad)
  if sym == :L2
    L2Conformity()
  elseif sym in h1
    H1Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a ModalC0RefFE with H1 conformity.
    Possible values of conformity for this reference fe are $((:L2, h1...)).
    """
  end
end

function get_face_own_dofs(
  reffe::GenericRefFE{ModalC0},conf::GradConformity)
  lagrangian_reffe = reffe.metadata
  get_face_own_dofs(lagrangian_reffe,conf)
end

function get_face_own_dofs_permutations(
  reffe::GenericRefFE{ModalC0},conf::GradConformity)
  lagrangian_reffe = reffe.metadata
  get_face_own_dofs_permutations(lagrangian_reffe,conf)
end

function compute_shapefun_bboxes!(
  a::Vector{Point{D,V}},
  b::Vector{Point{D,V}},
  bboxes::Vector{Point{D,V}},
  face_own_dofs) where {D,V}
  for i in 1:length(face_own_dofs)
    a[face_own_dofs[i]] .= bboxes[2*i-1]
    b[face_own_dofs[i]] .= bboxes[2*i]
  end
end

function compute_cell_to_modalC0_reffe(
  p::Polytope{D},
  ncells::Int,
  ::Type{T},
  orders::Union{Integer,NTuple{D,Int}},
  bboxes;
  space::Symbol=_default_space(p)) where {T,D} # type-stability?

  @notimplementedif ! is_n_cube(p)
  @notimplementedif minimum(orders) < one(eltype(orders))
  @assert ncells == length(bboxes)

  ndofs, predofs, lag_reffe, face_dofs = compute_reffe_data(T,p,orders,space=space)
  face_own_dofs = get_face_own_dofs(lag_reffe,GradConformity())

  filter = space == :Q ? _q_filter : _s_filter_mc0

  sh(bbs) = begin
    a = fill(Point{D,eltype(T)}(tfill(zero(eltype(T)),Val{D}())),ndofs)
    b = fill(Point{D,eltype(T)}(tfill(zero(eltype(T)),Val{D}())),ndofs)
    compute_shapefun_bboxes!(a,b,bbs,face_own_dofs)
    ModalC0Basis{D}(T,orders,a,b,filter=filter)
  end

  reffe(sh) = GenericRefFE{ModalC0}(ndofs,
                                    p,
                                    predofs,
                                    GradConformity(),
                                    lag_reffe,
                                    face_dofs,
                                    sh)

  reffes = [ reffe(sh(bbs)) for bbs in bboxes ]
  CompressedArray(reffes,1:ncells)
end

function compute_cell_to_modalC0_reffe(
  p::Polytope{D},
  ncells::Int,
  ::Type{T},
  orders::Union{Integer,NTuple{D,Int}};
  space::Symbol=_default_space(p)) where {T,D} # type-stability?

  @notimplementedif ! is_n_cube(p)
  @notimplementedif minimum(orders) < one(eltype(orders))

  filter = space == :Q ? _q_filter : _s_filter_mc0

  ndofs, predofs, lag_reffe, face_dofs = compute_reffe_data(T,p,orders,space=space)
  reffe = GenericRefFE{ModalC0}(ndofs,
                                p,
                                predofs,
                                GradConformity(),
                                lag_reffe,
                                face_dofs,
                                ModalC0Basis{D}(T,orders,filter=filter))

  Fill(reffe,ncells)
end