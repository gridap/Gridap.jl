# This source file is though to put that code required in order to
# customize ConformingFESpaces.jl to H(div)-conforming global FE Spaces built
# out of RaviartThomas FEs. In particular, this customization is in the
# definition of the shape functions (get_shapefuns) and the DoFs
# (get_dof_basis) of the **global** FE space, which requires a sign flip
# for those sitting on facets of the slave cell of the facet.

function get_dof_basis(model::DiscreteModel,
                       cell_reffe::AbstractArray{<:GenericRefFE{RaviartThomas}})
    sign_flip = get_sign_flip(model, cell_reffe)
    lazy_map(get_dof_basis,
             cell_reffe,
             get_cell_map(Triangulation(model)),
             sign_flip)
end

function get_shapefuns(model::DiscreteModel,
                       cell_reffe::AbstractArray{<:GenericRefFE{RaviartThomas}})
    sign_flip = get_sign_flip(model, cell_reffe)
    lazy_map(get_shapefuns,
             cell_reffe,
             get_cell_map(Triangulation(model)),
             lazy_map(Broadcasting(ConstantField), sign_flip))
end

function get_sign_flip(model::DiscreteModel,
                       cell_reffe::AbstractArray{<:GenericRefFE{RaviartThomas}})
    lazy_map(get_sign_flip,
            Fill(model, length(cell_reffe)),
            cell_reffe,
            get_cell_to_bgcell(model))
end

function get_sign_flip(model::DiscreteModel,
                       reffe::GenericRefFE{RaviartThomas},
                       cell_id::Integer)

    cache = return_cache(get_sign_flip,
                         reffe,
                         cell_id)
    evaluate!(cache,
              get_sign_flip,
              reffe,
              cell_id)
end

function return_cache(::typeof(get_sign_flip),
                      model::DiscreteModel,
                      reffe::GenericRefFE{RaviartThomas},
                      cell_id::Integer)
    D = num_cell_dims(model)
    gtopo = get_grid_topology(model)

    # Extract composition among cells and facets
    cell_wise_facets_ids = get_faces(gtopo, D, D - 1)
    cache_cell_wise_facets_ids = array_cache(cell_wise_facets_ids)

    # Extract cells around facets
    cells_around_facets = get_faces(gtopo, D - 1, D)
    cache_cells_around_facets = array_cache(cells_around_facets)

    (cell_wise_facets_ids,
     cache_cell_wise_facets_ids,
     cells_around_facets,
     cache_cells_around_facets,
     CachedVector(Bool))
end

function evaluate!(cache,
                   ::typeof(get_sign_flip),
                  model::DiscreteModel,
                  reffe::GenericRefFE{RaviartThomas},
                  cell_id::Integer)

    cell_wise_facets_ids,
    cache_cell_wise_facets_ids,
    cells_around_facets,
    cache_cells_around_facets,
    sign_flip = cache

    setsize!(sign_flip, (num_dofs(reffe),))
    sign_flip .= false

    D = num_dims(reffe)
    face_own_dofs = get_face_own_dofs(reffe)
    facet_lid = get_offsets(get_polytope(reffe))[D] + 1
    cell_facets_ids = getindex!(cache_cell_wise_facets_ids,
                                cell_wise_facets_ids,
                                cell_id)
    for facet_gid in cell_facets_ids
        facet_cells_around = getindex!(cache_cells_around_facets,
             cells_around_facets,
             facet_gid)
        is_slave = (findfirst((x) -> (x == cell_id), facet_cells_around) == 2)
        if is_slave
            for dof in face_own_dofs[facet_lid]
                sign_flip[dof] = true
            end
        end
        facet_lid = facet_lid + 1
    end
    sign_flip
end
