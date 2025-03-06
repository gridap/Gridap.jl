module HHOConvgTests

  using Gridap
  using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
  using Gridap.CellData, Gridap.Fields, Gridap.Helpers
  using Gridap.ReferenceFEs

  function projection_operator(V, Ω, dΩ)
    Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
    mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
    P = FESpaces.LocalOperator(
      FESpaces.LocalSolveMap(), V, mass, mass; trian_out = Ω
    )
    return P
  end

  function reconstruction_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)

    Λ = FESpaces.PolytopalFESpace(Ω, Float64, 0; space=:P)

    nrel = get_normal_vector(Γp)
    Πn(v) = ∇(v)⋅nrel
    Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
    lhs((u,λ),(v,μ))   = ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
    rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
    
    mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
    W = MultiFieldFESpace([L,Λ];style=mfs)
    R = FESpaces.LocalOperator(
      FESpaces.HHO_ReconstructionOperatorMap(), ptopo, W, X, lhs, rhs; space_out = L
    )
    return R
  end

  function patch_dof_ids_w_bcs(ptopo,M::FESpace,MD::FESpace)
    ids_D = get_cell_dof_ids(MD)
    ids = get_cell_dof_ids(M)
    
    id_map = zeros(Int32,num_free_dofs(M))
    for (I,ID) in zip(ids,ids_D)
      id_map[I] = ID
    end
    
    patch_dofs = FESpaces.get_patch_assembly_ids(M,ptopo)
    patch_dofs_D = lazy_map(Broadcasting(Reindex(id_map)),patch_dofs)

    free_values = zero_free_values(MD)
    dirichlet_values = get_dirichlet_dof_values(MD)
    cell_values = lazy_map(Broadcasting(Gridap.Arrays.PosNegReindex(free_values,dirichlet_values)),patch_dofs_D)
    cell_mask = collect(Bool,lazy_map(dofs -> any(dof -> dof < 0, dofs),patch_dofs_D))

    return patch_dofs_D, cell_values, cell_mask
  end

  function patch_dof_ids_w_bcs(ptopo,X::MultiFieldFESpace,XD::MultiFieldFESpace)
    nfields = MultiField.num_fields(X)

    pids, cv, cm = [], [], []
    map(X,XD) do M, MD
      a, b, c = patch_dof_ids_w_bcs(ptopo,M,MD)
      push!(pids,a)
      push!(cv,b)
      push!(cm,c)
    end

    patch_dofs = lazy_map(BlockMap(nfields,collect(1:nfields)),pids...)
    cell_values = lazy_map(BlockMap(nfields,collect(1:nfields)),cv...)
    cell_mask = map(any,zip(cm...))

    return patch_dofs, cell_values, cell_mask
  end

  ##############################################################
  #   HHO: Mixed order variant
  ##############################################################
  function solve_Poisson_HHO(domain,nc,order,uex,f)
    model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,nc))
    vmodel = Gridap.Geometry.voronoi(simplexify(model))

    D = num_cell_dims(vmodel)
    Ω = Triangulation(ReferenceFE{D}, vmodel)
    Γ = Triangulation(ReferenceFE{D-1}, vmodel)

    ptopo = Geometry.PatchTopology(vmodel)
    Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
    Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

    qdegree = 2*(order+1)
    dΩp = Measure(Ωp,qdegree)
    dΓp = Measure(Γp,qdegree)

    V = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P) # Bulk space
    reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)  
    M = FESpace(Γ, reffe_M; conformity=:L2)                      # Skeleton space
    L = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P) # Reconstruction space

    mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
    X   = MultiFieldFESpace([V, M];style=mfs)

    PΓ = projection_operator(M, Γp, dΓp)
    R = reconstruction_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)

    M0 = TestFESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")  # For assembly and imposing boundary conditions
    MD = TrialFESpace(M0, uex)

    Y0 = MultiFieldFESpace([V, M0];style=mfs)
    YD = MultiFieldFESpace([V, MD];style=mfs)

    function a(u,v)
      Ru_Ω, Ru_Γ = R(u)
      Rv_Ω, Rv_Γ = R(v)
      return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
    end

    polys = get_polytopes(vmodel)
    hT = map(x -> FESpaces.get_facet_diameter(x,D), polys)

    hTinv = 1 ./ CellField(hT,Ωp)
    function s(u,v)
      function S(u)
        u_Ω, u_Γ = u
        return PΓ(u_Ω) - u_Γ
      end
      return ∫(hTinv * (S(u)⋅S(v)))dΓp
    end

    l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

    function patch_weakform(u,v)
      assem = FESpaces.PatchAssembler(ptopo,YD,YD)

      data = FESpaces.collect_patch_cell_matrix_and_vector(assem,YD,YD,s(u,v),l(v),zero(YD))
      matvecdata, matdata, vecdata = data

      dofs_a, cellvals, cellmask = patch_dof_ids_w_bcs(ptopo,X,YD)
      cell_matvec = attach_dirichlet(get_array(a(u,v)),cellvals,cellmask)
      p = Geometry.get_patch_faces(ptopo,num_cell_dims(vmodel))
      q = Base.OneTo(num_cells(vmodel))

      rows_a = FESpaces.attach_patch_map(FESpaces.AssemblyStrategyMap{:rows}(assem.strategy),[dofs_a],[q])[1]
      cols_a = FESpaces.attach_patch_map(FESpaces.AssemblyStrategyMap{:cols}(assem.strategy),[dofs_a],[q])[1]
      push!(matvecdata[1],cell_matvec)
      push!(matvecdata[2],rows_a)
      push!(matvecdata[3],cols_a)
      push!(matvecdata[4],p)

      data = (matvecdata,matdata,vecdata)
      return assemble_matrix_and_vector(assem,data)
    end

    function sc_assembly(u,v)
      full_matvecs = patch_weakform(u,v)
      sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)

      assem = FESpaces.PatchAssembler(ptopo,YD,YD)
      patch_rows = assem.strategy.array.array[2,2].patch_rows
      patch_cols = assem.strategy.array.array[2,2].patch_cols
      matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
      matdata = ([],[],[]) # dummy matdata
      vecdata = ([],[],[]) # dummy vecdata
      data = (matvecdata, matdata, vecdata)
      A, b = assemble_matrix_and_vector(SparseMatrixAssembler(MD,M0),data)
      return A, b
    end

    function backward_sc(xb)
      assem = FESpaces.PatchAssembler(ptopo,YD,YD)
      patch_rows_bb = assem.strategy.array.array[2,2].patch_rows
    
      patchwise_xb = map(patch_rows_bb) do patch_row
        patchwise_xb = xb[patch_row]
        return patchwise_xb
      end
    
      full_matvecs = patch_weakform(get_trial_fe_basis(YD),get_fe_basis(YD))
      patchwise_xi_vec = lazy_map(Gridap.FESpaces.BackwardStaticCondensationMap(),full_matvecs,patchwise_xb)
    
      patch_rows_ii = assem.strategy.array.array[1,1].patch_rows
      patch_cols_ii = assem.strategy.array.array[1,1].patch_cols
      patchwise_xi_vecdata = ([patchwise_xi_vec,],[patch_rows_ii,],[patch_cols_ii,])
    
      xi = assemble_vector(SparseMatrixAssembler(V,V), patchwise_xi_vecdata)
      wh = FEFunction(V,xi)
      return wh
    end  

    Asc, bsc = sc_assembly(get_trial_fe_basis(YD),get_fe_basis(YD))

    solver = LUSolver()
    ns = numerical_setup(symbolic_setup(solver,Asc),Asc)
    xb = zeros(size(bsc))
    solve!(xb,ns,bsc)
    uΓ = FEFunction(MD,xb)

    uΩ = backward_sc(xb)

    eu = uex - uΩ
    l2u= sqrt( sum(∫(eu ⋅ eu)dΩp) )
    h1u= l2u + sqrt( sum(∫(∇(eu) ⋅ ∇(eu))dΩp) )

    h = maximum(hT)

    return l2u, h1u, h
  end

  function convg_test(domain,ncs,order,uex,f)
      
    el2 = Float64[]
    eh1 = Float64[]
    hs = Float64[]
    for nc in ncs
      l2, h1, h = solve_Poisson_HHO(domain,nc,order,uex,f)
      push!(el2,l2)
      push!(eh1,h1)
      push!(hs,h)
    end
    println(el2)
    println(eh1)
    println(hs)
    el2, eh1, hs
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x), 1), x) \ y
    linreg[2]
  end

  uex(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
  f(x) = -Δ(uex)(x)

  domain = (0,1,0,1)
  ncs = [(2,2),(4,4),(8,8),(16,16),(32,32),(64,64),(128,128),(256,256)]
  order = 0

  el2s, eh1s, hs = convg_test(domain,ncs,order,uex,f)
  println("Slope L2-norm u: $(slope(hs,el2s))")
  println("Slope H1-norm u: $(slope(hs,eh1s))")

end # module