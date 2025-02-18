module PolytopalHDGConvgTests

  using Gridap
  using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
  using Gridap.CellData, Gridap.Fields, Gridap.Helpers


  function statically_condensed_assembly(ptopo,X,Y,M,M_test,a,l)
    # Lazily assemble the global patch-systems and 
    # perform the static condensation
    assem = FESpaces.PatchAssembler(ptopo,X,Y)
    full_matvecs = assemble_matrix_and_vector(a,l,assem,X,Y)
    sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)

    # Regular assembly of the statically-assembled systems
    patch_rows = assem.strategy.array.array[2,2].patch_rows
    patch_cols = assem.strategy.array.array[2,2].patch_cols
    matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
    matdata = ([],[],[]) # dummy matdata
    vecdata = ([],[],[]) # dummy vecdata
    data = (matvecdata, matdata, vecdata)
    A, b = assemble_matrix_and_vector(SparseMatrixAssembler(M,M_test),data)
    return A, b
  end

  function backward_static_condensation(ptopo,X,Y,a,l,xb)
    assem = FESpaces.PatchAssembler(ptopo,X,Y)
    patch_rows_bb = assem.strategy.array.array[2,2].patch_rows
    patch_cols_bb = assem.strategy.array.array[2,2].patch_cols

    # To discuss: can patch_rows and patch_cols be different from each other for a patch?
    if patch_rows_bb != patch_cols_bb
      @notimplemented
    else
      patchwise_xb = map(patch_rows_bb) do patch_row
        patchwise_xb = xb[patch_row]
        return patchwise_xb
      end
    end

    full_matvecs = assemble_matrix_and_vector(a,l,assem,X,Y)
    patchwise_xi_vec = lazy_map(Gridap.FESpaces.BackwardStaticCondensationMap(),full_matvecs,patchwise_xb)

    patch_rows_ii = assem.strategy.array.array[1,1].patch_rows
    patch_cols_ii = assem.strategy.array.array[1,1].patch_cols
    patchwise_xi_vecdata = ([patchwise_xi_vec,],[patch_rows_ii,],[patch_cols_ii,])

    Xi = MultiFieldFESpace([X[1],X[2]])
    Yi = MultiFieldFESpace([Y[1],Y[2]])
    assem = SparseMatrixAssembler(Xi,Yi)
    xi = assemble_vector(assem, patchwise_xi_vecdata)
    wh = FEFunction(Yi,xi)
    return wh
  end

  function solve_Poisson_HDG(domain,nc,order,u,f)
    model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,nc))
    vmodel = Gridap.Geometry.voronoi(simplexify(model))

    D = num_cell_dims(vmodel)
    Ω = Triangulation(ReferenceFE{D}, vmodel)
    Γ = Triangulation(ReferenceFE{D-1}, vmodel)
    polys = get_polytopes(vmodel)
    h = maximum( map(p->FESpaces.get_facet_diameter(p,D),polys) )

    ptopo = Geometry.PatchTopology(vmodel)
    Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
    Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

    # HDG bulk spaces
    V = FESpaces.PolytopalFESpace(Ω, VectorValue{D, Float64}, order; space=:P)
    Q = FESpaces.PolytopalFESpace(Ω, Float64, order; space=:P)
    
    # HDG skeleton space(s)
    reffeM = ReferenceFE(lagrangian, Float64, order; space=:P) 
    M_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")
    M = TrialFESpace(M_test, u)

    mfs = MultiField.BlockMultiFieldStyle(2,(2,1))
    Y = MultiFieldFESpace([V, Q, M_test];style=mfs)
    X = MultiFieldFESpace([V, Q, M];style=mfs)

    τ = 1.0 # HDG stab parameter

    degree = 2*(order+1)
    dΩp = Measure(Ωp,degree)
    dΓp = Measure(Γp,degree)
    n = get_normal_vector(Γp)
    
    Πn(u) = u⋅n
    Π(u) = change_domain(u,Γp,DomainStyle(u))
    a((qh,uh,sh),(vh,wh,lh)) = ∫( qh⋅vh - uh*(∇⋅vh) - qh⋅∇(wh) )dΩp + ∫(sh*Πn(vh))dΓp +
                              ∫((Πn(qh) + τ*(Π(uh) - sh))*(Π(wh) + lh))dΓp
    l((vh,wh,hatmh)) = ∫( f*wh )*dΩp

    Asc, bsc = statically_condensed_assembly(ptopo,X,Y,M,M_test,a,l)

    solver = LUSolver()
    ns = numerical_setup(symbolic_setup(solver,Asc),Asc)
    xb = zeros(size(bsc))
    solve!(xb,ns,bsc)
    sh = FEFunction(M_test,xb)

    xi = backward_static_condensation(ptopo,X,Y,a,l,xb)
    qh, uh = xi

    eu = u - uh
    l2u = sqrt( sum( ∫( eu ⋅ eu )dΩp ) ) # TO discuss: Is this integration consistent?

    return l2u, h
  end

  function convg_test(domain,ncs,order,u,f)
    
    el2 = Float64[]
    hs = Float64[]
    for nc in ncs
      l2u, h = solve_Poisson_HDG(domain,nc,order,u,f)
      println(l2u)
      push!(el2,l2u)
      push!(hs,h)
    end
    println(el2)
    println(hs)
    el2, hs
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x), 1), x) \ y
    linreg[2]
  end


  u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
  q(x) = -∇(u)(x)
  f(x) = (∇ ⋅ q)(x)

  domain = (0,1,0,1)
  ncs = [(2,2),(4,4),(8,8),(16,16),(32,32),(64,64),(128,128)]
  order = 4

  el, hs = convg_test(domain,ncs,order,u,f)
  println("Slope L2-norm u: $(slope(hs,el))")

end # module