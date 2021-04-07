"""
    struct DuffyQuadrature{D,T} <: Quadrature{D,T}
      coordinates::Vector{Point{D,T}}
      weights::Vector{T}
    end

Duffy quadrature for simplices in [0,1]^D
"""
struct DuffyQuadrature{D,T} <: Quadrature{D,T}
  coordinates::Vector{Point{D,T}}
  weights::Vector{T}
end

get_coordinates(q::DuffyQuadrature) = q.coordinates

get_weights(q::DuffyQuadrature) = q.weights


"""
    DuffyQuadrature{D}(degree::Integer) where D
"""
function DuffyQuadrature{D}(order::Integer) where D
  x,w = _duffy_quad_data(order,D)
  DuffyQuadrature(x,w)
end

function _duffy_quad_data(order::Integer,D::Int)

  beta = 0
  dim_to_quad_1d = [
    _gauss_jacobi_in_0_to_1(order,(D-1)-(d-1),beta) for d in 1:(D-1) ]

  quad_1d = _gauss_legendre_in_0_to_1(order)
  push!(dim_to_quad_1d,quad_1d)

  x_pos = 1
  w_pos = 2
  dim_to_xs_1d = [quad_1d[x_pos] for quad_1d in dim_to_quad_1d]
  dim_to_ws_1d = [quad_1d[w_pos] for quad_1d in dim_to_quad_1d]

  a = 0.5
  for d in (D-1):-1:1
    ws_1d = dim_to_ws_1d[d]
    ws_1d[:] *= a
    a *= 0.5
  end

  x,w = _tensor_product_duffy(dim_to_xs_1d,dim_to_ws_1d)

  (_duffy_map.(x),w)

end

# Duffy map from the n-cube in [0,1]^d to the n-simplex in [0,1]^d
function _duffy_map(q::Point{D,T}) where {D,T}
  m = zero(Mutable(Point{D,T}))
  m[1] = q[1]
  a = one(T)
  for i in 2:D
    a *= (1-q[i-1])
    m[i] = a*q[i]
  end
  Point{D,T}(m)
end

_duffy_map(q::Point{1,T}) where T = q

function _gauss_jacobi_in_0_to_1(order,alpha,beta)
  n = _npoints_from_order(order)
  x,w = gaussjacobi(n,alpha,beta)
  _map_to(0,1,x,w)
end

function _gauss_legendre_in_0_to_1(order)
  n = _npoints_from_order(order)
  x,w = gausslegendre(n)
  _map_to(0,1,x,w)
end

# Transforms a 1-D quadrature from `[-1,1]` to `[a,b]`, with `a<b`.
function _map_to(a,b,points,weights)
  points_ab = 0.5*(b-a)*points .+ 0.5*(a+b)
  weights_ab = 0.5*(b-a)*weights
  (points_ab, weights_ab)
end

function _npoints_from_order(order)
  ceil(Int, (order + 1.0) / 2.0 )
end

function _tensor_product_duffy(
  dim_to_xs_1d::Vector{Vector{T}},
  dim_to_ws_1d::Vector{Vector{W}}) where {T,W}

  D = length(dim_to_ws_1d)
  @assert D == length(dim_to_xs_1d)
  dim_to_n = [length(ws_1d) for ws_1d in dim_to_ws_1d]
  n = prod(dim_to_n)
  xs = zeros(Point{D,T},n)
  ws = zeros(W,n)
  cis = CartesianIndices(tuple(dim_to_n...))
  m = zero(Mutable(Point{D,T}))
  _tensor_product_duffy!(xs,ws,dim_to_xs_1d,dim_to_ws_1d,cis,m)
  (xs,ws)
end

function _tensor_product_duffy!(
  xs,ws,dim_to_xs_1d,dim_to_ws_1d,cis::CartesianIndices{D},m) where D
  k = 1
  for ci in cis
    w = 1.0
    for d in 1:D
      xs_1d = dim_to_xs_1d[d]
      ws_1d = dim_to_ws_1d[d]
      i = ci[d]
      xi = xs_1d[i]
      wi = ws_1d[i]
      w *= wi
      m[d] = xi
    end
    xs[k] = m
    ws[k] = w
    k += 1
  end
end

# Hard-coded quadratures for simplices

function simplex_quadrature(::Val{D},degree) where D
  if D == 2
    tri_quadrature(degree)
  elseif D == 3
    tet_quadrature(degree)
  else
    GenericQuadrature(DuffyQuadrature{D}(degree))
  end
end

function tri_quadrature(degree)
  GenericQuadrature(DuffyQuadrature{2}(degree))
end

function _symmetric_gauss_legendre_tri(degree)
  k2n = Dict(1=>1,2=>3,3=>4,4=>6,5=>7,7=>13,9=>19,11=>28)
  if ! haskey(k2n,degree)
    msg = """\n
      _symmetric_gauss_legendre_tri not implemented for degree = $degree
      Implemented degrees are in $(collect(keys(k2n))).
    """
    error(msg)
  end
  n = k2n[degree]
  x = Vector{VectorValue{2,Float64}}(undef,n)
  w = Vector{Float64}(undef,n)
  if degree == 1
    a = 1.0/3.0
    b = 1.0/2.0
    x[1] = (a,a)
    w[1] = b
  elseif degree == 2
    a = 2.0/3.0
    b = 1.0/6.0
    x[1] = (b,b)
    x[2] = (a,b)
    x[3] = (b,a)
    w[1] = b
    w[2] = b
    w[3] = b
  elseif degree == 3
    a = 1.0/3.0
    b = 1.0/5.0
    c = 3.0/5.0
    d = -27.0/96.0
    e = 25.0/96.0
    x[1] = (a,a)
    x[2] = (b,b)
    x[3] = (c,b)
    x[4] = (b,c)
    w[1] = d
    w[2] = e
    w[3] = e
    w[4] = e
  elseif degree == 4
    ex1 = 0.816847572980459
    et1 = 0.091576213509771
    ez1 = 0.091576213509771
    ex2 = 0.108103018168070
    et2 = 0.445948490915965
    ez2 = 0.445948490915965
    a = 0.054975870996713638
    b = 0.1116907969117165    
    x[1] = (et1,ez1)
    x[2] = (ez2,ex2)
    x[3] = (ex1,et1)
    x[4] = (ex2,et2)
    x[5] = (et2,ez2)
    x[6] = (ez1,ex1)
    w[1]  = a
    w[2]  = b
    w[3]  = a
    w[4]  = b
    w[5]  = b
    w[6]  = a
  elseif degree == 5
    a = 1.0 / 3.0
    b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0
    c = ( 6.0 -          sqrt ( 15.0 ) ) / 21.0
    d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0
    e = ( 6.0 +          sqrt ( 15.0 ) ) / 21.0
    w1 = 0.1125
    w2 = ( 155.0 - sqrt ( 15.0 ) ) / 2400.0
    w3 = ( 155.0 + sqrt ( 15.0 ) ) / 2400.0
    x[1] = (a,a)
    x[2] = (b,c)
    x[3] = (c,b)
    x[4] = (c,c)
    x[5] = (d,e)
    x[6] = (e,d)
    x[7] = (e,e)
    w[1] = w1
    w[2] = w2
    w[3] = w2
    w[4] = w2
    w[5] = w3
    w[6] = w3
    w[7] = w3
  elseif degree == 7
    a = 0.333333333333333
    b = 0.479308067841923
    c = 0.869739794195568
    d = 0.638444188569809
    e = 0.260345966079038
    f = 0.065130102902216
    g = 0.312865496004875
    h = 0.048690315425316
    w1=-0.149570044467670/2.0
    w2= 0.175615257433204/2.0
    w3= 0.053347235608839/2.0
    w4= 0.077113760890257/2.0
    x[ 1] = (a,a)
    x[ 2] = (e,e)
    x[ 3] = (b,e)        
    x[ 4] = (e,b)        
    x[ 5] = (f,f)        
    x[ 6] = (c,f)        
    x[ 7] = (f,c)        
    x[ 8] = (d,g)        
    x[ 9] = (d,h)        
    x[10] = (g,d)        
    x[11] = (g,h)        
    x[12] = (h,d)        
    x[13] = (h,g)
    w[ 1] = w1
    w[ 2] = w2
    w[ 3] = w2
    w[ 4] = w2
    w[ 5] = w3
    w[ 6] = w3
    w[ 7] = w3
    w[ 8] = w4
    w[ 9] = w4
    w[10] = w4
    w[11] = w4
    w[12] = w4
    w[13] = w4
  elseif degree == 9
    a = 1.0 / 3.0
    b = 0.02063496160252593
    c = 0.4896825191987370
    d = 0.1258208170141290
    e = 0.4370895914929355
    f = 0.6235929287619356
    g = 0.1882035356190322
    r = 0.9105409732110941
    s = 0.04472951339445297
    t = 0.7411985987844980
    u = 0.03683841205473626
    v = 0.22196298916076573
    w1 = 0.09713579628279610/2.0
    w2 = 0.03133470022713983/2.0
    w3 = 0.07782754100477543/2.0
    w4 = 0.07964773892720910/2.0
    w5 = 0.02557767565869810/2.0
    w6 = 0.04328353937728940/2.0
    x[ 1] = (a,a)
    x[ 2] = (b,c)
    x[ 3] = (c,b)
    x[ 4] = (c,c)
    x[ 5] = (d,e)
    x[ 6] = (e,d)
    x[ 7] = (e,e)
    x[ 8] = (f,g)
    x[ 9] = (g,f)
    x[10] = (g,g)
    x[11] = (r,s)
    x[12] = (s,r)
    x[13] = (s,s)
    x[14] = (t,u)
    x[15] = (t,v)
    x[16] = (u,t)
    x[17] = (u,v)
    x[18] = (v,t)
    x[19] = (v,u)
    w[ 1] = w1
    w[ 2] = w2
    w[ 3] = w2
    w[ 4] = w2
    w[ 5] = w3
    w[ 6] = w3
    w[ 7] = w3
    w[ 8] = w4
    w[ 9] = w4
    w[10] = w4
    w[11] = w5
    w[12] = w5
    w[13] = w5
    w[14] = w6
    w[15] = w6
    w[16] = w6
    w[17] = w6
    w[18] = w6
    w[19] = w6
  elseif degree == 11
    a = 1.0 / 3.0
    b = 0.9480217181434233
    c = 0.02598914092828833
    d = 0.8114249947041546
    e = 0.09428750264792270
    f = 0.01072644996557060
    g = 0.4946367750172147
    p = 0.5853132347709715
    q = 0.2073433826145142
    r = 0.1221843885990187
    s = 0.4389078057004907
    t = 0.6779376548825902
    u = 0.04484167758913055
    v = 0.27722066752827925
    w = 0.8588702812826364
    x = 0.0
    y = 0.1411297187173636
    w1 = 0.08797730116222190/2.0
    w2 = 0.008744311553736190/2.0
    w3 = 0.03808157199393533/2.0
    w4 = 0.01885544805613125/2.0
    w5 = 0.07215969754474100/2.0
    w6 = 0.06932913870553720/2.0
    w7 = 0.04105631542928860/2.0
    w8 = 0.007362383783300573/2.0
    x[ 1] = (a,a)
    x[ 2] = (b,c)
    x[ 3] = (c,b)
    x[ 4] = (c,c)
    x[ 5] = (d,e)
    x[ 6] = (e,d)
    x[ 7] = (e,e)
    x[ 8] = (f,g)
    x[ 9] = (g,f)
    x[10] = (g,g)
    x[11] = (p,q)
    x[12] = (q,p)
    x[13] = (q,q)
    x[14] = (r,s)
    x[15] = (s,r)
    x[16] = (s,s)
    x[17] = (t,u)
    x[18] = (t,v)
    x[19] = (u,t)
    x[20] = (u,v)
    x[21] = (v,t)
    x[22] = (v,u)
    x[23] = (w,x)
    x[24] = (w,y)
    x[25] = (x,w)
    x[26] = (x,y)
    x[27] = (y,w)
    x[28] = (y,x)
    w[ 1] = w1
    w[ 2] = w2
    w[ 3] = w2
    w[ 4] = w2
    w[ 5] = w3
    w[ 6] = w3
    w[ 7] = w3
    w[ 8] = w4
    w[ 9] = w4
    w[10] = w4
    w[11] = w5
    w[12] = w5
    w[13] = w5
    w[14] = w6
    w[15] = w6
    w[16] = w6
    w[17] = w7
    w[18] = w7
    w[19] = w7
    w[20] = w7
    w[21] = w7
    w[22] = w7
    w[23] = w8
    w[24] = w8
    w[25] = w8
    w[26] = w8
    w[27] = w8
    w[28] = w8
  else
    @unreachable
  end
  GenericQuadrature(x,w)


if(ngaus==1) then
     else if(ngaus==3) then
     else if(ngaus==4) then
     else if(ngaus==6) then
     else if(ngaus==7) then
     else if(ngaus==13) then
     else if(ngaus==19) then
     else if(ngaus==28) then
     else
        write(*,*) __FILE__,__LINE__,'ERROR:: Quadrature not defined',ndime,ngaus
        check(.false.)
     end if
end

function _symmetric_gauss_legendre_tet(degree)
  k2n = Dict(1=>1,2=>4,3=>5,4=>11,5=>14)
  if ! haskey(k2n,degree)
    msg = """\n
      _symmetric_gauss_legendre_tet not implemented for degree = $degree
      Implemented degrees are in $(collect(keys(k2n))).
    """
    error(msg)
  end
  n = k2n[degree]
  x = Vector{VectorValue{3,Float64}}(undef,n)
  w = Vector{Float64}(undef,n)
  if degree==1
    a = 1.0/4.0
    b = 1.0/6.0
    x[1] = (a,a,a)
    w[1] = b
  elseif degree==2
    a = 0.5854101966249685
    b = 0.1381966011250105
    c = 1.0/24.0
    x[1] = (b,b,b)
    x[2] = (a,b,b)
    x[3] = (b,a,b)
    x[4] = (b,b,a)
    w[1] = c
    w[2] = c
    w[3] = c
    w[4] = c
  elseif degree==3
    a = 1.0/4.0
    b = 1.0/6.0
    c = 1.0/2.0
    d = -2.0/15.0
    e = 1.5/20.0
    x[1] = (a,a,a)
    x[2] = (b,b,b)
    x[3] = (c,b,b)
    x[4] = (b,c,b)
    x[5] = (b,b,c)
    w[1] = d
    w[2] = e
    w[3] = e
    w[4] = e
    w[5] = e
  elseif degree== 4
    a = 0.3994035761667992
    b = 0.1005964238332008
    c = (343.0/7500.0)/6.0
    d = (56.0/375.0)/6.0
    e = 1.0/4.0
    f = 11.0/14.0
    g = 1.0/14.0
    h = (-148.0/1875.0)/6.0
    x[ 1] = (e,e,e)
    x[ 2] = (f,g,g)
    x[ 3] = (g,f,g)
    x[ 4] = (g,g,f)
    x[ 5] = (g,g,g)
    x[ 6] = (a,a,b)
    x[ 7] = (a,b,a)
    x[ 8] = (a,b,b)
    x[ 9] = (b,a,a)
    x[10] = (b,a,b)
    x[11] = (b,b,a)
    w[ 1] = h
    w[ 2] = c
    w[ 3] = c
    w[ 4] = c
    w[ 5] = c
    w[ 6] = d
    w[ 7] = d
    w[ 8] = d
    w[ 9] = d
    w[10] = d
    w[11] = d
  elseif degree == 5
    a = 0.0673422422100983
    b = 0.3108859192633005
    c = 0.7217942490673264
    d = 0.0927352503108912
    e = 0.4544962958743506
    f = 0.0455037041256494
    p = 0.1126879257180162/6.0
    q = 0.0734930431163619/6.0
    r = 0.0425460207770812/6.0
    x[ 1] = (a,b,b)
    x[ 2] = (b,a,b)
    x[ 3] = (b,b,a)
    x[ 4] = (b,b,b)
    x[ 5] = (c,d,d)
    x[ 6] = (d,c,d)
    x[ 7] = (d,d,c)
    x[ 8] = (d,d,d)
    x[ 9] = (e,e,f)
    x[10] = (e,f,e)
    x[11] = (e,f,f)
    x[12] = (f,e,e)
    x[13] = (f,e,f)
    x[14] = (f,f,e)
    w[ 1] = p
    w[ 2] = p
    w[ 3] = p
    w[ 4] = p
    w[ 5] = q
    w[ 6] = q
    w[ 7] = q
    w[ 8] = q
    w[ 9] = r
    w[10] = r
    w[11] = r
    w[12] = r
    w[13] = r
    w[14] = r
  else
    @unreachable
  end
  GenericQuadrature(x,w)
end



