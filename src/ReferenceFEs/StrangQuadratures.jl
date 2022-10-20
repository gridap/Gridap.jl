
struct Strang <: QuadratureName end
const strang = Strang()

function Quadrature(p::Polytope,::Strang,degree::Integer;T::Type{<:AbstractFloat}=Float64)
  msg = """\n
  `strang` quadrature rule only available for simplices.
  Use `tensor_product` for n-cubes.
  """
  @assert is_simplex(p) msg
  if num_dims(p) == 2
    quad = _strang_quad_tri(degree;T=T)
  elseif num_dims(p) == 3
    quad = _strang_quad_tet(degree,T=T)
  else
    msg = """\n
    `strang` quadrature rule only available for tris and tets.
    Use `duffy` for other simplices.
    """
    @unreachable msg
  end
  quad
end

const _strang_tri_k2n = Dict(1=>1,2=>3,3=>4,4=>6,5=>7,7=>13,9=>19,11=>28)

function _strang_quad_tri(degree;T::Type{<:AbstractFloat}=Float64)
  if ! haskey(_strang_tri_k2n,degree)
    msg = """\n
      `strang` quadrature rule not implemented for degree = $degree on a triangle.
      Implemented degrees are in $(sort(collect(keys(_strang_tri_k2n)))).
      Use `duffy` instead.
    """
    error(msg)
  end
  n = _strang_tri_k2n[degree]
  x = Vector{VectorValue{2,T}}(undef,n)
  w = Vector{T}(undef,n)
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
    b = (9.0+2.0*sqrt(15.0))/21.0
    c = (6.0-sqrt(15.0))/21.0
    d = (9.0-2.0*sqrt(15.0))/21.0
    e = (6.0+sqrt(15.0))/21.0
    w1 = 0.1125
    w2 = (155.0-sqrt(15.0))/2400.0
    w3 = (155.0+sqrt(15.0))/2400.0
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
    h = 0.8588702812826364
    z = 0.0
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
    x[23] = (h,z)
    x[24] = (h,y)
    x[25] = (z,h)
    x[26] = (z,y)
    x[27] = (y,h)
    x[28] = (y,z)
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
  GenericQuadrature(x,w,"Strang quadrature of degree $degree on the reference triangle.")
end

const _strang_tet_k2n = Dict(1=>1,2=>4,3=>5,4=>11,5=>14)

function _strang_quad_tet(degree;T::Type{<:AbstractFloat}=Float64)
  if ! haskey(_strang_tet_k2n,degree)
    msg = """\n
      `strang` quadrature rule not implemented for degree = $degree on a tet.
      Implemented degrees are in $(sort(collect(keys(_strang_tet_k2n)))).
      Use `duffy` instead.
    """
    error(msg)
  end
  n = _strang_tet_k2n[degree]
  x = Vector{VectorValue{3,T}}(undef,n)
  w = Vector{T}(undef,n)
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
  GenericQuadrature(x,w,"Strang quadrature of degree $degree on the reference tet.")
end
