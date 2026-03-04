struct Strang <: QuadratureName end
const strang = Strang()

# Partially adapted from
# https://github.com/FEniCS/basix/blob/main/cpp/basix/quadrature.cpp

function Quadrature(
  p::Polytope, ::Strang, degree::Integer; T::Type{<:AbstractFloat}=Float64
)
  if p == TRI
    wx = _strang_quad_tri(degree; T)
    polyname = "triangle"
  elseif p == TET
    wx = _strang_quad_tet(degree; T)
    polyname = "tetrahedron"
  else
    msg = """\n
    `strang` quadrature rule only available for triangles and tetrahedra.
    Use `duffy` for other simplices, and `tensor_product` for n-cubes.
    """
    error(msg)
  end
  coords, weights = _weightcoords_to_coord_weights(p, wx, T)
  name = "Strang quadrature of degree $degree on the reference $(polyname)."
  GenericQuadrature(coords, weights, name)
end

function maxdegree(p::Polytope, ::Strang)
  if p == TRI
    11
  elseif p == TET
    5
  else
    0
  end
end

#######
# TRI #
#######
function _strang_quad_tri(degree; T::Type{<:AbstractFloat}=Float64)
  maxdegree = 11
  if degree > maxdegree
    msg = """\n
    `strang` quadrature rule only available up to degree $(maxdegree) on a triangle.
    Use `duffy` instead.
    """
    error(msg)
  end

  data1 = zero(T) # w (1/3, 1/3)
  data2 = Tuple{T,T,T}[] # w (s, t) (t, s) (t, t)
  data3 = Tuple{T,T,T,T}[] # w (s, t) (t, s) (s, u) (u, s) (t, u) (u, t)

  if degree in (0, 1)
    data1 = 1 / 2
  elseif degree == 2
    data2 = [
      (1 / 6, 2 / 3, 1 / 6)
    ]
  elseif degree == 3
    data1 = -27 / 96
    data2 = [
      (25 / 96, 3 / 5, 1 / 5)
    ]
  elseif degree == 4
    data2 = [
      (0.054975871827661, 0.816847572980459, 0.091576213509771),
      (0.1116907948390055, 0.108103018168070, 0.445948490915965),
    ]
  elseif degree == 5
    s = sqrt(15)
    data1 = 0.1125
    data2 = [
      ((155 - s) / 2400, (9 + 2 * s) / 21, (6 - s) / 21),
      ((155 + s) / 2400, (9 - 2 * s) / 21, (6 + s) / 21),
    ]
  elseif degree == 6
    data2 = [
      (0.0254224531851035, 0.873821971016996, 0.063089014491502)
      (0.0583931378631895, 0.501426509658179, 0.249286745170910)
    ]
    data3 = [
      (0.041425537809187, 0.053145049844816, 0.310352451033785, 0.636502499121399)
    ]
  elseif degree == 7
    data1 = -0.149570044467670 / 2
    data2 = [
      (0.175615257433204 / 2, 0.479308067841923, 0.260345966079038),
      (0.053347235608839 / 2, 0.869739794195568, 0.065130102902216)
    ]
    data3 = [
      (0.077113760890257 / 2, 0.638444188569809, 0.312865496004875, 0.048690315425316)
    ]
  elseif degree in (8, 9)
    data1 = 0.09713579628279610 / 2
    data2 = [
      (0.03133470022713983 / 2, 0.02063496160252593, 0.4896825191987370),
      (0.07782754100477543 / 2, 0.1258208170141290, 0.4370895914929355),
      (0.07964773892720910 / 2, 0.6235929287619356, 0.1882035356190322),
      (0.02557767565869810 / 2, 0.9105409732110941, 0.04472951339445297),
    ]
    data3 = [
      (0.04328353937728940 / 2, 0.7411985987844980, 0.03683841205473626, 0.22196298916076573)
    ]
  elseif degree in (10, 11)
    data1 = 0.08797730116222190 / 2.0
    data2 = [
      (0.008744311553736190 / 2, 0.9480217181434233, 0.02598914092828833),
      (0.03808157199393533 / 2, 0.8114249947041546, 0.09428750264792270),
      (0.01885544805613125 / 2, 0.01072644996557060, 0.4946367750172147),
      (0.07215969754474100 / 2, 0.5853132347709715, 0.2073433826145142),
      (0.06932913870553720 / 2, 0.1221843885990187, 0.4389078057004907),
    ]
    data3 = [
      (0.04105631542928860 / 2, 0.6779376548825902, 0.04484167758913055, 0.27722066752827925),
      (0.007362383783300573 / 2, 0.8588702812826364, 0.0, 0.1411297187173636)
    ]
  else
    @unreachable
  end

  n = !iszero(data1) + 3 * length(data2) + 6 * length(data3)
  wx = Array{Float64,2}(undef, n, 3)

  r = 1
  if !iszero(data1)
    wx[r, 1] = data1
    wx[r, 2] = 1 / 3
    wx[r, 3] = 1 / 3
    r += 1
  end
  for (w, s, t) in data2
    for (x, y) in ((s, t), (t, s), (t, t))
      wx[r, 1] = w
      wx[r, 2] = x
      wx[r, 3] = y
      r += 1
    end
  end
  for (w, s, t, u) in data3
    for (x, y) in ((s, t), (t, s), (s, u), (u, s), (t, u), (u, t))
      wx[r, 1] = w
      wx[r, 2] = x
      wx[r, 3] = y
      r += 1
    end
  end

  wx
end

#######
# TET #
#######
function _strang_quad_tet(degree; T::Type{<:AbstractFloat}=Float64)
  maxdegree = 5
  if degree > maxdegree
    msg = """\n
    `strang` quadrature rule only available up to degree $(maxdegree) on a tetrahedron.
    Use `duffy` instead.
    """
    error(msg)
  end

  data1 = zero(T) # w (1/4, 1/4, 1/4)
  data2 = Tuple{T,T,T}[] # w (s, t, t) (t, s, t) (t, t, s) (t, t, t)
  data3 = Tuple{T,T,T}[] # w (s, t, t) (t, s, t) (t, t, s) (t, s, s) (s, t, s) (s, s, t)

  if degree in (0, 1)
    data1 = 1 / 6
  elseif degree == 2
    data2 = [
      (1 / 24, 0.5854101966249685, 0.1381966011250105)
    ]
  elseif degree == 3
    data1 = -2 / 15
    data2 = [
      (3 / 40, 1 / 2, 1 / 6)
    ]
  elseif degree == 4
    data1 = (-148 / 1875) / 6
    data2 = [
      ((343 / 7500) / 6, 11 / 14, 1 / 14)
    ]
    data3 = [
      ((56 / 375) / 6, 0.3994035761667992, 0.1005964238332008)
    ]
  elseif degree == 5
    data2 = [
      (0.1126879257180162 / 6, 0.0673422422100983, 0.3108859192633005),
      (0.0734930431163619 / 6, 0.7217942490673264, 0.0927352503108912),
    ]
    data3 = [
      (0.0425460207770812 / 6, 0.4544962958743506, 0.0455037041256494)
    ]
  else
    @unreachable
  end

  n = !iszero(data1) + 4 * length(data2) + 6 * length(data3)
  wx = Array{Float64,2}(undef, n, 4)

  r = 1
  if !iszero(data1)
    wx[r, 1] = data1
    wx[r, 2] = 1 / 4
    wx[r, 3] = 1 / 4
    wx[r, 4] = 1 / 4
    r += 1
  end
  for (w, s, t) in data2
    for (x, y, z) in ((s, t, t), (t, s, t), (t, t, s), (t, t, t))
      wx[r, 1] = w
      wx[r, 2] = x
      wx[r, 3] = y
      wx[r, 4] = z
      r += 1
    end
  end
  for (w, s, t) in data3
    for (x, y, z) in ((s, t, t), (t, s, t), (t, t, s), (t, s, s), (s, t, s), (s, s, t))
      wx[r, 1] = w
      wx[r, 2] = x
      wx[r, 3] = y
      wx[r, 4] = z
      r += 1
    end
  end

  wx
end
