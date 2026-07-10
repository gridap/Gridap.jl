using Gridap, Gridap.MultiField
using FiniteDiff
using Test


model = CartesianDiscreteModel((0,1,0,1),(8,8))
order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

V = FESpace(model,reffe;dirichlet_tags="boundary",vector_type=Vector{ComplexF64})

xh = interpolate(x -> x[1]*x[2],V)

using Gridap.FESpaces, Gridap.CellData
using Gridap.FESpaces: _change_argument, _compute_cell_ids, GridapADTag, autodiff_array_gradient, get_cell_dof_values, get_domains, lazy_map, add_contribution!, get_free_dof_values
using Gridap.CellData: get_ad_level

function _change_argument_real(f,cell_u)
  y = lazy_map(imag,cell_u)
  g_at_complex = x -> f(lazy_map((xi,yi) -> xi + im * yi, x, y))
  f₁ = real ∘ g_at_complex
  f₂ = imag ∘ g_at_complex
  f₁,f₂
end

function Gridap.FESpaces._gradient(f,uh,fuh::DomainContribution;tag::GridapADTag=get_ad_level(fuh)+1)
  terms = DomainContribution(;ad_level = tag)
  for trian in get_domains(fuh)
    g = _change_argument(gradient,f,trian,uh)
    cell_id = _compute_cell_ids(uh,trian)
    cell_u = get_cell_dof_values(uh)
    # We compute the derivative of f as f'(u) = ∂ᵣf₁ + i*∂ᵣf₂ via Cauchy-Riemann equations
    # where f(u) = f₁(r,s) + i*f₂(r,s) and u = r + i*s. I.e., in this case we perturb only the
    # real part and keep the imaginary part fixed.
    r = lazy_map(real,cell_u)
    f₁,f₂ = _change_argument_real(g,cell_u)
    ∂ᵣf₁ = autodiff_array_gradient(f₁,r,cell_id;tag)
    ∂ᵣf₂ = autodiff_array_gradient(f₂,r,cell_id;tag)
    add_contribution!(terms,trian, lazy_map((u,v)-> u + im*v,∂ᵣf₁,∂ᵣf₂))
  end
  terms
end

_dj = gradient(u->∫(u)dΩ,xh)
_dj_vec = assemble_vector(_dj,V)
_dj_vec_analytic = assemble_vector(v-> ∫(v)dΩ,V)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(u*u)dΩ,xh)
_dj_vec = assemble_vector(_dj,V)
_dj_vec_analytic = assemble_vector(v-> ∫(2*xh*v)dΩ,V)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(im*u*u)dΩ,xh)
_dj_vec = assemble_vector(_dj,V)
_dj_vec_analytic = assemble_vector(v-> ∫(2im*xh*v)dΩ,V)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(u + im*u*u)dΩ,xh)
_dj_vec = assemble_vector(_dj,V)
_dj_vec_analytic = assemble_vector(v-> ∫(v + 2im*xh*v)dΩ,V)
_dj_vec≈_dj_vec_analytic

xh2 = interpolate(x -> im*x[1]*x[2],V)
_dj = gradient(u->∫(u)dΩ,xh2)
_dj_vec = assemble_vector(_dj,V)
_dj_vec_analytic = assemble_vector(v-> ∫(v)dΩ,V)
_dj_vec≈_dj_vec_analytic

xh2 = interpolate(x -> im*x[1]*x[2],V)
_dj = gradient(u->∫(im*u)dΩ,xh2)
_dj_vec = assemble_vector(_dj,V)
_dj_vec_analytic = assemble_vector(v-> ∫(im*v)dΩ,V)
_dj_vec≈_dj_vec_analytic

xh2 = interpolate(x -> x[1] + im*x[1]*x[2],V)
_dj = gradient(u->∫(u + im*u*u)dΩ,xh2)
_dj_vec = assemble_vector(_dj,V)
_dj_vec_analytic = assemble_vector(v-> ∫(v + 2im*xh2*v)dΩ,V)
_dj_vec≈_dj_vec_analytic