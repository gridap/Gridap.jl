module TestGridapMakieExt

using Gridap
using GLMakie
using Test

import FileIO

const OUTDIR = joinpath(@__DIR__, "output")
rm(OUTDIR, force=true, recursive=true)
mkpath(OUTDIR)

function savefig(f, suffix::String)
    fig = f()
    println("*"^80)
    filename = "$(suffix).png"
    @show filename
    path = joinpath(OUTDIR, filename)
    FileIO.save(path, fig)
    return true
end

@testset "Tests 1D" begin
    model = CartesianDiscreteModel((0,15),90)
    Ω = Triangulation(model)
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = FESpace(model, reffe)
    u=x->cos(π*x[1])*exp(-x[1]/10)
    uh = interpolate(u, V)
    Γ = BoundaryTriangulation(model)

    @test savefig("1d_Fig1") do
        fig = lines(Ω,u)
        scatter!(uh,color=:red)
        fig
    end
    @test savefig("1d_Fig2") do
        fig=scatterlines(uh)
        plot!(Γ,uh,color=:red)
        fig
    end
    @test savefig("1d_Fig3") do
        fig=plot(Ω,color=:black)
        plot!(Γ,color=:red,marker=:xcross,markersize=15)
        fig
    end

end

@testset "Tests 2D" begin

    domain = (0,1,0,1)
    cell_nums = (10,10)
    model = CartesianDiscreteModel(domain,cell_nums) |> simplexify
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)
    n_Λ = get_normal_vector(Λ)
    reffe = ReferenceFE(lagrangian,Float64,1)
    V = FESpace(model,reffe)
    uh = interpolate(x->sin(π*(x[1]+x[2])),V)
    celldata = rand(num_cells(Ω))

    @test savefig("2d_Fig1") do
        fig = plot(Ω)
        wireframe!(Ω, linestyle=:dash, color=:black)
        plot!(Γ,color=:blue,linewidth=5)
        scatter!(Ω, markersize=10, marker=:diamond, color=:brown2)
        fig
    end
    @test savefig("2d_Fig2") do
        fig, _ , sc = plot(Ω, uh)
        Colorbar(fig[1,2], sc)
        fig
    end
    @test savefig("2d_Fig3") do
        fig,_,sc = plot(Γ, uh, colormap=:algae, linewidth=10)
        Colorbar(fig[1,2],sc)
        fig
    end
    @test savefig("2d_Fig4") do
        fig = plot(Ω, color=:green)
        wireframe!(Ω, color=:red, linewidth=2.5)
        fig
    end
    @test savefig("2d_Fig5") do
        fig, _ , plt = plot(Ω, color=3*celldata, colormap=:heat)
        Colorbar(fig[1,2], plt)
        fig
    end
    @test savefig("2d_Fig6") do
        fig, _ , plt = plot(Λ, jump(n_Λ⋅∇(uh)),linewidth=4)
        Colorbar(fig[1,2], plt)
        fig
    end
    @test savefig("2d_fig7") do 
        fig, _ , plt = plot(uh, colormap=:Spectral)
        Colorbar(fig[1,2], plt)
        fig   
    end
end

@testset "Tests 3D" begin

    domain = (0,1,0,1,0,1)
    cell_nums = (10,10,10)
    model = CartesianDiscreteModel(domain,cell_nums) |> simplexify
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)
    n_Λ = get_normal_vector(Λ)
    reffe = ReferenceFE(lagrangian,Float64,1)
    V = FESpace(model,reffe)
    uh = interpolate(x->sin(π*(x[1]+x[2]+x[3])),V)
    celldata = rand(num_cells(Ω))

    @test savefig("3d_Fig1") do
        fig = plot(Ω)
        wireframe!(Ω)
        scatter!(Ω)
        fig
    end
    @test savefig("3d_Fig2") do
        fig, _ , sc = plot(Ω, uh, colorrange=(0,1))
        Colorbar(fig[1,2],sc)
        fig
    end
    @test savefig("3d_Fig3") do
        fig, _ , sc = plot(Γ, uh, colormap=:algae)
        Colorbar(fig[1,2],sc)
        fig
    end
    @test savefig("3d_Fig4") do
        fig = plot(Ω, color=:green)
        wireframe!(Ω, color=:red, linewidth=2.5)
        fig
    end
    @test savefig("3d_Fig5") do
        fig, _ , plt = plot(Ω, color=3*celldata, colormap=:heat)
        Colorbar(fig[1,2], plt)
        fig
    end
    @test savefig("3d_Fig6") do
        fig, _ , plt = plot(Λ, jump(n_Λ⋅∇(uh)))
        Colorbar(fig[1,2],plt)
        fig
    end
    @test savefig("3d_fig7") do 
        fig, _ , plt = plot(uh, colormap=:Spectral, colorrange=(0,1))
        Colorbar(fig[1,2], plt)
        fig   
    end
end

end #module
