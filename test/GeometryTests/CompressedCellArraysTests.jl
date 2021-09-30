module CompressedCellArraysTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields

u(x) = x[1]^2 + x[2]

domain= (0,1,0,1)
cells = (3,3)
model = CartesianDiscreteModel(domain,cells)

Ω = Triangulation(model)
cell_xs = get_cell_coordinates(Ω)
cell_mask = lazy_map(cell_xs) do xs
  R = 0.7
  n = length(xs)
  x = (1/n)*sum(xs)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end
Ω1 = Triangulation(model,cell_mask)
Γ = Boundary(model)

tcell_num = [ 1 for cell in 1:num_cells(Γ)]
tcell_mat = [ones(3,3) for cell in 1:num_cells(Γ)]
tcell_vec = [ones(3) for cell in 1:num_cells(Γ)]
tcell_matvec = pair_arrays(tcell_mat,tcell_vec)
tcell_block = lazy_map(BlockMap(2,2),tcell_mat)

mcell_num = move_contributions(tcell_num,Γ,Ω)
@test length(mcell_num) == num_cells(Ω)
@test sum(mcell_num) == sum(tcell_num)

mcell_mat = move_contributions(tcell_mat,Γ,Ω)
@test length(mcell_mat) == num_cells(Ω)
@test sum(mcell_mat[length.(mcell_mat) .!= 0]) == sum(tcell_mat)

mcell_vec = move_contributions(tcell_vec,Γ,Ω)
@test length(mcell_vec) == num_cells(Ω)
@test sum(mcell_vec[length.(mcell_vec) .!= 0]) == sum(tcell_vec)

mcell_matvec = move_contributions(tcell_matvec,Γ,Ω)
@test length(mcell_matvec) == num_cells(Ω)
mcell_mat,mcell_vec = unpair_arrays(mcell_matvec)
@test sum(mcell_mat[length.(mcell_mat) .!= 0]) == sum(tcell_mat)
@test sum(mcell_vec[length.(mcell_vec) .!= 0]) == sum(tcell_vec)

mcell_block = move_contributions(tcell_block,Γ,Ω)
@test length(mcell_block) == num_cells(Ω)
mcell_mat = lazy_map(b->b.array[2],mcell_block)
@test length(mcell_mat) == num_cells(Ω)
@test sum(mcell_mat[length.(mcell_mat) .!= 0]) == sum(tcell_mat)

tcell_num = [ 1 for cell in 1:num_cells(Ω1)]
tcell_mat = [ones(3,3) for cell in 1:num_cells(Ω1)]
tcell_vec = [ones(3) for cell in 1:num_cells(Ω1)]
tcell_matvec = pair_arrays(tcell_mat,tcell_vec)
tcell_block = lazy_map(BlockMap(2,2),tcell_mat)

mcell_num = move_contributions(tcell_num,Ω1,Ω)
@test length(mcell_num) == num_cells(Ω)
@test sum(mcell_num) == sum(tcell_num)

mcell_mat = move_contributions(tcell_mat,Ω1,Ω)
@test length(mcell_mat) == num_cells(Ω)
@test sum(mcell_mat[length.(mcell_mat) .!= 0]) == sum(tcell_mat)

mcell_vec = move_contributions(tcell_vec,Ω1,Ω)
@test length(mcell_vec) == num_cells(Ω)
@test sum(mcell_vec[length.(mcell_vec) .!= 0]) == sum(tcell_vec)

mcell_matvec = move_contributions(tcell_matvec,Ω1,Ω)
@test length(mcell_matvec) == num_cells(Ω)
mcell_mat,mcell_vec = unpair_arrays(mcell_matvec)
@test sum(mcell_mat[length.(mcell_mat) .!= 0]) == sum(tcell_mat)
@test sum(mcell_vec[length.(mcell_vec) .!= 0]) == sum(tcell_vec)

mcell_block = move_contributions(tcell_block,Ω1,Ω)
@test length(mcell_block) == num_cells(Ω)
mcell_mat = lazy_map(b->b.array[2],mcell_block)
@test length(mcell_mat) == num_cells(Ω)
@test sum(mcell_mat[length.(mcell_mat) .!= 0]) == sum(tcell_mat)

end # module
