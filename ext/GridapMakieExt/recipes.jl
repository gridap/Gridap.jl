setup_color(color::Union{Symbol, Makie.Colorant}, ::Grid) = color

function setup_color(color::AbstractArray, grid::Grid)
    color = if length(color) == num_nodes(grid)
                to_dg_node_values(grid, color)
            elseif length(color) == num_cells(grid)
                to_dg_cell_values(grid, color)
            else
                @unreachable
            end
end

setup_face_color(color::Union{Symbol, Makie.Colorant}, ::Grid, ::Any) = color

function setup_face_color(color::AbstractArray, grid::Grid, face_to_cell)
    color = if length(color) == num_nodes(grid)
                color
            elseif length(color) == num_cells(grid)
                color[face_to_cell]
            else
                @unreachable
            end
end

@Makie.recipe MeshField begin
    linewidth = @inherit linewidth
    linecap = @inherit linecap
    joinstyle = @inherit joinstyle
    linestyle = nothing
    miter_limit = @inherit miter_limit
    cycle      = nothing
    uv_transform = Makie.automatic
    matcap = nothing
    marker = @inherit marker
    markersize = @inherit markersize
    strokecolor = @inherit markerstrokecolor
    strokewidth = @inherit markerstrokewidth
    glowcolor = (:black, 0.0)
    glowwidth = 0.0
    rotation = Makie.Billboard()
    marker_offset = Makie.Vec3f(0)
    font = @inherit markerfont
    distancefield = nothing
    interpolate = true
    font = "default"
    uv_offset_width = (0.0, 0.0, 0.0, 0.0)
    markerspace = :pixel
    depthsorting = false
    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_colormap_attributes()...
    Makie.mixin_shading_attributes()...
    fxaa = false
    shading    = Makie.NoShading
    colormap   = :bluesreds
    color = :pink
end

Makie.plottype(::Triangulation{<:Any,1}) = Makie.Scatter
Makie.plottype(::Triangulation{<:Any,1}, ::Any) = Makie.Scatter
Makie.plottype(::Triangulation) = MeshField
Makie.plottype(::Triangulation,::Any) = MeshField
Makie.plottype(c::CellField) = Makie.plottype(get_triangulation(c),c)
Makie.args_preferred_axis(::Triangulation{<:Any,Dp}) where Dp = Dp<=2 ? Makie.Axis : Makie.LScene

function Makie.convert_arguments(::Type{<:Makie.Wireframe}, trian::Triangulation{<:Any,1})
    x=to_point1D(trian)
    return (x[1],x[2],zeros(length(x[1]),length(x[1])))
end

function Makie.convert_arguments(::Type{<:Makie.Wireframe}, trian::Triangulation)
    grid = to_grid(trian)
    mesh = to_plot_mesh(grid)
    (mesh, )
end

function Makie.convert_arguments(::Type{<:Makie.Scatter}, trian::Triangulation)
    grid = to_grid(trian)
    node_coords = get_node_coordinates(grid)
    x = map(to_point, node_coords)
    (x, )
end

function Makie.convert_arguments(::Type{<:MeshField}, c::CellField)
    (get_triangulation(c), c)
end

function Makie.convert_arguments(t::Union{Type{Makie.Lines},Type{Makie.ScatterLines},Type{Makie.Scatter}}, c::CellField)
    Makie.convert_arguments(t,get_triangulation(c), c)
end

function Makie.convert_arguments(::Union{Type{Makie.Lines},Type{Makie.ScatterLines},Type{Makie.Scatter}}, trian::Triangulation{<:Any,1}, uh::Any)
    return to_point1D(trian, uh)
end

function Makie.convert_arguments(::Union{Type{Makie.Lines},Type{Makie.ScatterLines},Type{Makie.Scatter}}, trian::Triangulation{<:Any,1})
    return to_point1D(trian)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{0,<:Any}, Any}})
    map!(p.attributes,[:converted_1,:converted_2],[:newcolor,:mesh]) do trian,uh
        grid_and_data=to_grid(trian,uh)
        mesh=to_plot_dg_mesh(grid_and_data[1])
        color0=grid_and_data[2]
        color=setup_color(color0,grid_and_data[1])
        return (color,mesh)
    end
    Makie.scatter!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{1,<:Any}, Any}})
    map!(p.attributes,[:converted_1,:converted_2],[:newcolor,:mesh]) do trian,uh
        grid_and_data=to_grid(trian,uh)
        mesh=to_plot_dg_mesh(grid_and_data[1])
        color0=grid_and_data[2]
        color=setup_color(color0,grid_and_data[1])
        return (color,mesh)
    end
    Makie.linesegments!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{2,<:Any}, Any}})
    map!(p.attributes,[:converted_1,:converted_2],[:newcolor,:mesh]) do trian,uh
        grid_and_data=to_grid(trian,uh)
        mesh=to_plot_dg_mesh(grid_and_data[1])
        color0=grid_and_data[2]
        color=setup_color(color0,grid_and_data[1])
        return (color,mesh)
    end
    Makie.mesh!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{3,<:Any}, Any}})
    map!(p.attributes,[:converted_1,:converted_2],[:newcolor,:mesh]) do trian,uh
        grid_and_data=to_grid(trian,uh)
        color0=grid_and_data[2]
        grid = grid_and_data[1]
        face_grid_and_map = to_face_grid_with_map(grid)
        face_grid = face_grid_and_map[1]
        face_to_cell = face_grid_and_map[2]
        face_color = setup_face_color(color0, grid, face_to_cell)
        mesh = face_grid |> to_plot_dg_mesh |> GeometryBasics.normal_mesh
        color = setup_color(face_color, face_grid)
        return (color,mesh)
    end
    Makie.mesh!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{0,<:Any}}})
    map!(p.attributes,[:converted_1,:color],[:newcolor,:mesh]) do trian,c
        grid=to_grid(trian)
        mesh=to_plot_dg_mesh(grid)
        color=setup_color(c,grid)
        return (color,mesh)
    end
    Makie.scatter!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{1,<:Any}}})
    map!(p.attributes,[:converted_1,:color],[:newcolor,:mesh]) do trian,c
        grid=to_grid(trian)
        mesh=to_plot_dg_mesh(grid)
        color=setup_color(c,grid)
        return (color,mesh)
    end
    Makie.linesegments!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{2,<:Any}}})
    map!(p.attributes,[:converted_1,:color],[:newcolor,:mesh]) do trian,c
        grid=to_grid(trian)
        mesh=to_plot_dg_mesh(grid)
        color=setup_color(c,grid)
        return (color,mesh)
    end
    Makie.mesh!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(p::MeshField{<:Tuple{Triangulation{3,<:Any}}})
    map!(p.attributes,[:converted_1,:color],[:newcolor,:mesh]) do trian,c
        grid=to_grid(trian)
        face_grid_and_map = to_face_grid_with_map(grid)
        face_grid_and_map = to_face_grid_with_map(grid)
        face_grid = face_grid_and_map[1]
        face_to_cell = face_grid_and_map[2]
        face_color = setup_face_color(c, grid, face_to_cell)
        mesh = face_grid |> to_plot_dg_mesh |> GeometryBasics.normal_mesh
        color = setup_color(face_color, face_grid)
        return (color,mesh)
    end
    Makie.mesh!(p, p.attributes, p.mesh,color=p.newcolor)
end

function Makie.plot!(::MeshField{<:Tuple{Triangulation, Any}})
    @unreachable
end

function Makie.plot!(::MeshField{<:Tuple{Triangulation}})
    @unreachable
end

