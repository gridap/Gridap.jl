#!/usr/bin/env julia
#
# Examples:  
#      julia create_gridap_image.jl -h 
#      julia create_gridap_image.jl -g v0.10.4

using Pkg
Pkg.add("ArgParse")
Pkg.add("PackageCompiler")
using ArgParse
using LibGit2
using PackageCompiler
const do_not_clone_gridap_flag="--do-not-clone-gridap"
const do_not_clone_tutorials_flag="--do-not-clone-tutorials"

function parse_commandline()
    s = ArgParseSettings()
    do_not_clone_gridap_flag="--do-not-clone-gridap"
    do_not_clone_tutorials_flag="--do-not-clone-tutorials"
    @add_arg_table! s begin
        "--image-name", "-n"
        help = "The name of the Gridap.jl custom system image that will be created"
        arg_type = String
        default = "Gridap.so"
        "--image-path", "-p"
        help = "The relative or absolute PATH where the Gridap.jl custom system image will be created"
        arg_type = String
        default = "./"
        "--gridap-tag", "-g"
        help = """Gridap.jl Git repo Tag string with the source code from which the Gridap.jl custom system 
                  image will be created. This option is ignored if $(do_not_clone_gridap_flag) flag IS PASSED.
               """
        arg_type = String
        default = "master"
        "--tutorials-tag", "-t"
        help = "Tutorials Git repo Tag string with the base code line that will be executed in order to generate 
                the Gridap.jl custom system image. This option is ignored if $(do_not_clone_tutorials_flag) flag IS PASSED"
        arg_type = String
        default = "master"
        "--gridap-path"
        help = """If $(do_not_clone_gridap_flag) flag IS NOT PASSED, the relative or absolute PATH where to 
                  clone the Gridap.jl Git repo (Warning: Removed if it exists!).
                  If $(do_not_clone_gridap_flag) flag IS PASSED, the relative or absolute PATH where an existing
                  Gridap.jl source directory tree can be found.
               """
        arg_type = String
        default = "/tmp/Gridap.jl/"
        "--tutorials-path"
        arg_type = String
        default = "/tmp/Tutorials/"
        help = """If $(do_not_clone_tutorials_flag) flag IS NOT PASSED, the relative or absolute PATH where to 
                  clone the Tutorials Git repo (Warning: Removed if it exists!).
                  If $(do_not_clone_tutorials_flag) flag IS PASSED, the relative or absolute PATH where an existing
                  Tutorials source directory tree can be found.
               """
        "--do-not-clone-gridap"
        help = "Do not clone the Gridap.jl Git repo, but instead use an existing source directory."
        action = :store_true
        "--do-not-clone-tutorials"
        help = "Do not clone the Tutorials Git repo, but instead use an existing source directory."
        action = :store_true
    end
    return parse_args(s)
end

function clone_and_checkout_tag(
    repo_url::AbstractString,
    repo_path::AbstractString,
    git_tag::AbstractString,
)
    rm(repo_path; force = true, recursive = true)
    repo = LibGit2.clone(repo_url, repo_path)
    if (lowercase(git_tag) != "master")
        @assert git_tag in LibGit2.tag_list(repo) "Wrong Git Tag $(git_tag). Available choices for $(repo_url) are $(LibGit2.tag_list(repo))"
        gt=LibGit2.GitTag(repo, git_tag)
        commit_hash=LibGit2.target(gt)
        LibGit2.checkout!(repo, string(commit_hash))
    end
end


function main()
    parsed_args = parse_commandline()
    image_name = parsed_args["image-name"]
    image_path = parsed_args["image-path"]
    gridap_tag = parsed_args["gridap-tag"]
    tutorials_tag = parsed_args["tutorials-tag"]
    gridap_path = parsed_args["gridap-path"]
    tutorials_path = parsed_args["tutorials-path"]
    clone_and_checkout_tag(
        "https://github.com/gridap/Gridap.jl",
        gridap_path,
        gridap_tag,
    )
    clone_and_checkout_tag(
       "https://github.com/gridap/Tutorials",
       tutorials_path,
       tutorials_tag,
    )

    @info "Creating system image for Gridap.jl#$(gridap_tag) object file at: '$(image_path)'"
    @info "Building Gridap.jl#$(gridap_tag) into system image: $(joinpath(image_path,image_name))"

    start_time = time()

    Pkg.activate(gridap_path)
    Pkg.instantiate(verbose = true)

    pkgs = Symbol[]
    push!(pkgs, :Gridap)
    #if climatemachine_pkg
    #    push!(pkgs, :ClimateMachine)
    #else
        # TODO: reorg project structure (Project.toml) to use
        # Pkg.dependencies() with PackageCompiler
        #if VERSION >= v"1.4"
        #    append!(pkgs, [Symbol(v.name) for v in values(Pkg.dependencies())])
        #end
        append!(pkgs, [Symbol(name) for name in keys(Pkg.installed())])
    #end

    # use package compiler
    PackageCompiler.create_sysimage(
        pkgs,
        sysimage_path = joinpath(image_path,image_name),
        precompile_execution_file = joinpath(
            tutorials_path,
            "src",
            "validation.jl",
        ),
    )

    tot_secs = Int(floor(time() - start_time))
    @info "Created system image object file at: $(joinpath(image_path,image_name))"
    @info "System image build time: $tot_secs sec"
end

main()
