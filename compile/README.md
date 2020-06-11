The Julia script available in this folder (`create_gridap_image.jl`) lets one to create a custom sysimage of `Gridap.jl` using the so-called [`PackageCompiler.jl`](https://github.com/JuliaLang/PackageCompiler.jl) Julia package. See [documentation page](https://julialang.github.io/PackageCompiler.jl/dev/sysimages/) for more details.

The Julia script at hand provides a set of command-line-arguments (CLAs) in order to customize the creation process of the sysimage. You can see the list of options available using the following command:

```bash
$ julia create_gridap_image.jl -h 
```

In order to create the sysimage, the script clones the Git repository of `Gridap.jl` and that of its [Tutorials](https://github.com/gridap/Tutorials), and (currently) runs the `validation.jl` tutorial. By default, the script will generate the sysimage from the code corresponding to the last commit pushed to the `master` branch of both Git repositories. However, via appropriate CLAs (see help message), one may select the particular version of both projects to be cloned by specifying the corresponding Git release tags. For example, if one wants to generate a sysimage of Gridap `v0.10.4` from Tutorials `v0.10.0`, one has to execute the following command:


```bash
$ julia create_gridap_image.jl -g v0.10.4 -t v0.10.0 -n Gridapv0.10.4.so
```

Obviously, one has to select a version of the Tutorials package compatible with the version selected for Gridap, **otherwise the image creation process will fail**. Note that in the above command, we are also specifying the name of the image to be generated, which by default, is created in the directory from which the script is executed
(the path where the image is generated can be customized as well via the corresponding CLA).

Once the image is generated, one may use it to execute a Julia script that uses Gridap as follows:

```bash
$ julia -J Gridapv0.10.4.so gridap_script.jl
```

or run an interactive Julia session as:

```bash
$ julia -J Gridapv0.10.4.so
```

### Known issues, TO-DO, misc learned lessons, etc.

1. **The sysimage works and is useful when `Gridap.jl` is being used but not being modified**, e.g., while developing `GridapDistributed.jl` or `GridapODEs.jl`. If you want to develop `Gridap.jl`, use [`Revise.jl`](https://github.com/gridap/Gridap.jl/wiki/REPL-based-workflow#editing-small-parts-of-the-project).  

2. Also note the following from the [documentation page](https://julialang.github.io/PackageCompiler.jl/dev/sysimages/) of `PackageCompiler.jl`: *It should be clearly stated that there are some drawbacks to using a custom sysimage, thereby sidestepping the standard Julia package precompilation system. The biggest drawback is that packages that are compiled into a sysimage (including their dependencies!) are "locked" to the version they where at when the sysimage was created. This means that no matter what package version you have installed in your current project, the one in the sysimage will take precedence. This can lead to bugs where you start with a project that needs a specific version of a package, but you have another one compiled into the sysimage.*

3. We have to decide which code to execute in order to create the image. At present, `validation.jl` is chosen just for testing purposes. Ideally it should force the compilation of a large part of the project, but the problem that it solves should take small amount of time to execute once the code is compiled. 

4. We have observed that the resulting image does not contain all the code whose compilation is triggered from `validation.jl`. The code that is compiled during execution time can be spotted with
    ```bash
    $ julia --trace-compile=stderr -q -J Gridapv0.10.4.so validation.jl
    ```
    At present, we have the following:
    ```
    $ julia --trace-compile=stderr ../../Tutorials/src/validation.jl 2>&1 | wc -l 
    2920   
    $ julia --trace-compile=stderr -J./Gridap.so ../../Tutorials/src/validation.jl 2>&1 | wc -l
    825
    ```
   An hypothesis (not confirmed): the list of packages provided to `PackageCompiler.create_sysimage` should not only contain direct dependencies of `Gridap.jl`. This is related to the following comment in the script:
    ```
    # TODO: reorg project structure (Project.toml) to use
        # Pkg.dependencies() with PackageCompiler
        #if VERSION >= v"1.4"
        #    append!(pkgs, [Symbol(v.name) for v in values(Pkg.dependencies())])
        #end
    ```
    `Pkg.dependencies()` is a brand new method in Julia `v1.4`. From the resulting data structure that it provides one may obtain all the nodes of the dependency graph of the currently loaded project.

5. See also issue [#276](https://github.com/gridap/Gridap.jl/issues/276).




