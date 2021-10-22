module DensifyInnerMostBlockLevelMapsTests
   using Test
   using Gridap
   using Gridap.Fields
   using Gridap.Arrays

   m=DensifyInnerMostBlockLevelMap()

   # VectorBlock{Vector} test
   y=Vector{Vector{Float64}}(undef,3)
   touched=Vector{Bool}(undef,3)
   touched .= true
   x=rand(3)
   y[1]=x
   y[2]=x.+1
   y[3]=x.+2
   by=ArrayBlock(y,touched)
   cache=return_cache(m,by)
   @test vcat(y...)==evaluate!(cache,m,by)

   # MatrixBlock{Vector} test
   x=rand(3)
   y=Array{Vector{Float64},2}(undef,(1,4))
   y[1,1]=x
   y[1,3]=x.+2
   y[1,4]=x.+3
   touched  = Array{Bool,2}(undef,(1,4))
   touched .= true
   touched[1,2]=false
   by=ArrayBlock(y,touched)
   cache=return_cache(m,by)
   @test hcat(y[1,1],zeros(3),y[1,3],y[1,4]) == evaluate!(cache,m,by)

   # VectorBlock{Matrix} test
   x=rand(1,3)
   y=Array{Matrix{Float64}}(undef,(4,))
   y[1]=x
   y[2]=x.+1
   y[3]=x.+2
   y[4]=x.+3
   touched  = Array{Bool,1}(undef,(4,))
   touched .= true
   by=ArrayBlock(y,touched)
   cache=return_cache(m,by)
   @test vcat(y...)==evaluate!(cache,m,by)

   # MatrixBlock{Matrix} test
   x1=rand(2,3)
   x2=rand(3,3)
   x3=rand(4,3)
   y=Array{Matrix{Float64}}(undef,(3,3))
   y[1,1]=x
   y[2,1]=x2
   y[3,1]=x3
   y[1,2]=x.+3
   y[1,3]=x.+5
   touched  = Array{Bool,2}(undef,(3,3))
   touched .= true
   touched[2:3,2:3].=false
   by=ArrayBlock(y,touched)
   cache=return_cache(m,ArrayBlock(y,touched))
   @test vcat(hcat(y[1,1:3]...),
              hcat(y[2,1],zeros(3,3),zeros(3,3)),
              hcat(y[3,1],zeros(4,3),zeros(4,3)))==evaluate!(cache,m,by)

   # Nested blocking test
   touched_parent      = Array{Bool,2}(undef,(3,3))
   touched_parent     .= false
   touched_parent[2,3] = true
   y_parent            = Array{ArrayBlock{Matrix{Float64}},2}(undef,(3,3))
   y_parent[2,3]       = by
   parent              = ArrayBlock(y_parent,touched_parent)
   cache=return_cache(m,parent)

   result_touched  = Array{Bool,2}(undef,(3,3))
   result_array    = Array{Matrix{Float64}}(undef,(3,3))
   result_touched .= touched_parent
   result_array[2,3] = vcat(hcat(y[1,1:3]...),
                            hcat(y[2,1],zeros(3,3),zeros(3,3)),
                            hcat(y[3,1],zeros(4,3),zeros(4,3)))

   @test result_array[2,3] == evaluate!(cache,m,parent).array[2,3]
end
