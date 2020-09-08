for op in (:+,:-,:tr, :transpose, :adjoint, :symmetric_part)
  @eval begin

    function ($op)(f::Field)
      operate_fields($op,f)
    end

    #function ($op)(f::AbstractArray{<:Field})
    #  operate_arrays_of_fields($op,f)
    #end

  end
end

for op in (:+,:-,:*,:inner,:outer,:dot)
  @eval begin

    function ($op)(f::Field,g::Field)
      operate_fields($op,f,g)
    end

    #function ($op)(f::AbstractArray{<:Field},g::AbstractArray{<:Field})
    #  operate_arrays_of_fields($op,f,g)
    #end

  end
end

function operate_fields(op::Function,args...)
  k = FieldOpKernel(op)
  apply_kernel_to_field(k,args...)
end
