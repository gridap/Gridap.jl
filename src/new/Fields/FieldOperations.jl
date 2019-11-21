
# Unary operations on fields and arrays of fields

for op in (:+,:-,:tr, :transpose, :adjoint, :symmetic_part)
  @eval begin

    function ($op)(f::Field)
      apply_kernel_to_field(bcast($op),f)
    end

    function ($op)(f::AbstractArray{<:Field})
      apply_to_field_array(bcast($op),f)
    end

  end

end

# Binary operations on fields and arrays of fields

for op in (:+,:-,:*,:inner,:outer)
  @eval begin

    function ($op)(f::Field, g::Field)
      apply_kernel_to_field(bcast($op),f,g)
    end

    function ($op)(f::AbstractArray{<:Field},g::AbstractArray{<:Field})
      apply_to_field_array(bcast($op),f,g)
    end

  end

end

# Define the gradient of some operations

@inline apply_kernel_gradient(k::BCasted{typeof(+)},a) = field_gradient(a)

@inline apply_kernel_gradient(k::BCasted{typeof(-)},a) = apply_kernel_to_field(k,field_gradient(a))

for op in (:+,:-)
  @eval begin
    @inline function apply_kernel_gradient(k::BCasted{typeof($op)},a,b)
      ga = field_gradient(a)
      gb = field_gradient(b)
      apply_kernel_to_field(k,ga,gb)
    end
    @inline function apply_kernel_gradient(k::BCasted{typeof($op)},f...)
      apply_kernel_to_field(k,field_gradients(f...)...)
    end
  end
end

for op in (:+,:-)
  @eval begin
    function apply_gradient(k::Valued{BCasted{typeof($op)}},f...)
      g = field_array_gradients(f...)
      apply(k,g...)
    end
  end
end

