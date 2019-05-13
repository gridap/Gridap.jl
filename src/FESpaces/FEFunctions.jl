"""
FE Function
"""
abstract type FEFunction{D,Z,T,E} end

free_dofs(::FEFunction) = @abstractmethod
fixed_dofs(::FEFunction) = @abstractmethod
# @santiagobadia : Do we want to make a difference between Dirichlet and constrained ?

struct ConformingFEFunction{D,Z,T,E} <: FEFunction{D,Z,T,E}
	fesp::FESpace{D,Z,T,E}
	dof_values::CellFieldFromExpand{D,T,E,T}
end

free_dofs(this::ConformingFEFunction) = this.dof_values.coeffs.gid_to_val_pos
fixed_dofs(this::ConformingFEFunction) = this.dof_values.coeffs.gid_to_val_neg
