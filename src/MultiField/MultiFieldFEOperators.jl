
function AffineFEOperator(
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace},
  assem::Assembler,
  terms::FETerm...)

  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  AffineFEOperator(_test,_trial,assem,terms...)
end

function AffineFEOperator(
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace},
  terms::FETerm...)

  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  AffineFEOperator(_test,_trial,terms...)
end

function AffineFEOperator(
  mat::Type{<:AbstractSparseMatrix},
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace},
  terms::FETerm...)

  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  AffineFEOperator(mat,_test,_trial,terms...)
end

function FEOperator(
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace},
  assem::Assembler,
  terms::FETerm...)

  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  FEOperator(_test,_trial,assem,terms...)
end

function FEOperator(
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace},
  terms::FETerm...)

  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  FEOperator(_test,_trial,terms...)
end

function FEOperator(
  mat::Type{<:AbstractSparseMatrix},
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace},
  terms::FETerm...)

  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  FEOperator(mat,_test,_trial,terms...)
end


