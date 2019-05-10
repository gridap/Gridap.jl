abstract type MeshConformity end
struct ConformingMesh <: MeshConformity end
struct NonConformingMesh <: MeshConformity end

abstract type str{D,Z,T,E} end
MeshConformity(::FESpace)::MeshConformity = error("Not defined")

struct ConformingFESpace{D,Z,T,E} <: FESpace{D,Z,T,E}
	a::D,Z,T,E
end
MeshConformity(::Type{ConformingFESpace{D,Z,T,E}}) where {D,Z,T,E}= ConformingMesh

struct NonConformingFESpace{D,Z,T,E} <: FESpace{D,Z,T,E}
	a::D,Z,T,E
end
MeshConformity(::Type{NonConformingFESpace{D,Z,T,E}}) where {D,Z,T,E} = NonConformingMesh

struct FESpaceWithDirichletData{D,Z,T,E,S<:FESpace{D,Z,T,E}}
	a::D,Z,T,E
	s::S
end
MeshConformity(::Type{FESpaceWithDirichletData{D,Z,T,E,S}}) where {D,Z,T,E,S<:ConformingFESpace{D,Z,T,E}} = ConformingMesh()
MeshConformity(::Type{FESpaceWithDirichletData{D,Z,T,E,S}}) where {D,Z,T,E,S<:NonConformingFESpace{D,Z,T,E}} = NonConformingMesh()

s1 = ConformingFESpace(10)
s2 = NonConformingFESpace(10)
s1c = FESpaceWithDirichletData(20,s1)
s2c = FESpaceWithDirichletData(20,s2)



MeshConformity(typeof(s1))
MeshConformity(typeof(s2))
MeshConformity(typeof(s1c))
MeshConformity(typeof(s2c))
