export AbstractFEBasis
export FEBasis
export TestFEBasis, TrialFEBasis
export ∇, GradFEBasis
export assemblesystem

# To be eliminated
# export NewTrialFEBasis
# export allocatecellmatrix
# export gettestbasis
# export gettrialbasis
# export getfesp
# export RequiredOpDict
# export requireddiffops!



struct RefFEValues
	id::Array{Float64}
	grad::Array{Float64}
end

mutable struct RequiredOpDict
	reqop::Dict
	function RequiredOpDict()
		opdict = new()
		opdict.reqop=Dict("Id"=>false, "Grad"=>false, "Curl"=>false, "Div"=>false)
		return opdict
	end
end

abstract type AbstractFEBasis end

mutable struct FEBasis <: AbstractFEBasis
	fespace::FESpace
	values::RefFEValues
	istest::Bool
	function FEBasis(fesp::FESpace, istest::Bool)
		febasis=new()
		febasis.fespace=fesp
		febasis.istest=istest
		return febasis
	end
end
TestFEBasis(fesp::FESpace) = FEBasis(fesp::FESpace, true)
TrialFEBasis(fesp::FESpace) = FEBasis(fesp::FESpace, false)

mutable struct GradFEBasis <: AbstractFEBasis
	febasis::FEBasis
end

abstract type AbstractBilinearForm end

mutable struct BilinearForm <: AbstractBilinearForm
    testbasis::AbstractFEBasis
	trialbasis::AbstractFEBasis
end

mutable struct SumBilinearForm <: AbstractBilinearForm
	cellm1::AbstractBilinearForm
	cellm2::AbstractBilinearForm
end

# FUNCTIONS
(∇)(febasis::FEBasis) = GradFEBasis(febasis)

import Base.*
function (*)(febasis1::AbstractFEBasis, febasis2::AbstractFEBasis)
	cfebasis1 = getfebasis(febasis1)
	cfebasis2 = getfebasis(febasis2)
	if ( cfebasis1.istest == cfebasis2.istest ) return 0 end
	( cfebasis1.istest ) ? bilf = BilinearForm(febasis1, febasis2) : bilf = BilinearForm(febasis2, febasis1)
	return bilf
end

import Base.+
function (+)(cellm1::AbstractBilinearForm, cellm2::AbstractBilinearForm)
	return SumBilinearForm(cellm1,cellm2)
end

getfesp(febasis::FEBasis) = febasis.fespace
getfesp(febasis::GradFEBasis) = febasis.febasis.fespace

getfebasis(febasis::FEBasis) = febasis
getfebasis(febasis::GradFEBasis) = febasis.febasis

getvalues(febasis::FEBasis) = febasis.values.id
getvalues(febasis::GradFEBasis) = febasis.febasis.values.grad

getrank(febasis::FEBasis) = febasis.fespace.reffe.rank
getrank(febasis::GradFEBasis) = getrank(febasis.febasis)+1

getspdim(febasis::FEBasis) = dim(febasis.fespace.reffe.polytope)
getspdim(febasis::GradFEBasis) = getspdim(febasis.febasis)

getnumshf(febasis::FEBasis) = febasis.fespace.reffe.numdof
getnumshf(febasis::GradFEBasis) = getnumshf(febasis.febasis)

gettestbasis(cellm::BilinearForm) = getfebasis(cellm.testbasis)
gettrialbasis(cellm::BilinearForm) = getfebasis(cellm.trialbasis)
gettestbasis(cellm::SumBilinearForm) = gettestbasis(cellm.cellm1)
gettrialbasis(cellm::SumBilinearForm) = gettrialbasis(cellm.cellm1)

function requireddiffops!(febasis::FEBasis, dict::RequiredOpDict)
	dict.reqop["Id"]=true
end

function requireddiffops!(febasis::GradFEBasis, dict::RequiredOpDict)
	dict.reqop["Grad"]=true
end

function requireddiffops!(feform::BilinearForm, reqoptest::RequiredOpDict, reqoptrial::RequiredOpDict)
	requireddiffops!(feform.testbasis, reqoptest)
	requireddiffops!(feform.trialbasis, reqoptrial)
	nothing
end

function requireddiffops!(feform::SumBilinearForm, reqoptest::RequiredOpDict, reqoptrial::RequiredOpDict)
	requireddiffops!(feform.cellm1, reqoptest, reqoptrial)
	requireddiffops!(feform.cellm2, reqoptest, reqoptrial)
end

function refcelldiffops(reffe::RefFE, quad::Quadrature, dict::RequiredOpDict)
	refval = [0]; refgrad = [0]
	if get(dict.reqop, "Id", false)
		refval = shfsps(reffe,quad.points)
	end
	if get(dict.reqop, "Grad", false)
		refgrad = gradshfsps(reffe,quad.points)
	end
	return RefFEValues(refval, refgrad)
end

function allocatecellmatrix(cellm::BilinearForm)
	size1 = getnumshf(cellm.testbasis)
	array = zeros(Float64, size1, size1)
	return array
end

allocatecellmatrix(cellm::SumBilinearForm) = allocatecellmatrix(cellm.cellm1)

function computecellmatrix!(matrix::Array{Float64,2}, cellm::BilinearForm)
	testvals = getvalues(cellm.testbasis)
	trialvals = getvalues(cellm.trialbasis)
	rank1 = getrank(cellm.testbasis)
	rank2 = getrank(cellm.trialbasis)
	@assert (rank1==rank2)
	spdim = getspdim(cellm.testbasis)
	flattenfieldarray!(testvals, spdim, rank1)
	flattenfieldarray!(trialvals, spdim, rank2)
	contracttesttrialterms(testvals,trialvals,matrix)
	nothing
end

function computecellmatrix!(matrix::Array{Float64,2},cellm::SumBilinearForm)
	computecellmatrix!(matrix,cellm.cellm1)
	computecellmatrix!(matrix,cellm.cellm2)
end

function flattenfieldarray!(array::Array{Float64}, spdims::Int64, rank::Int64)
	ncoef = spdims^rank
	array = reshape(array, size(array,1), size(array,2), ncoef)
	return array
end

function unflattenfieldarray!(array::Array{Float64,N}, spdims::Int64, rank::Int64) where N
	aux1=spdims*ones(Int64,rank)
	array = reshape(array, size(array,1), size(array,2), aux1...)
	return array
end

function contracttesttrialterms(u, v, a)
	# @assert(size(u)==size(v))
	# @assert(size(a,1)==size(u,1))
	# @assert(size(a,2)==size(u,2))
	@inbounds for l=1:size(u, 3)
		for k=1:size(u, 2)
			for i=1:size(u, 1)
				for j=1:size(u, 1)
					a[i, j] += u[i, k, l] * v[j, k, l]
				end
			end
		end
	end
    return a
end

function assemblecellmatrix!(elmat::Array{Float64,2}, sysmat::Array{Float64,2}, l2gid::Array{Int64,1})
	vsysm = view(sysmat,l2gid,l2gid)
	vsysm+=elmat
end

function assemblesystem(feform::AbstractBilinearForm,quad::Quadrature)
	elmat = allocatecellmatrix(feform)
	testbasis = gettestbasis(feform)
	trialbasis = gettrialbasis(feform)
	testfesp = getfesp(testbasis)
	trialfesp = getfesp(trialbasis)
	mesh = testfesp.mesh
	sysmat = zeros(Float64,testfesp.numgdof,trialfesp.numgdof)
	reqoptest = RequiredOpDict(); reqoptrial = RequiredOpDict()
	requireddiffops!(feform, reqoptest, reqoptrial)
	testbasis.values = refcelldiffops(testfesp.reffe, quad, reqoptest)
	trialbasis.values = refcelldiffops(trialfesp.reffe, quad, reqoptrial)
	for icell=1:length(mesh.cellvefs)
		# Update in physical space
		# updatecelldiffops()
		elmat.=0.0
		computecellmatrix!(elmat,feform)
		l2gid = testfesp.l2giddof[icell]
		sysmat[l2gid,l2gid] += elmat
	end
	return sysmat
end
