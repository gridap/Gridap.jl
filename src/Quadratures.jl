using QuadGK

export Quadrature
"""
Tensor product quadrature rule (nodes and weights) integrating exactly 2´order´-1 polynomials
"""
struct Quadrature
	points::Array{Array{Float64,1},1}
	weights::Array{Array{Float64,1},1}
	tppoints::Array{Float64,2}
	tpweights::Array{Float64,1}
	function Quadrature(gps::Array{Int64,1})
		spdims = length(gps)
		points=Array{Array{Float64,1}}(undef,length(gps))
		weights=Array{Array{Float64,1}}(undef,length(gps))
		points = [gauss(gps[i]) for i=1:length(gps)]
		a = [points[i][1] for i=1:length(gps)]
		b = [points[i][2] for i=1:length(gps)]
		tpa = Array{Float64,spdims+1}(undef,tuple(gps...,spdims))
		tensorfill!(tpa,a)
		tpa = reshape(tpa,:,spdims)
		tpb = Array{Float64,spdims+1}(undef,tuple(gps...,spdims))
		tensorfill!(tpb,b)
		tpb = reshape(tpb,:,spdims)
		tpb = [prod(tpb[i,:]) for i=1:size(tpb,1)]
		return new(a,b,tpa,tpb)
	end
end

@generated function tensorfill!(A::Array{T,N},c) where {T,N}
    quote
        @nloops $N i A begin
            t = @ntuple $N i; d = length(t)
            dim = t[d]
            nodedim = t[dim]
            (@nref $N A i) = c[dim][nodedim]
        end
    end
end
