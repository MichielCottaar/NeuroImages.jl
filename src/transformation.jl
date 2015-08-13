using AffineTransforms
import Base.inv

type ReduceDimensions{T, N, M}
    scalefwd :: Matrix{T}
    offset :: Vector{T}
end

ReduceDimensions{T, N}(aff :: AffineTransform{T, N}) = ReduceDimensions{T, N, M}(aff.tscalefwd, aff.offset)

function ReduceDimensions(index :: Integer, dimension :: Integer, ndims :: Integer)
    scalefwd = zeros(Int, ndims, ndims-1)
    for ix in 1:ndims-1
        scalefwd[ix + (ix >= dimension), ix] = 1
    end
    offset = zeros(Int, ndims)
    offset[dimension] = index
    ReduceDimensions{Int, ndims, ndims-1}(scalefwd, offset)
end

ndims_in{T, N, M}(::ReduceDimensions{T, N, M}) = M
ndims_out{T, N, M}(::ReduceDimensions{T, N, M}) = N
ndims_in{T, N}(::AffineTransform{T, N}) = N
ndims_out{T, N}(::AffineTransform{T, N}) = N

*{T, S, N, M}(rd :: ReduceDimensions{T, N, M}, aff :: AffineTransform{S, M}) = ReduceDimensions{typeof(one(T) + one(S)), N, M}(rd.scalefwd * aff.scalefwd, rd.offset + rd.scalefwd * aff.offset)
*{T, S, N, M}(aff :: AffineTransform{T, N}, rd :: ReduceDimensions{S, N, M}) = ReduceDimensions{typeof(one(T) + one(S)), N, M}(aff.scalefwd * rd.scalefwd, aff.offset + aff.scalefwd * rd.offset)
*{T, S, N, M, L}(rd_f :: ReduceDimensions{T, N, M}, rd_s :: ReduceDimensions{S, M, L}) = ReduceDimensions{typeof(one(T) + one(S)), N, L}(rd_f.scalefwd * rd_s.scalefwd, rd_f.offset + rd_f.scalefwd * rd_s.offset)
*{T, N, M}(rd :: ReduceDimensions{T, N, M}, i :: Number) = rd * tformscale(i, M)
*{T, N, M}(i :: Number, rd :: ReduceDimensions{T, N, M}) = tformscale(i, N) * rd
*(rd :: ReduceDimensions, x :: AbstractVector) = rd.scalefwd * x + rd.offset

inv(aff :: AffineTransform) = AffineTransform(inv(aff.scalefwd), inv(aff.scalefwd) * -aff.offset)
*{T, N}(aff :: AffineTransform{T, N}, i :: Number) = aff * AffineTransforms.tformscale(i, N)
*{T, N}(i :: Number, aff :: AffineTransform{T, N}) = tformscale(i, N) * aff

type ComposeTransform
    first
    second
end

inv(comp :: ComposeTransform) = inv(comp.second) * inv(con.first)
*(comp :: ComposeTransform, x :: AbstractVector) = comp.first * (comp.second * x)

typealias AnyTransform Union(AffineTransform, Number, ComposeTransform)

# composing transforms
*(comp :: ComposeTransform, other :: AffineTransform) = ComposeTransform(comp.first, (comp.second * other))
*(other :: AnyTransform, comp :: ComposeTransform) = ComposeTransform((other * comp.first), comp.second)
*(a :: AnyTransform, b :: AnyTransform) = ComposeTransform(a, b)

index_to_transform(ndims, dimension, index :: Integer) = index_to_transform(ndims, dimension, index:index)
function index_to_transform(ndims, dimension, index :: Range)
    offset = zeros(Float64, ndims)
    offset[dimension] = start(index) - step(index)
    aff = eye(Float64, ndims)
    aff[dimension, dimension] =  step(index)
    AffineTransform(aff, offset)
end

function inds_to_transform(ndims, inds...)
    trans = tformeye(ndims)
    for dim in 1:ndims
        if dim <= size(inds, 1)
            trans = trans * index_to_transform(ndims_in(trans), dim + ndims_in(trans) - ndims, inds[dim])
        end
    end
    trans
end
